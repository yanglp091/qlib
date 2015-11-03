function [eField, hField] = scattedwavefunction( obj, x, y, z,a,b,p,q )
%scattedwavefunction Summary of this function goes here
%   Detailed explanation goes here
% as it's changed from lg1.focBeam.wavefunction,it would be better be an formal
% function.
import model.phy.PhysicalObject.LaserBeam.assist.gammaMN
import model.math.misc.Cart2Sph

eField=zeros(1, 3); hField=zeros(1, 3);
[r, theta, phi]=Cart2Sph(x, y, z);
kr=obj.k*r;
if kr==0
    kr=eps;
end

Z=obj.medium.Z;
for kk=1:length(obj.aNNZ)
    aTerm=obj.aNNZ(kk,:);
    n=aTerm(1); m=aTerm(2); amn=aTerm(3);
    
    prefact = (1.j)^(n-1)* sqrt(4.0*pi);% See Doc: VSWF.m
    [M_mode, N_mode] = ott13.vswfcart(n,m,kr,theta,phi,3);
    eField=eField+ amn* prefact *N_mode;
    hField=hField+ amn* prefact *M_mode / (1.j*Z);
end

for kk=1:length(obj.bNNZ)
    bTerm=obj.bNNZ(kk,:);
    n=bTerm(1);  m=bTerm(2); bmn=bTerm(3);
    prefact = (1.j)^(n-1)* sqrt(4.0*pi);% See Doc: VSWF.m
    [M_mode,N_mode] = ott13.vswfcart(n,m,kr,theta,phi,3);
    eField=eField+ bmn* prefact *M_mode;
    hField=hField+ bmn* prefact *N_mode / (1.j*Z);
end
eField=eField*obj.AmplitudeFactor;
hField=hField*obj.AmplitudeFactor;



end

