function SetCentralSpins(obj,p)
%set central spin
% get the parameters for setting the central spin e.g. ZFS, eta, principle_axis, coordinate
    import model.phy.PhysicalObject.NV  
     
    obj.parameters.NVStates=p.get_parameter('SetCentralSpin', 'CentralSpinStates');
    depth=p.get_parameter('SetCentralSpin', 'Depth');obj.parameters.NVDepth=depth;
    dist=p.get_parameter('SetCentralSpin', 'Distance');obj.parameters.NVDistance=dist;

    para1.orientation=[0,0,1];para1.coordinate=[-dist/2,0,-depth];
    para2.orientation=[0,0,1];para2.coordinate=[dist/2,0,-depth];
    
    obj.parameters.NVCenter1=NV(para1);
    obj.parameters.NVCenter2=NV(para2);
 end