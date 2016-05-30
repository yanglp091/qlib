function CentralSpinParameters(obj,p)
     %set central spin
     % get the parameters for setting the central spin e.g. ZFS, eta, principle_axis, coordinate
     
    obj.parameters.SetCentralSpin.name=p.get_parameter('SetCentralSpin', 'Name');
    obj.parameters.SetCentralSpin.CentralSpinStates=p.get_parameter('SetCentralSpin', 'CentralSpinStates');
    obj.parameters.SetCentralSpin.SetSpin=p.get_parameter('SetCentralSpin', 'SetSpin');
    Orientation=p.get_parameter('SetCentralSpin', 'Orientation');
    Isotope=p.get_parameter('SetCentralSpin', 'Isotope');
    Coordinate=p.get_parameter('SetCentralSpin', 'Coordinate');
    if obj.parameters.SetCentralSpin.SetSpin
        if ~strcmp(Orientation,'default') 
            obj.parameters.SetCentralSpin.orientation=Orientation;
        end
        if ~strcmp(Isotope,'default') 
            obj.parameters.SetCentralSpin.isotope=Isotope;
        end
        if ~strcmp(Coordinate,'default') 
            obj.parameters.SetCentralSpin.coordinate=Coordinate;
        end
    end

%    center_spin_name=obj.parameters.SetCentralSpin.name;
%    para_central_spin=obj.parameters.SetCentralSpin; 
%    center_spin=eval(strcat(center_spin_name,'(','para_central_spin',')'));
%    obj.keyVariables('center_spin')=center_spin;
end