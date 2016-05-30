function BathSpinParameters(obj, p)
     %set bath spin
     % get the parameters for setting the bath spins e.g. ZFS, eta, principle_axis, coordinate
    import model.phy.SpinCollection.Strategy.FromFile
    import model.phy.SpinCollection.Strategy.FromSpinList
    import model.phy.SpinCollection.SpinCollection
    
    obj.parameters.SetBathSpins.SetSpin=p.get_parameter('SetBathSpins','SetSpin');
    if obj.parameters.SetBathSpins.SetSpin
        SpeciesNumber=p.get_parameter('SetBathSpins','SpeciesNumber');
        BathSpinsSettingCell=cell(1,4);
        NameList=p.get_parameter('SetBathSpins','Name');
        ZFSList=p.get_parameter('SetBathSpins','ZFS');
        etaList=p.get_parameter('SetBathSpins','eta');
        paxisList=p.get_parameter('SetBathSpins','principle_axis');
        for k=1:SpeciesNumber
            BathSpinsSettingCell{k}.name=NameList(k,:);
            BathSpinsSettingCell{k}.ZFS=ZFSList(k,:);
            if ~strcmp(etaList,'default')                   
                BathSpinsSettingCell{k}.eta=etaList(k,:);
            end
            if ~strcmp(paxisList,'default')                   
                BathSpinsSettingCell{k}.principle_axis=paxisList(k,:);
            end                  
        end
        obj.parameters.SetBathSpins.BathSpinsSettingCell=BathSpinsSettingCell;               
    end
    
    spin_collection=SpinCollection( FromFile([INPUT_FILE_PATH, obj.parameters.InputFile]));
    % reset the gyromagnetic ratio
    obj.parameters.reset_gyromagnetic_ratio=p.get_parameter('SetBathSpins','ResetGyroMagneticRatio');
    if obj.parameters.reset_gyromagnetic_ratio
        gamma2reset=p.get_parameter('SetBathSpins','Gamma2Reset');
        obj.parameters.gamma2reset=gamma2reset;
        nspin=spin_collection.getLength;
        for kk=1:nspin
            spin_collection.spin_list{kk}.gamma=gamma2reset;
        end        
    end

   if obj.parameters.SetBathSpins.SetSpin;
       paraCell=obj.parameters.SetBathSpins.BathSpinsSettingCell;
       spin_collection.set_spin(paraCell);
   else
       spin_collection.set_spin();
   end
   set_spin_qAxis(obj,spin_collection);
   obj.keyVariables('spin_collection')=spin_collection;
end

function set_spin_qAxis(obj,spin_collection)
    magnetic_field=obj.parameters.MagneticField;
    mf_len=norm(magnetic_field);
    mf_z=[0,0,1]*magnetic_field';
    if mf_len>mf_z
        vec_z=magnetic_field/mf_len;
        vec_y=cross([0,0,1],vec_z);vec_y=vec_y/norm(vec_y);
        vec_x=cross(vec_y,vec_z);vec_x=vec_x/norm(vec_x);
        Axis=[vec_x;vec_y;vec_z];
    else
        Axis=eye(3);
    end
    
    nspin=spin_collection.getLength;
    for kk=1:nspin
        spin=spin_collection.spin_list{kk};
        spin.qAxis=Axis;
    end
end