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
    
%     spin_collection=SpinCollection( FromFile([INPUT_FILE_PATH, obj.parameters.InputFile]));
    dist=p.get_parameter('SetBathSpins','SquareLength');% the distance between the two closest 2Fe2S clusters
    obj.parameters.SquareLength=dist;
    seed=p.get_parameter('SetBathSpins','Seed');%set the rand seed
    obj.parameters.seed=seed;
    nseed=p.get_parameter('SetBathSpins','NSeed');% the number of the copies of the MagR with rand displacement
    obj.parameters.nseed=nseed;
    rRand=p.get_parameter('SetBathSpins','rRandom');% the displacement of the MagR is rRand*vec, where vec is a normalized randdom vector
    obj.parameters.rRand=rRand;
    
    SpinListTotal=get_spin_list(dist,seed,nseed,rRand);
    spin_collection_tot=SpinCollection(FromSpinList(SpinListTotal));
    % reset the gyromagnetic ratio
    obj.parameters.reset_gyromagnetic_ratio=p.get_parameter('SetBathSpins','ResetGyroMagneticRatio');
    if obj.parameters.reset_gyromagnetic_ratio
        gamma2reset=p.get_parameter('SetBathSpins','Gamma2Reset');
        obj.parameters.gamma2reset=gamma2reset;
        nspin=length(SpinListTotal);
        for kk=1:nspin
            SpinListTotal{kk}.gamma=gamma2reset;
        end        
    end

   if obj.parameters.SetBathSpins.SetSpin;
       paraCell=obj.parameters.SetBathSpins.BathSpinsSettingCell;
       spin_collection_tot.set_spin(paraCell);
   else
       spin_collection_tot.set_spin();
   end
   set_spin_qAxis(obj,spin_collection_tot);
   obj.keyVariables('SpinListTotal')=SpinListTotal;
end

function spin_list=get_spin_list(dist,seed,nseed,rRand)
    coord_cell=generate_MagR_coordinate(dist,seed,nseed,rRand);
    nspin=length(coord_cell);
    spin_list=cell(1,nspin);
    for n=1:nspin
        name='Fe2S2';
        coord=coord_cell{n};
        spin_list{n}=model.phy.PhysicalObject.Spin(name,coord);
    end
end

function coord_cell=generate_MagR_coordinate(dist,seed,nseed,rRand)
    nspin=4*nseed;
    coord_square=[-0.5*dist,0.5*dist,0;...
                  0.5*dist,0.5*dist,0;...
                  0.5*dist,-0.5*dist,0;...
                  -0.5*dist,-0.5*dist,0];
    rng(seed);
    rand_radius=rRand*rand(nseed,1);
%     rand_radius=rRand*(0:0.5:2);
    rotation_angles=pi/4*rand(nseed,1);
%     rotation_angles=(0:1:45)';
    rand_vec_list=randn(nseed,3);
%     philist=-pi:0.1*pi:pi;
%     rand_vec_list=[cos(philist'),sin(philist'),zeros(nseed,1)];
    coord_cell=cell(1,nspin);
    pos=1;
    for kk=1:nseed
        vec=rand_vec_list(kk,:);
        vec=rand_radius(kk)*vec/norm(vec);
        R=rotz(rotation_angles(kk));
        rotated_square=R*(coord_square');
        rotated_square=rotated_square';
        coord_cell{pos}=rotated_square(1,:)+vec;pos=pos+1;
        coord_cell{pos}=rotated_square(2,:)+vec;pos=pos+1;
        coord_cell{pos}=rotated_square(3,:)+vec;pos=pos+1;
        coord_cell{pos}=rotated_square(4,:)+vec;pos=pos+1;
    end
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