function [hami_matix_list,NVDipIntStrength]=GetHamiltonian(obj,spin_collection)
    import model.phy.SpinCollection.Strategy.FromSpinList
    import model.phy.SpinCollection.SpinCollection
    import model.phy.SpinInteraction.SpinChainInteraction.XXInteraction
    import model.phy.SpinInteraction.SpinChainInteraction.YYInteraction
    import model.phy.SpinInteraction.SpinChainInteraction.ZZInteraction
    
%     %add dipolar interaction between the NV and the MagR
%     spin_index_list=cell(2*nMagRSpin,1);
%     for kk=1:nMagRSpin
%        spin_index_list{2*kk-1}=[kk,nMagRSpin+1];
%        spin_index_list{2*kk}=[kk,nMagRSpin+2];
%     end
%     iter=model.phy.SpinCollection.Iterator.SpinIterator(spin_collection,spin_index_list);
%     dip_NV_MagR=model.phy.SpinInteraction.DipolarInteraction(spin_collection,iter);
%     hamiltonian.addInteraction(dip_NV_MagR);


    hami_matix_list=cell(1,4);
    nMagRSpin=obj.parameters.MagRSpinNumber;
 %% Hamiltonian of the interaction between the electron spins in the MagR  
    if  nMagRSpin>1
        hami_inner=model.phy.QuantumOperator.SpinOperator.Hamiltonian(spin_collection);
        AddDipInt=obj.parameters.AddDipInt;
        AddContact=obj.parameters.AddContact;
        % dipoalr interaction between spins in the MagR
        idx=1:nMagRSpin;
        idx_list=nchoosek(idx,2);
        npair=size(idx_list,1);
        spin_index_list=cell(npair,1);
        for kk=1:npair
           spin_index_list{kk}=idx_list(kk,:);
        end

        if AddDipInt
            iter=model.phy.SpinCollection.Iterator.SpinIterator(spin_collection,spin_index_list);
            dip_MagR=model.phy.SpinInteraction.DipolarInteraction(spin_collection,iter);
            hami_inner.addInteraction(dip_MagR);
        end
        % contact interaction between spin in the MagR
        if AddContact    
            IntPara.interaction=obj.parameters.contactInt;
            IntPara.AddIndexList=1;
            IntPara.IndexList=spin_index_list;

            hami_inner.addInteraction( XXInteraction(spin_collection, IntPara));
            hami_inner.addInteraction( YYInteraction(spin_collection, IntPara));
            hami_inner.addInteraction( ZZInteraction(spin_collection, IntPara));
        end

        hm_int=hami_inner.getMatrix;
    else
        hm_int=0;
    end
    
%%  interaction of the Zeeman splitting of the spins in the MagR under the dipolar field of NV centers
    %calcualte the NV-NV dipolar interaction strength
    nv1=obj.parameters.NVCenter1.espin;
    nv2=obj.parameters.NVCenter2.espin;
    sc=SpinCollection( FromSpinList([INPUT_FILE_PATH, [{nv1},{nv2}] ]));
    dip=model.phy.SpinInteraction.DipolarInteraction(sc);
    coef_mat=dip.calculate_coeff([{nv1},{nv2}]);
    NVDipIntStrength=coef_mat(3,3);

    state_list=[1,1;2,2;2,1;1,2];
    for kk=1:4
        nv_state=state_list(kk,:);    
        %add dipolar field provided by NV centers
        obj.set_local_dipolar_field(spin_collection,nv_state);
        hami0_redued=model.phy.QuantumOperator.SpinOperator.Hamiltonian(spin_collection);
        zee_interaction=model.phy.SpinInteraction.ZeemanInteraction(spin_collection);
        hami0_redued.addInteraction(zee_interaction);
        hm0_red=hami0_redued.getMatrix;
        hami_matix_list{kk}=hm_int+hm0_red;
        
    end
    
end