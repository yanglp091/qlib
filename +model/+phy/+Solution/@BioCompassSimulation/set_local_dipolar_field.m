function  set_local_dipolar_field(obj,spin_collection,nv_state)
    import model.phy.SpinCollection.Strategy.FromSpinList
    import model.phy.SpinCollection.SpinCollection
    
    nMagRSpin=obj.parameters.MagRSpinNumber;
    if nMagRSpin<1
       return; 
    end
    
    nv1=obj.parameters.NVCenter1.espin;
    nv2=obj.parameters.NVCenter2.espin;
    st1_idx=nv_state(1);st1=nv1.eigen_vect(:,st1_idx);
    st2_idx=nv_state(2);st2=nv2.eigen_vect(:,st2_idx);
    sz1=st1'*Sz(3)*st1;
    sz2=st2'*Sz(3)*st2;
    svec1=[0,0,sz1];
    svec2=[0,0,sz2];
    
    for kk=1:nMagRSpin
        spin=spin_collection.spin_list{kk};
        sc1=SpinCollection( FromSpinList([INPUT_FILE_PATH, [{nv1},{spin}]]));
        sc2=SpinCollection( FromSpinList([INPUT_FILE_PATH, [{nv2},{spin}]]));
        dip1=model.phy.SpinInteraction.DipolarInteraction(sc1);
        dip2=model.phy.SpinInteraction.DipolarInteraction(sc2);
        coeff_mat1=dip1.calculate_coeff([{nv1},{spin}]);
        coeff_mat2=dip2.calculate_coeff([{nv2},{spin}]);
        local_field=-svec1*coeff_mat1/spin.gamma-svec2*coeff_mat2/spin.gamma;%-s_val*s_vec*coeff/spin2.gamma
        spin.local_field=local_field;
    end
end