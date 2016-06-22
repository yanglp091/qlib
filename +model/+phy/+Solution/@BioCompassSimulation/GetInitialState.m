function [initial_st,spin_collection]=GetInitialState(obj)
% Here, the initial state of the Radical pair is the singlet state.
% We use the random single sample states to simulation the states of the MagR...
% It should be noted that the NV centers is placed at the end of the spin lists ...
%     import model.phy.SpinCollection.Strategy.FromSpinList
%     import model.phy.SpinCollection.SpinCollection 
% 
%     nv1=obj.parameters.NVCenter1.espin;
%     nv2=obj.parameters.NVCenter2.espin;
%     MagR_spins=obj.parameters.MagRSpinCollection.spin_list;
%     spin_list_tot=[MagR_spins,{nv1},{nv2}];
%     spin_collection=SpinCollection( FromSpinList(spin_list_tot) );
%     obj.keyVariables('spin_collection')=spin_collection;
% 
%     nv_st_idx=obj.parameters.NVStates;
%     nv1_st1=nv1.eigen_vect(:, nv_st_idx(1) );nv1_st2=nv1.eigen_vect(:, nv_st_idx(2) );
%     nv2_st1=nv2.eigen_vect(:, nv_st_idx(1) );nv2_st2=nv2.eigen_vect(:, nv_st_idx(2) );
%     vec1=(nv1_st1+nv1_st2)/sqrt(2);
%     vec2=(nv2_st1+nv2_st2)/sqrt(2);
%     nv_init_st=kron(vec1,vec2)*kron(vec1,vec2)';

    spin_collection=obj.parameters.MagRSpinCollection;
    dim_MagR=spin_collection.getDim;
    initial_st=speye(dim_MagR);
end