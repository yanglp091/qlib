function obs=GetObservables(obj,spin_collection,singlet_state)
    import model.phy.QuantumOperator.SpinOperator.Observable
    MagR_sc=obj.parameters.MagRSpinCollection;
    dim_magR=MagR_sc.getDim;
    
    obs=Observable(spin_collection,'SingletStateProjector');
    obs.setName('singletStateProbability');
    obs.setMatrix( kron(speye(dim_magR), singlet_state*singlet_state') );
end