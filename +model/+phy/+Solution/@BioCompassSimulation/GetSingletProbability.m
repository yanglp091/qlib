function PS=GetSingletProbability(obj,hamiltonian,initial_state,obs,initial_state_type)
    TimeList=obj.parameters.TimeList;
    ntime=length(TimeList);
    
    switch initial_state_type
        
        case 'PureState'
            hm=full(hamiltonian.getMatrix);
            [V,D]=eig(hm);
            d_vec=diag(D);

            PS=zeros(1,ntime);
            PS(1)=1;
            for kk=2:ntime
                time=TimeList(kk);
                d_vec_k=exp(-1i*time*d_vec);
                state_k=V*diag(d_vec_k)*V'*initial_state;
                PS(kk)=state_k'*obs.getMatrix*state_k;
            end
            
        case 'MixedState'
            
            
    end 
end