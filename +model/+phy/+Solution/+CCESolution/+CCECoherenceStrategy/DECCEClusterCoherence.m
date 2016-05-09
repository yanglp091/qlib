classdef DECCEClusterCoherence < model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence
    %ECCECLUSTERCOHERENCE 
                %calculate the coherence of single cluster for ensemble CCE 
    properties
      coherence_tilde
      vertical_decay_rates
      parallel_decay_rates
    end
    
    methods
        function obj=DECCEClusterCoherence(cluster_spin_index,cluster_parameters)
            obj@model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence();
            if nargin >0
               obj.generate(cluster_spin_index,cluster_parameters); 
            end                    
        end
        
        function coh=calculate_cluster_coherence(obj,evolution_para)
            obj.npulse=evolution_para.npulse;
            center_spin_states=evolution_para.center_spin_states;
            is_secular=evolution_para.is_secular;
            obj.timelist=evolution_para.timelist;
            
            obj.vertical_decay_rates=evolution_para.vertical_decay_rates;
            obj.parallel_decay_rates=evolution_para.parallel_decay_rates;
            
            %generate the spin_collection for this cluster including the central spin
            obj.spin_collection= model.phy.SpinCollection.SpinCollection();
            obj.spin_collection.spin_source=model.phy.SpinCollection.Strategy.FromSpinList([{obj.center_spin},obj.cluster_bath_spin]);
            obj.spin_collection.generate();
             
            hamiCell=obj.gen_reduced_hamiltonian(center_spin_states,is_secular);
            [h_list,hami_prefactor]=obj.gen_hami_list(hamiCell);
            
            [bath_cluster_sc,denseMat,initial_state_type]=obj.set_initial_state;
            
            % get evolution kernal
            [liouList,prefactors]=obj.gen_liouvillian_list(bath_cluster_sc,h_list,hami_prefactor);
            
            %Observable
            obs=model.phy.QuantumOperator.SpinOperator.Observable(bath_cluster_sc,'IdentityMatrix');
            dim=obs.dim;
            obs.setMatrix(speye(dim));
              
            coh=obj.calculate_coherence_liouville(liouList,prefactors,obs,denseMat,initial_state_type);
        end
        function [bath_cluster_sc,denseMat,initial_state_type]=set_initial_state(obj)
            %generate a SpinCollection of bath spins in this cluster
            bath_cluster_sc= model.phy.SpinCollection.SpinCollection();
            bath_cluster_sc.spin_source=model.phy.SpinCollection.Strategy.FromSpinList(obj.cluster_bath_spin);
            bath_cluster_sc.generate();            
            
            % DensityMatrix
            denseMat=model.phy.QuantumOperator.SpinOperator.DensityMatrix(bath_cluster_sc,'IdentityMatrix');
            dim=denseMat.dim;
            denseMat.setMatrix(eye(dim)/dim);
                      
            initial_state_type='MixedState';
        end
        
        function [liouList,prefactor]=gen_liouvillian_list(obj,bath_cluster,hami_list,hami_prefactor)
            % Here, we want to add decoherence operator for bath spins
            % The Liouvillian is generated as the order: first get the
            % direct product of the operator for different spins, then use
            % the rule: A*rho -->E\otimes A.';rho*B -->B\otimes E.';A*rho*B -->B\otimes A.';
            % to transform a operator (A,B, ...) in the Hilbert space to Liouville space;
            % Transformation of density matrix to a vector in Liouville space must be take the following method:
            % first get the matrix of the density matrix rho, the vec_rho=rho(:);
            import model.phy.QuantumOperator.SpinOperator.DecoherenceSuperOperator.CNMDecohSuperOperator

            %setting the parameters for decoherence supper operator
            nspin=bath_cluster.getLength;
            para.AddVerticalDecay=0;
            para.AddParallelDecay=0;
            vertical_decay_rate=obj.vertical_decay_rates;            
            if vertical_decay_rate>0
                para.AddVerticalDecay=1;
                para.VerticalDecayRateList=vertical_decay_rate*ones(1,nspin);
            end            
            parallel_decay_rate=obj.parallel_decay_rates;
            if parallel_decay_rate>0
                para.AddParallelDecay=1;
                para.ParallelDecayRateList=parallel_decay_rate*ones(1,nspin);
            end

            L_decay=CNMDecohSuperOperator(bath_cluster,para);
            L_decay_mat=L_decay.getMatrix;
            noperator=length(hami_list);
            if mod(noperator,2)==1
                error('The number of the cce hamiltonians is not a even number.');
            end
            liouList=cell(1,noperator/2);
            for kk=1:noperator/2
                hami1=hami_list{kk};
                hami2=hami_list{noperator-kk+1};   
                
                Amat=hami1.getMatrix(); Bmat=hami2.getMatrix(); eyeMat=speye(hami1.dim);
                L_kk=model.phy.QuantumOperator.MultiSpinSuperOperator(bath_cluster);
                Lmat=kron(eyeMat, Amat)-kron(Bmat.', eyeMat)+1i*L_decay_mat;
                L_kk.setMatrix(Lmat);
                liouList{noperator/2-kk+1}=L_kk;
            end
            
            prefactor=-1i*abs(hami_prefactor(1:noperator/2));           
        end                       
    end
    
end
