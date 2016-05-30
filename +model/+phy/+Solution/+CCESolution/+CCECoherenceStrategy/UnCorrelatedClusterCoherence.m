classdef UnCorrelatedClusterCoherence < model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence
    %UNCORRELATEDCLUSTERCOHERENCE :: This class is to calculate the
    %coherence of a big cluster. This cluster is a whole and it is
    %uncorrelated with other bath spins.
    
    properties
        vertical_decay_rates
        parallel_decay_rates
    end
    
    methods
        function obj=UnCorrelatedClusterCoherence(cluster_spin_index,cluster_parameters)
            obj@model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence();
            if nargin >0
               obj.generate(cluster_spin_index,cluster_parameters); 
            end                    
        end 
        
       function coh=calculate_cluster_coherence(obj,evolution_para,varargin)
            obj.npulse=evolution_para.npulse;
            center_spin_states=evolution_para.center_spin_states;
            is_secular=evolution_para.is_secular;
            obj.timelist=evolution_para.timelist;
            
            IntPara.AddContact=evolution_para.AddContact;
            IntPara.contactInt=evolution_para.contactInt;
            IntPara.AddDipInt=evolution_para.AddDipInt;
                        
            obj.spin_collection=obj.bath_spin_collection;
             
            hamiCell=obj.gen_reduced_hamiltonian(center_spin_states,is_secular,IntPara);
            [h_list,hami_prefactor]=obj.gen_hami_list(hamiCell);            
            [denseMat,initial_state_type]=obj.set_initial_state;
            
            AddDecay=evolution_para.AddDecay;
            if AddDecay
                obj.vertical_decay_rates=evolution_para.vertical_decay_rates;
                obj.parallel_decay_rates=evolution_para.parallel_decay_rates;
                obs=model.phy.QuantumOperator.SpinOperator.Observable(obj.spin_collection,'IdentityMatrix');
                dim=obs.dim;
                obs.setMatrix(speye(dim));
                [liouList,prefactors]=obj.gen_liouvillian_list(h_list,hami_prefactor);
                coh=obj.calculate_coherence_liouville(liouList,prefactors,obs,denseMat,initial_state_type);
            else
                obs=model.phy.QuantumOperator.SpinOperator.Observable(obj.spin_collection,'IdentityMatrix');
                obs.setMatrix(1);
                coh=obj.calculate_coherence_hilbert(h_list,hami_prefactor,obs,denseMat,initial_state_type);
            end
       end
       
       function reduced_hami = gen_reduced_hamiltonian(obj,center_spin_states,is_secular,IntPara)
            import model.phy.SpinInteraction.SpinChainInteraction.XXInteraction
            import model.phy.SpinInteraction.SpinChainInteraction.YYInteraction
            import model.phy.SpinInteraction.SpinChainInteraction.ZZInteraction

            cluster=obj.spin_collection;

            obj.set_local_field(center_spin_states(1));
            hami1=model.phy.QuantumOperator.SpinOperator.Hamiltonian(cluster);
            zee_interaction=model.phy.SpinInteraction.ZeemanInteraction(cluster);
            hami1.addInteraction(zee_interaction);
            if IntPara.AddDipInt
                dip_interaction=model.phy.SpinInteraction.DipolarInteraction(cluster);
                hami1.addInteraction(dip_interaction);
            end
            if IntPara.AddContact
                para.interaction=IntPara.contactInt;
                hami1.addInteraction( XXInteraction(cluster, para));
                hami1.addInteraction( YYInteraction(cluster, para));
                hami1.addInteraction( ZZInteraction(cluster, para));
            end
            hami1.getMatrix;

            obj.set_local_field(center_spin_states(2));
            hami2=model.phy.QuantumOperator.SpinOperator.Hamiltonian(cluster);                
            zee_interaction=model.phy.SpinInteraction.ZeemanInteraction(cluster);
            hami2.addInteraction(zee_interaction);
            if IntPara.AddDipInt
                dip_interaction=model.phy.SpinInteraction.DipolarInteraction(cluster);
                hami2.addInteraction(dip_interaction);
            end
            if IntPara.AddContact
                para.interaction=IntPara.contactInt;
                hami2.addInteraction( XXInteraction(cluster, para));
                hami2.addInteraction( YYInteraction(cluster, para));
                hami2.addInteraction( ZZInteraction(cluster, para));
            end
            hami2.getMatrix;


            if is_secular && hami1.spin_collection.getLength>1
                approx1=model.phy.SpinApproximation.SpinSecularApproximation(hami1.spin_collection);
                approx2=model.phy.SpinApproximation.SpinSecularApproximation(hami2.spin_collection);
                hami1.apply_approximation(approx1);
                hami2.apply_approximation(approx2);
            end

            reduced_hami=cell(1,2);
            reduced_hami{1,1}=hami1;
            reduced_hami{1,2}=hami2;
       end
       
       function set_local_field(obj,center_spin_state)
            nspin=obj.spin_collection.getLength;% get the number of the total bath spin 
            center_spin=obj.center_spin;
            for kk=1:nspin
               target_spin=obj.spin_collection.spin_list{kk};
               local_field=obj.calculate_dipolar_field(target_spin,center_spin,center_spin_state);% this function is defined in AbstractClusterCoherence              
               target_spin.local_field=local_field;
            end            
       end
        
       function [denseMat,initial_state_type]=set_initial_state(obj)
                      
            % DensityMatrix
            denseMat=model.phy.QuantumOperator.SpinOperator.DensityMatrix(obj.spin_collection,'IdentityMatrix');
            dim=denseMat.dim;
            denseMat.setMatrix(eye(dim)/dim);
                      
            initial_state_type='MixedState';
        end
        
        function [liouList,prefactor]=gen_liouvillian_list(obj,hami_list,hami_prefactor)
            % Here, we want to add decoherence operator for bath spins
            % The Liouvillian is generated as the order: first get the
            % direct product of the operator for different spins, then use
            % the rule: A*rho -->E\otimes A.';rho*B -->B\otimes E.';A*rho*B -->B\otimes A.';
            % to transform a operator (A,B, ...) in the Hilbert space to Liouville space;
            % Transformation of density matrix to a vector in Liouville space must be take the following method:
            % first get the matrix of the density matrix rho, the vec_rho=rho(:);
            import model.phy.QuantumOperator.SpinOperator.DecoherenceSuperOperator.CNMDecohSuperOperator

            %setting the parameters for decoherence supper operator
            bath_cluster=obj.spin_collection;
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

