classdef ECCEClusterCoherence < model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence
    %ECCECLUSTERCOHERENCE 
                %calculate the coherence of single cluster for ensemble CCE 
    properties
      coherence_tilde
    end
    
    methods
        function obj=ECCEClusterCoherence(cluster_spin_index,cluster_parameters)
            obj@model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence();
            if nargin >0
               obj.generate(cluster_spin_index,cluster_parameters); 
            end                    
        end
        
        function coh=calculate_cluster_coherence(obj,evolution_para,varargin)
            obj.npulse=evolution_para.npulse;
            center_spin_states=evolution_para.center_spin_states;
            is_secular=evolution_para.is_secular;
            timelist=evolution_para.timelist;
            
            %generate the spin_collection for this cluster including the central spin
            obj.spin_collection= model.phy.SpinCollection.SpinCollection();
            obj.spin_collection.spin_source=model.phy.SpinCollection.Strategy.FromSpinList([{obj.center_spin},obj.cluster_bath_spin]);
            obj.spin_collection.generate();
             
            hamiCell=obj.gen_reduced_hamiltonian(center_spin_states,is_secular);
            [~]=obj.gen_hami_list(hamiCell);
            
            
            %generate a SpinCollection of bath spins in this cluster
            bath_cluster_sc= model.phy.SpinCollection.SpinCollection();
            bath_cluster_sc.spin_source=model.phy.SpinCollection.Strategy.FromSpinList(obj.cluster_bath_spin);
            bath_cluster_sc.generate();
            
            % DensityMatrix
            denseMat=model.phy.QuantumOperator.SpinOperator.DensityMatrix(bath_cluster_sc,'IdentityMatrix');
            dim=denseMat.dim;
            denseMat.setMatrix(eye(dim)/dim);
            

            if dim<2000
                coh=obj.calculate_coherence_hilbert(bath_cluster_sc,denseMat,'MixedState',timelist);
            else
                coh=obj.calculate_coherence_liouville(bath_cluster_sc,denseMat,'MixedState',timelist);
            end
            obj.coherence=coh;
        end
                       
    end
    
end
