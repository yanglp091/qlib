classdef SSCCEClusterCoherence < model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence
    %ECCECLUSTERCOHERENCE 
                %calculate the coherence of single cluster for single sample CCE 
    properties
      bath_spin_state
    end
    
    methods
        function obj=SSCCEClusterCoherence(cluster_spin_index,cluster_parameters)
            obj@model.phy.Solution.CCESolution.CCECoherenceStrategy.AbstractClusterCoherence();
            if nargin >0
               obj.generate(cluster_spin_index,cluster_parameters);
               obj.bath_spin_state=cluster_parameters.bath_spin_state;
            end                    
        end
        
        function coh=calculate_cluster_coherence(obj,center_spin_states,timelist,varargin)
            p = inputParser;
            addRequired(p,'center_spin_states');
            addRequired(p,'timelist');
            addOptional(p,'npulse',0,@isnumeric);
            addOptional(p,'is_secular',0,@isnumeric); 

            parse(p,center_spin_states,timelist,varargin{:});

            obj.npulse=p.Results.npulse;
            center_spin_states=p.Results.center_spin_states;
            is_secular=p.Results.is_secular;
            timelist=p.Results.timelist;
            
            %generate the spin_collection for this cluster including the central spin
            obj.spin_collection= model.phy.SpinCollection.SpinCollection();
            obj.spin_collection.spin_source=model.phy.SpinCollection.Strategy.FromSpinList([{obj.central_spin},obj.cluster_bath_spin]);
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
        
        function set_local_field(obj)
            
            
        end
                       
    end
    
end
