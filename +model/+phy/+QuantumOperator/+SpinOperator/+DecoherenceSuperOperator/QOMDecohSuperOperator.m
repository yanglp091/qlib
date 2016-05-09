classdef QOMDecohSuperOperator < model.phy.QuantumOperator.SpinOperator.DecoherenceSuperOperator.AbstractDecoherenceOperator
    %QOMDecohSuperOperator generate a superoperator to describe the decoherence of the multi-spin system for its quantum optical master equation.
    % the total decoherence operator in the Liouville space is the summation of the decoherence operator 
    % for each single spin. But if the spin is not spin half, the corresponding decoherence operator is 
    % replaced by the zero matrix. Here, we just consider the independent bath  case. %%%

    %   The decoherence operator of the single spin is of the Lindblad form:
    %  L rho = (1/2)*Gamma_vertical*[f(omega)+1](2sigma_{-}*rho*sigma_{+} - rho*sigma_{+}sigma_{-} - sigma_{+}sigma_{-}*rho)
    %          +(1/2)*Gamma_vertical*f(omega)*(2sigma_{+}*rho*sigma_{-} - rho*sigma_{-}sigma_{+} - sigma_{-}sigma_{+}*rho)
    %          +(1/2)*Gamma_parallel*[sigma_{z}*rho*sigma_{z} - rho]
    % where sigma is the Pauli operator of the spin, Gamma_vertical (Gamma_vertical) describe the vertical (parallel) coupling 
    % between the spin and its bath, omega is the Zeeman splitting of the spin, and f(omega)=1/[exp(hbar*omega/kT)] is the 
    % occupation number of the bosonic mode with frequency omega.
    % parameter contants: 
    % parameter.AddVerticalDecay :: Option to determine add  vertical decay or not %%%%% 
    % parameter.VerticalDecayRateList :: A list of vertical decay rates for the spins respectively  %%%%% 
    % parameter.AddParallelDecay :: Option to determine add  vertical decay or not %%%%% 
    % parameter.ParallelDecayRateList :: A list of parallel decay rates for the spins respectively  %%%%% 
    properties        
        parameter
    end
    
    methods
        function obj=QOMDecohSuperOperator(spin_collection,para,matrix_strategy)
            obj@model.phy.QuantumOperator.SpinOperator.DecoherenceSuperOperator.AbstractDecoherenceOperator(spin_collection);
            obj.name='decoherence_super_operator';            
            obj.parameter=para;
            
            if nargin > 2
                obj.matrix_strategy = matrix_strategy;
            else
                obj.matrix_strategy = [];
            end
            obj.name='decoherence_super_operator'; 
            obj.generate_interaction_list;
        end
        
        function generate_interaction_list(obj)
            if obj.parameter.AddVerticalDecay
                obj.add_local_vertical_decay;
            end
            if obj.parameter.AddParallelDecay
                obj.add_local_paralle_decay;
            end
            
        end
        
       function add_local_vertical_decay(obj)           
            iter=model.phy.SpinCollection.Iterator.SingleSpinIterator(obj.spin_collection);
            nItem=iter.getLength;
            decay_rate=obj.parameter.VerticalDecayRateList;
            para.DecayRateList=cell(nItem,1);
            para.FactorList1=cell(nItem,1);
            para.FactorList2=cell(nItem,1);
            for kk=1:nItem
                para.DecayRateList{kk}=decay_rate(kk);
                % ocupation number of the bath mode                
                condition=model.phy.LabCondition.getCondition;
                temperature=condition.getValue('temperature');
                if temperature==0
                    factor=0;
                else
                    spin=iter.getItem(kk);
                    eigen_vals=spin{1}.eigen_val;
                    omega=abs(eigen_vals(1)-eigen_vals(2));% the spin has only two eigen values
                    factor=1/(exp(hbar*omega/kB/temperature)-1);
                end
                para.FactorList1{kk}=factor+1;
                para.FactorList2{kk}=factor;
            end
            
            local_decay_term=model.phy.SpinInteraction.DecoherenceInteraction.LocalVerticalDecayInteraction(para,iter);
            obj.addInteraction(local_decay_term);
            
        end
        
        function add_local_paralle_decay(obj)           
            iter=model.phy.SpinCollection.Iterator.SingleSpinIterator(obj.spin_collection);
            nItem=iter.getLength;
            decay_rate=obj.parameter.ParallelDecayRateList;
            para.DecayRateList=cell(nItem,1);
            for kk=1:nItem
                para.DecayRateList{kk}=decay_rate(kk);
            end
            
            local_decay_term=model.phy.SpinInteraction.DecoherenceInteraction.LocalParallelDecayInteraction(para,iter);
            obj.addInteraction(local_decay_term);
            
        end
               
        
    end
    
end

