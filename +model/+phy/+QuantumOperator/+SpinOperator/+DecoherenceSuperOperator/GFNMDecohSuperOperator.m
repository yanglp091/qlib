classdef GFNMDecohSuperOperator < model.phy.QuantumOperator.SpinOperator.DecoherenceSuperOperator.AbstractDecoherenceOperator
    %GFNMDecohSuperOperator generate a superoperator to describe the decoherence of the multi-spin system.
    % The master equation is generated from a gradient field noise, i.e.,
    % all the spin experience the same noise field with different
    % amplitude: h_j(t)= B^z(x_j)*f(t). Here, B^z(x_j) is the amplitude of
    % the noise field along the z direction, f(t) is a stochastic function, whoes spectrum is normalized:
    % S(omega)=\int_{-\infty}^{\infty} ds <f(t)f(t-s)> exp(-i*omega*s);
    % If we assume <f(t)f(t-s)> = exp(-|s|/ tau_z)/(2*pi), the corresponding spectrum is
    % S(omega) = (1/pi)*tau_z/(omega^2*tau_z^2+1); % tau_z is the correlation time of the noise field.
    % \int_{-\infty}^{\infty} domega S(omega)=1; The integral of this spectrum is 1.
 
    % In this case, the nose fields exerted on different spins are
    % correlated, i.e., <h_i(t)h_j(t')> \neq = 0;
    % the total decoherence operator in the Liouville space is the summation of the decoherence operator 
    % for each single spin and all pairs. But if the spin is not spin half, the corresponding decoherence operator is 
    % replaced by the zero matrix. 

    %   The decoherence operator of the single spin is of the Lindblad form:
    %  L rho(t) = \sum_{j} Gamma_{j}*[sigma_{j}^{z}*rho*sigma_{j}^{z} - rho(t)] % the diagonal part
    %              \sum_{i<j} Gamma_{ij}*[sigma_{i}^{z}*rho(t)*sigma_{j}^{z} + sigma_{j}^{z}*rho(t)*sigma_{i}^{z}
    %              -sigma_{i}^{z}*sigma_{j}^{z}*rho(t) - rho(t)*sigma_{i}^{z}*sigma_{j}^{z}] 
    % where sigma is the Pauli operator of the spin, 
    % Gamma_{j}= [gamm_j*B^z(x_j)]^2 * S(0)/2 = |B^z(x_j)|^2 * tau_z/(2*pi);
    % Gamma_{ij}= [gamm_i*B^z(x_i)*gamm_j*B^z(x_j)] * S(0)/4 = [gamm_i*B^z(x_i)*gamm_j*B^z(x_j)] * tau_z/(2*pi);
    % gamma_j  is the gyromagnetic ratio of the jth spin
    
        % parameter contants: 
    % parameter.CorrelationTime :: The correlation time of the noise field %%%%% 
   
    properties
        parameter
        axis
    end
    
    methods
        function obj=GFNMDecohSuperOperator(spin_collection,para, matrix_strategy)
            obj@model.phy.QuantumOperator.SpinOperator.DecoherenceSuperOperator.AbstractDecoherenceOperator(spin_collection);
            obj.parameter=para;%The correlation time of the noise field;
            if nargin > 2
                obj.matrix_strategy = matrix_strategy;
            else
                obj.matrix_strategy = [];
            end

            obj.name='decoherence_super_operator'; 
            obj.generate_interaction_list;
        end
        
        function generate_interaction_list(obj)
            obj.add_local_paralle_decay;
            obj.add_cross_paralle_decay;
            
        end
        
        function add_local_paralle_decay(obj)           
            iter=model.phy.SpinCollection.Iterator.SingleSpinIterator(obj.spin_collection);
            nItem=iter.getLength;
            para.DecayRateList=cell(nItem,1);
            for kk=1:nItem
                spins=iter.getItem(kk);
                correlation_time=obj.parameter.CorrelationTime;
                gradient_field=spins{1}.local_field;
                gammaj=spins{1}.gamma;%the gyromagnetic ratio
                decay_rate=gradient_field*gradient_field'*(gammaj^2)*correlation_time/pi;%
                para.DecayRateList{kk}=decay_rate;
            end
            local_decay_term=model.phy.SpinInteraction.DecoherenceInteraction.LocalParallelDecayInteraction(para,iter);
            obj.addInteraction(local_decay_term);
            
        end
               
        function add_cross_paralle_decay(obj)
            iter=model.phy.SpinCollection.Iterator.PairSpinIterator(obj.spin_collection);
            nItem=iter.getLength;
            para.DecayRateList=cell(nItem,1);
            for kk=1:nItem
                spins=iter.getItem(kk);
                
                 % gamma_{ij}=[gamm_i*B^z(x_i)*gamm_j*B^z(x_j)] * tau_z/(2*pi);
                correlation_time=obj.parameter.CorrelationTime;
                g_field1=spins{1}.local_field;
                g_field2=spins{2}.local_field;
                gamma1=spins{1}.gamma;%the gyromagnetic ratio
                gamma2=spins{2}.gamma;%the gyromagnetic ratio
                decay_rate=g_field1*g_field2'*gamma1*gamma2*correlation_time/pi/2;% the parallel decay rate
                para.DecayRateList{kk}=decay_rate;                
            end

            cross_decay_term=model.phy.SpinInteraction.DecoherenceInteraction.CrossParallelDecayInteraction(para,iter);
            obj.addInteraction(cross_decay_term);
        end
        
    
    end
    
end

