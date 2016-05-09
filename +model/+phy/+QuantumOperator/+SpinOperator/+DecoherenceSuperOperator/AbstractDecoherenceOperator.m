classdef AbstractDecoherenceOperator < model.phy.QuantumOperator.MultiSpinSuperOperator
    %DECOHERENCESUPEROPERATOR generate a superoperator to describe the decoherence of the multi-spin system.
    % the total decoherence operator in the Liouville space is the summation of the decoherence operator 
    % for each single spin. But if the spin is not spin half, the corresponding decoherence operator is 
    % replaced by the zero matrix. 

    % There are two ways to generate a supper operator from an operator in the Hilbert space: 
    % The first one is : get the matrix of the operator in the
    % Hilbert space, then transform the matrix to the corresponding matrix
    % in the Liouville space. The rules are: 
    % A*rho <---> E \otimes A * rho_vec
    % rho*A <---> A^{T} \otimes E * rho_vec
    % A*rho*B <---> B^{T} \otimes A * rho_vec
    % E is the identity matrix with the same dimension of the total dimension of the system in the Hilbert space,
    % rho_vec=rho(:), "T" is the transpose operation.
    
    %The second one is: transfer the operator of each spin to the Liouville
    %space, the supper operator is the  Kronecker tensor product of the
    %matice of each spin. The rules are (e.g., two spin system):
    % A1*rho = (A1 \otimes E2)*rho_vec <---> (E1 \otimes A1) \times (E2 \otimes E2) * rho_vec
    % rho*A1 = rho_vec*(A1 \otimes E2) <---> (A1^T \otimes E1) \times (E2 \otimes E2) * rho_vec
    % A1rho*B2 = (A1 \otimes E2)*rho_vec*(E1 \otimes B2) <---> (E1 \otimes A1) \otimes (B2^T \otimes E2) * rho_vec
    % rho_vec=rho1(:) \otimes rho2(:);
    
    properties       
    end
    
    methods
        function obj=AbstractDecoherenceOperator(spin_collection)
            obj@model.phy.QuantumOperator.MultiSpinSuperOperator(spin_collection);
        end
                
        function generate_matrix(obj, mat)
            if nargin > 1
                obj.matrix=mat;
            elseif isempty(obj.matrix_strategy)
                %generate the supperoperator by the 1st method
                obj.matrix=obj.get_matrix_prod_space;                
            else
                %generate the supperoperator by the 2nd method
                obj.matrix_strategy.initialize(obj);
                obj.matrix=obj.matrix_strategy.calculate_matrix();
            end
            obj.hasMatrix=1;
        end
        
        function matrix=get_matrix_prod_space(obj)
            matrix=sparse(0);
            nInteraction=length(obj.interaction_list);
            
            for k=1:nInteraction
                interaction_k=obj.interaction_list{k};
                interaction_matrix=interaction_k.calculate_matrix();
                matrix=matrix+interaction_matrix;
            end
        end
    end
    
    
end

