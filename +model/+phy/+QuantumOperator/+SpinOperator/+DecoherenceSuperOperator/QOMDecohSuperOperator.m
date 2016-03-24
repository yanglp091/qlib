classdef QOMDecohSuperOperator < model.phy.QuantumOperator.MultiSpinSuperOperator.DecoherenceSuperOperator
    %QOMDecohSuperOperator generate a superoperator to describe the decoherence of the multi-spin system for its quantum optical master equation.
    % the total decoherence operator in the Liouville space is the summation of the decoherence operator 
    % for each single spin. But if the spin is not spin half, the corresponding decoherence operator is 
    % replaced by the zero matrix. 

    %   The decoherence operator of the single spin is of the Lindblad form:
    %  L rho = (1/2)*Gamma_vertical*[f(omega)+1](2sigma_{-}*rho*sigma_{+} - rho*sigma_{+}sigma_{-} - sigma_{+}sigma_{-}*rho)
    %          +(1/2)*Gamma_vertical*f(omega)*(2sigma_{+}*rho*sigma_{-} - rho*sigma_{-}sigma_{+} - sigma_{-}sigma_{+}*rho)
    %          +(1/2)*Gamma_parallel*[sigma_{z}*rho*sigma_{z} - rho]
    % where sigma is the Pauli operator of the spin, Gamma_vertical (Gamma_vertical) describe the vertical (parallel) coupling 
    % between the spin and its bath, omega is the Zeeman splitting of the spin, and f(omega)=1/[exp(hbar*omega/kT)] is the 
    % occupation number of the bosonic mode with frequency omega.
    
    
    properties
        decay_rate_list
        
    end
    
    methods
        function obj=QOMDecohSuperOperator(spin_collection,decay_rate_list)
            obj@model.phy.QuantumOperator.MultiSpinSuperOperator.DecoherenceSuperOperator(spin_collection);
            obj.name='decoherence_super_operator';            
            obj.decay_rate_list=decay_rate_list;
        end
        
        function generate_matrix(obj)
            obj.matrix=obj.calculate_matrix();
            obj.hasMatrix=1;
        end
        
        function L_matrix = calculate_matrix(obj)
            nspin=obj.spin_collection.getLength;
            L_matrix=0;   
            for kk=1:nspin
                spin=obj.spin_collection.spin_list{kk};
                if spin.dim>2
                    L_kk=0;
                    L_matrix=L_matrix+L_kk;
                elseif spin.dim==2
                    L_kk=obj.gen_superoperator_single(spin,kk);

                    L_matrix=sparse(L_matrix+L_kk);

                end
            end

        end
        
        function L_single=gen_superoperator_single(obj,spin,kk)
            Gamma_v=obj.decay_rate_list.Gamma_vertical_list(kk);% the vertical decay rate 
            Gamma_p=obj.decay_rate_list.Gamma_parallel_list(kk);% the parallel decay rate

            if Gamma_v >0 || Gamma_p > 0

                dim_tot=obj.spin_collection.getDim;
                spin_collection=obj.spin_collection;
                % generate ladder operators
                import model.phy.QuantumOperator.SpinOperator.Observable
                Sm=Observable(spin_collection, 'sigma-', ['1.0 * mat([0,0;1,0])_' num2str(kk)]);
                sigma_m=Sm.getMatrix;
                Sp=Observable(spin_collection, 'sigam+', ['1.0 * mat([0,1;0,0])_' num2str(kk)]);
                sigma_p=Sp.getMatrix;
                Sz=Observable(spin_collection, 'sigmaz', ['1.0 * mat([1,0;0,-1])_'  num2str(kk)]);
                sigma_z=Sz.getMatrix;


                % generate supperoperator for the current spin
                if Gamma_v>0            
                    % ocupation number of the bath mode
                    eigen_vals=spin.eigen_val;
                    omega=abs(eigen_vals(1)-eigen_vals(2));% the spin has only two eigen values
                    condition=model.phy.LabCondition.getCondition;
                    temperature=condition.getValue('temperature');
                    if temperature==0
                        f_omega=0;
                    else
                        f_omega=1/(exp(hbar*omega/kB/temperature)-1);
                    end

                    L_v=0.5*Gamma_v*(f_omega+1)*(2*kron(sigma_p',sigma_m)-kron(speye(dim_tot),sigma_p*sigma_m)-kron((sigma_p*sigma_m)',speye(dim_tot)))...
                        +0.5*Gamma_v*f_omega*(2*kron(sigma_m',sigma_p)-kron(speye(dim_tot),sigma_m*sigma_p)-kron((sigma_m*sigma_p)',speye(dim_tot)));
                    %In the kron product, we have used the Hermitian conjugation to replace the transpose operation. because all the matrix is real

                else
                   L_v=0; 
                end

                if Gamma_p>0
                    L_p=0.5*Gamma_p*(kron(sigma_z',sigma_z)-speye(dim_tot^2));
                else
                    L_p=0; 
                end
                L_single=sparse(L_v+L_p);
            else
               L_single=0;
            end

        end
    end
    
end

