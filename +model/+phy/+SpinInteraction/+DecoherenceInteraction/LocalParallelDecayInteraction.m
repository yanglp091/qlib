classdef LocalParallelDecayInteraction < model.phy.SpinInteraction.AbstractSpinInteraction
    %ABSTRACTDECOHINTERACTION Summary of this class goes here
    %   L rho(t) = \sum_{j} (1/2)*Gamma_{j}*[sigma_{j}^{z}*rho*sigma_{j}^{z} - rho(t)]
    % Gamma_j=1/T_{2,j} is the decay rate of the phase
    % parameter.DecayRateList :: gives the vertical decay rate list for all the spins in the spin collection %%%
    properties
        
    end
    
    methods
        function obj=LocalParallelDecayInteraction(para,iter)
            obj@model.phy.SpinInteraction.AbstractSpinInteraction(para, iter);
            obj.nbody=1;
        end
        
        function skp=single_skp_term(obj)
            spin=obj.iter.currentItem{1};
            idx=obj.iter.currentIndex();
            cursor=obj.iter.getCursor();            
            mat1=2*spin.sz;% we use sigma_z for spin-half instead of S_z
            dim=size(mat1,1);
            if dim ~=2
               error('This class can only generate the decoherence operator for spin-half system.'); 
            end
            coeff=obj.calculate_coeff(cursor);
            
            mat=0.5*coeff*(kron(mat1,mat1)-speye(dim^2));
            skp=obj.supper_kron_prod(1, idx, {mat});
        end
        
        
        function L_matrix = calculate_matrix(obj)
            nspin=obj.iter.spin_collection.getLength;
            L_matrix=0;   
            for kk=1:nspin
                spin=obj.iter.spin_collection.spin_list{kk};
                if spin.dim>2
                    L_kk=0;
                    L_matrix=L_matrix+L_kk;
                elseif spin.dim==2
                    L_kk=obj.gen_superoperator_single_term(spin,kk);

                    L_matrix=sparse(L_matrix+L_kk);

                end
            end

        end
        function L_p=gen_superoperator_single_term(obj,spin,idx)
            coeff=obj.calculate_coeff(idx);% the cursor is the same as the index
            dim=spin.dim;
            if dim ~=2
               error('This class can only generate the decoherence operator for spin-half system.'); 
            end
            
            if coeff>0
                dim_tot=obj.iter.spin_collection.getDim;
                sc=obj.iter.spin_collection;
                % generate ladder operators
                import model.phy.QuantumOperator.SpinOperator.Observable
                import model.phy.QuantumOperator.MatrixStrategy.FromKronProd
                
                strategy=FromKronProd();                
                Sz=Observable(sc,strategy, 'sigmaz', ['1.0 * mat([1,0;0,-1])_'  num2str(idx)]);
                sigma_z=full(Sz.getMatrix);
                %In the kron product, we have used the Hermitian conjugation to replace the transpose operation,
                %because all the matrix is real.
                L_p=0.5*coeff*(kron(sigma_z',sigma_z)-speye(dim_tot^2));
                L_p=sparse(L_p);
            else
                L_p=0; 
            end
            
        end
        function coeff=calculate_coeff(obj,idx)
            coeff=obj.parameter.DecayRateList{idx};
        end
        
        function dataCell=data_cell(obj)

        end
        
    end
    
end

