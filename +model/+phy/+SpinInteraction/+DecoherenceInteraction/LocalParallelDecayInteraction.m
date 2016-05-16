classdef LocalParallelDecayInteraction < model.phy.SpinInteraction.AbstractSpinInteraction
    %ABSTRACTDECOHINTERACTION Summary of this class goes here
    %   L rho(t) = \sum_{j} Gamma_{j}*[sigma_{j}^{z}*rho*sigma_{j}^{z} - rho(t)]
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
            
            mat1=2*spin.sqz;
            dim=size(mat1,1);
            if dim ~=2
               error('This class can only generate the decoherence operator for spin-half system.'); 
            end
            coeff=obj.calculate_coeff(cursor);
            
            mat=coeff*(kron(mat1.',mat1)-speye(dim^2));
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
                Sigmaz=Observable(sc,strategy, 'sigmaz', ['2.0 * sqz_'  num2str(idx)]);

                matrix_z=full(Sigmaz.getMatrix);
                L_p=coeff*(kron(matrix_z.',matrix_z)-speye(dim_tot^2));
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

