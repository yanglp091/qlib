classdef CrossParallelDecayInteraction < model.phy.SpinInteraction.AbstractSpinInteraction
    %PARALLELCROSSDECAYINTERACTION Summary of this class goes here
    %   L=\sum_{i<j} (1/2)*Gamma_{ij}*[sigma_{i}^{z}*rho(t)*sigma_{j}^{z} + sigma_{j}^{z}*rho(t)*sigma_{i}^{z}
    %              -sigma_{i}^{z}*sigma_{j}^{z}*rho(t) - rho(t)*sigma_{i}^{z}*sigma_{j}^{z}] 
    % parameter.DecayRateList :: gives the vertical decay rate list for all the spins in the spin collection %%%
    properties
    end
    
    methods
        function obj=CrossParallelDecayInteraction(para,iter)
            obj@model.phy.SpinInteraction.AbstractSpinInteraction(para, iter);
            obj.nbody=2;
        end
        
        function skp=single_skp_term(obj)
            spins=obj.iter.currentItem;
            idx=obj.iter.currentIndex();
            cursor=obj.iter.getCursor();
            spin1=spins{1}; spin2=spins{2};
            sigmaz1=2*spin1.sz;sigmaz2=2*spin2.sz;
            dim1=size(sigmaz1,1);dim2=size(sigmaz2,1);
            coeff=obj.calculate_coeff(cursor);
            
            mat1=0.5*coeff*kron(speye(dim1),sigmaz1); mat2=kron(sigmaz2,speye(dim2));
            Term1=obj.supper_kron_prod(1.0, idx, {mat1, mat2});
            
            mat1=0.5*coeff*kron(sigmaz1,speye(dim1)); mat2=kron(speye(dim2),sigmaz2);
            Term2=obj.supper_kron_prod(1.0, idx, {mat1, mat2});
            
            mat1=-0.5*coeff*kron(speye(dim1),sigmaz1); mat2=kron(speye(dim2),sigmaz2);
            Term3=obj.supper_kron_prod(1.0, idx, {mat1, mat2});
            
            mat1=-0.5*coeff*kron(sigmaz1,speye(dim1)); mat2=kron(sigmaz2,speye(dim2));
            Term4=obj.supper_kron_prod(1.0, idx, {mat1, mat2});
            
            skp=Term1+Term2+Term3+Term4;
        end
        
        function L_matrix=calculate_matrix(obj)
            nspin=obj.iter.spin_collection.getLength;
            L_matrix=0;
            if nspin<2
               % do nothing in this case.
            else               
                nPair=obj.iter.getLength;
                for kk=1:nPair
                    obj.iter.setCursor(kk);
                    spin_idx=obj.iter.currentIndex();
                    L_pair=obj.gen_superoperator_pair(spin_idx,kk);% the idx is the index of the spins, kk is the cursor of the iter
                    L_matrix=L_matrix+L_pair;                   
                end
            end
        end
        function L_pair=gen_superoperator_pair(obj,index,cursor)  
            coeff=obj.calculate_coeff(cursor); 
            if coeff>0
                dim_tot=obj.iter.spin_collection.getDim;
                spin_collection=obj.iter.spin_collection;
                % generate ladder operators
                import model.phy.QuantumOperator.SpinOperator.Observable
                import model.phy.QuantumOperator.MatrixStrategy.FromKronProd

                strategy=FromKronProd();
                S1z=Observable(spin_collection,strategy, 'sigmaz', ['1.0 * mat([1,0;0,-1])_'  num2str(index(1))]);
                sigma_1z=full(S1z.getMatrix);
                S2z=Observable(spin_collection,strategy, 'sigmaz', ['1.0 * mat([1,0;0,-1])_'  num2str(index(2))]);
                sigma_2z=full(S2z.getMatrix);
                S12z=Observable(spin_collection,strategy, 'sigmaz', ['1.0 * mat([1,0;0,-1])_'  num2str( index(1) )...
                    ' * mat([1,0;0,-1])_'  num2str(index(2))]);
                sigma_12z=full(S12z.getMatrix);


                L_p=0.5*coeff*( kron(sigma_1z',sigma_2z)+kron(sigma_2z',sigma_1z)...
                      -kron(sigma_12z',speye(dim_tot))-kron(speye(dim_tot),sigma_12z) );

                L_pair=sparse(L_p);
            else
                L_pair=0;
            end

        end
        
        function coeff=calculate_coeff(obj,idx)
            coeff=obj.parameter.DecayRateList{idx};
        end
        
        function dataCell=data_cell(obj)

        end
        
        
        
        
    end
    
end

