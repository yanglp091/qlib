classdef LocalVerticalDecayInteraction < model.phy.SpinInteraction.AbstractSpinInteraction
    %LOCALVERTICALDECAYINTERACTION Summary of this class goes here
    %   (1/2)*Gamma_vertical*factor1(2sigma_{-}*rho*sigma_{+} - rho*sigma_{+}sigma_{-} - sigma_{+}sigma_{-}*rho)
    %          +(1/2)*Gamma_vertical*factor2*(2sigma_{+}*rho*sigma_{-} - rho*sigma_{-}sigma_{+} - sigma_{-}sigma_{+}*rho)
    %
    % For quantum optical master equation, factor1=[f(omega)+1]; factor2=f(omega). %%%
    % Here, f(omega)=1/[exp(hbar*omega/kT)] is the 
    % occupation number of the bosonic mode with frequency omega.
    % But sometimes we want factor1 = factor2
    % Here, we want the parameter of this class must be the following form:
    % parameter.coeffs :: gives the vertical decay rate list for all the spins in the spin collection %%%
    % parameter.FactorList1 :: gives the factor1 list for all the spins in the spin collection%%%
    % parameter.FactorList1 :: gives the factor2 list for all the spins in the spin collection%%%
    
    
    properties
        
    end
    
    methods
        function obj=LocalVerticalDecayInteraction(para,iter)
            obj@model.phy.SpinInteraction.AbstractSpinInteraction(para, iter);
            obj.nbody=1;
        end
        
        function skp=single_skp_term(obj)
            spin=obj.iter.currentItem{1};
            dim=spin.dim;
            if dim ~=2
               error('This class can only generate the decoherence operator for spin-half system.'); 
            end
            
            idx=obj.iter.currentIndex();
            cursor=obj.iter.getCursor();            
            sigmap=spin.sqp;sigmam=spin.sqm;

           [Gamma_v,factor1,factor2]=obj.calculate_coeff(cursor);% the cursor is the same as the index
            
            mat=Gamma_v*factor1*kron(sigmap.',sigmam); 
            Term1=obj.supper_kron_prod(1, idx, {mat});
            
            mat=-0.5*Gamma_v*factor1*kron( (sigmap*sigmam).',speye(dim)); 
            Term2=obj.supper_kron_prod(1, idx, {mat});
            
            mat=-0.5*Gamma_v*factor1*kron(speye(dim), sigmap*sigmam); 
            Term3=obj.supper_kron_prod(1, idx, {mat});
            
            mat=Gamma_v*factor2*kron(sigmam.',sigmap);
            Term4=obj.supper_kron_prod(1, idx, {mat});
            
            mat=-0.5*Gamma_v*factor2*kron( (sigmam*sigmap).',speye(dim)); 
            Term5=obj.supper_kron_prod(1, idx, {mat});
            
            mat=-0.5*Gamma_v*factor2*kron(speye(dim), sigmam*sigmap);
            Term6=obj.supper_kron_prod(1, idx, {mat});
            
            skp=Term1+Term2+Term3+Term4+Term5+Term6;
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
                    L_kk=obj.gen_superoperator_single_term(kk);

                    L_matrix=sparse(L_matrix+L_kk);

                end
            end
        end
        
        function L_single=gen_superoperator_single_term(obj,idx)
            [Gamma_v,factor1,factor2]=obj.calculate_coeff(idx);% the cursor is the same as the index

            if Gamma_v >0
                dim_tot=obj.iter.spin_collection.getDim;
                sc=obj.iter.spin_collection;
                % generate ladder operators
                import model.phy.QuantumOperator.SpinOperator.Observable
                import model.phy.QuantumOperator.MatrixStrategy.FromKronProd
                
                strategy=FromKronProd();
                Sm=Observable(sc,strategy, 'sigma-', ['1.0 * sqm_' num2str(idx)]);
                sigma_m=full(Sm.getMatrix);
                Sp=Observable(sc,strategy, 'sigam+', ['1.0 * sqp_' num2str(idx)]);
                sigma_p=full(Sp.getMatrix);
          

                L_single=0.5*Gamma_v*factor1*(2*kron(sigma_p.',sigma_m)-kron(speye(dim_tot),sigma_p*sigma_m)-kron((sigma_p*sigma_m).',speye(dim_tot)))...
                    +0.5*Gamma_v*factor2*(2*kron(sigma_m.',sigma_p)-kron(speye(dim_tot),sigma_m*sigma_p)-kron((sigma_m*sigma_p).',speye(dim_tot)));
            else
               L_single=0;
            end
            
        end
        function [Gamma_v,factor1,factor2]=calculate_coeff(obj,idx)
            Gamma_v=obj.parameter.DecayRateList{idx};
            factor1=obj.parameter.FactorList1{idx};
            factor2=obj.parameter.FactorList2{idx};            
        end
        
        function dataCell=data_cell(obj)

        end
    end
    
end

