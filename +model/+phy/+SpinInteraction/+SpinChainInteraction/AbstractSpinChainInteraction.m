classdef AbstractSpinChainInteraction < model.phy.SpinInteraction.AbstractSpinInteraction
    %SPINCHAININTERACTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj=AbstractSpinChainInteraction(para, iter)
            obj@model.phy.SpinInteraction.AbstractSpinInteraction(para, iter);
        end
        
        function coeff=calculate_coeff(obj, spin_idx)
            if length(spin_idx)==2
                coeff=obj.parameter.interaction(spin_idx(1),spin_idx(2));
            elseif length(spin_idx)==1
                coeff=obj.parameter.interaction(spin_idx);
            else
                error('the number of the spins is larger than 2...');
            end
            
        end
        
        function mat=calculate_matrix(obj)
            spins=obj.iter.currentItem();
            spin1=spins{1}; spin2=spins{2};
            coeff=obj.calculate_coeff(spins);
            mat= coeff*kron(spin1.sx,spin2.sx)...
                +coeff*kron(spin1.sy,spin2.sy);
        end
        
    end
    
end

