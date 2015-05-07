classdef ZeemanInteraction < model.phy.SpinInteraction.AbstractSpinInteraction
    %ZEEMANINTERACTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj=ZeemanInteraction(iter, para)
            obj@model.phy.SpinInteraction.AbstractSpinInteraction(iter, para);
            obj.nbody=1;
        end
        
        function coeff=calculate_coeff(obj, spins)
            spin=spins{1};
            coeff=obj.para.B * spin.gamma;
        end
        function mat=calculate_matrix(obj, spins)
            spin=spins{1};
            coeff=obj.calculate_coeff(spins);
            mat=coeff*spin.sz;
        end
    end
    
end
