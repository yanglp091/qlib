classdef ZeemanInteraction < model.phy.SpinInteraction.AbstractSpinInteraction
    %ZEEMANINTERACTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj=ZeemanInteraction(spin_collection, varargin)
            iter_class=@model.phy.SpinCollection.Iterator.SingleSpinIterator;
            for k=1:length(varargin)
                if isa(varargin{k}, 'model.phy.SpinCollection.SpinCollectionIterator')
                    iter_class=varargin{k};
                end
            end
            iter=iter_class(spin_collection);
            
            condition=model.phy.LabCondition.getCondition;
            parameter.B=condition.getValue('magnetic_field');
            obj@model.phy.SpinInteraction.AbstractSpinInteraction(parameter, iter);
            obj.nbody=1;
        end
        
        function skp=single_skp_term(obj)
            spin=obj.iter.currentItem{1};
            idx=obj.iter.currentIndex();
            coeff=obj.calculate_coeff(spin);
            
            mat=coeff(1)*spin.sx + coeff(2)*spin.sy + coeff(3)*spin.sz;
            skp=obj.kron_prod(1, idx, {mat});
        end
        
        function coeff=calculate_coeff(obj, spin)
            coeff=-(obj.parameter.B+spin.local_field) * spin.gamma;
        end
        
        function mat=calculate_matrix(obj)
            spin=obj.iter.currentItem{1};
            coeff=obj.calculate_coeff(spin);
            mat=coeff(1)*spin.sx + coeff(2)*spin.sy + coeff(3)*spin.sz;
            
            % add the zero-field splitting term
            if spin.ZFS
                principle_axis=spin.principle_axis;
                principle_axis=principle_axis/norm(principle_axis);
                mat1=principle_axis(1)*spin.sx+principle_axis(2)*spin.sy+principle_axis(3)*spin.sz;
                mat2=spin.ZFS*mat1*mat1;
                mat=mat+mat2;
            end
        end
        
        function dataCell=data_cell(obj)
            nTerms=obj.iter.getLength;
            dataCell=cell(1, nTerms);
            
            obj.iter.setCursor(1);
            for ii=1:nTerms
                spin=obj.iter.currentItem{1};
                spin_idx=obj.iter.currentIndex;

                mat=obj.calculate_matrix();
                
                dataCell{ii}={1.0, obj.nbody, spin_idx, spin.dim, reshape(mat, [], 1)};
                obj.iter.nextItem;
            end
        end
        
    end
    
end

