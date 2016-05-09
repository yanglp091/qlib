classdef AbstractSpinInteraction < handle
    %SPININTERACTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nbody
        parameter
        iter
    end
    
    methods
        function obj=AbstractSpinInteraction(para, iter)
            obj.parameter=para;
            obj.iter=iter;
        end
        
        function res=isConsistent(obj, spin_collection)
            res=(spin_collection==obj.iter.spin_collection);
        end
        
        %This function generate a operator for multi-body system in Hilbert Space,
        % e.g., a single term C2*s_{z,2} is the angular operator of the 2nd spin in three-spin system
        % the full term should be C2*eye(s_1) \otimes S_{z,2} \otimes eye(s_3)
        % so we need to know the matrix of S_{z,2} and the dimenson of the
        % identity matrix before and after S_{z,2}
        function kp=kron_prod(obj, coeff, idx, mat)
            len=length(idx);          
            dim_list1= obj.iter.spin_collection.dim_compression(idx);
            mat_cell1 = num2cell(ones(1, 2*len+1));
            for ii=1:len
                mat_cell1{2*ii}=mat{ii};
            end

            dim_list=dim_list1(dim_list1>1);
            mat_cell=mat_cell1(dim_list1>1);
            kp=KronProd(mat_cell, fliplr(1:length(mat_cell)), fliplr(dim_list), coeff);            
        end
        
        %This function generate a supperoperator for multi-body system in Liouville Space,
        % e.g., a single term C2*s_{z,2} is the angular operator of the 2nd spin in three-spin system
        % the full term should be C2*eye(s_1) \otimes S_{z,2} \otimes eye(s_3)
        % so we need to know the matrix of S_{z,2} and the dimenson of the
        % identity matrix before and after S_{z,2} in the Liouville space
        function kp=supper_kron_prod(obj, coeff, idx, mat)
            len=length(idx);          
            dim_list1= obj.iter.spin_collection.dim_compression(idx);
            mat_cell1 = num2cell(ones(1, 2*len+1));
            for ii=1:len
                mat_cell1{2*ii}=mat{ii};
            end

            dim_list=dim_list1(dim_list1>1).^2;
            mat_cell=mat_cell1(dim_list1>1);
            kp=KronProd(mat_cell, fliplr(1:length(mat_cell)), fliplr(dim_list), coeff);            
        end
        
        function skp=skp_form(obj)
            obj.iter.setCursor(1);
            skp=obj.single_skp_term;
            while ~obj.iter.isLast()   
                obj.iter.moveForward();
                skp=skp+obj.single_skp_term();
            end
        end
    end
    
    methods (Abstract)
        calculate_coeff(obj, item);
        calculate_matrix(obj);
        single_skp_term(obj);
        %data_cell(obj);
    end
    
end

