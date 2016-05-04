classdef UCC_Iterator < model.phy.SpinCollection.SpinCollectionIterator
    %UCCEITERATOR : This iterator is used to generate uncorrelated large clusters
    % This iterator is used to solve the decoherence issue of the center
    % spin embeded in a spin bath
    
    properties
        cluster_size_list
        Split_Option=1;
        Combine_Option=1;
        CutOff
        cluster_size_uplimit
        cluster_size_lowlimit
    end
    
    methods
        function obj = UCC_Iterator(spin_collection,cutoff,upper_limit,lower_limit)
            obj@model.phy.SpinCollection.SpinCollectionIterator(spin_collection);
            obj.CutOff=cutoff;
            obj.cluster_size_uplimit=upper_limit;
            obj.cluster_size_lowlimit=lower_limit;
        end
        function index_list = index_gen(obj)
            nspin=obj.spin_collection.getLength;
            index_list={1:nspin};
            obj.cluster_size_list=nspin;
        end
        
        function sort_clusters(obj)
            [obj.cluster_size_list,orders]=sort(obj.cluster_size_list,'descend');
            obj.index_list=obj.index_list(orders);            
        end
        
    end
    
end

