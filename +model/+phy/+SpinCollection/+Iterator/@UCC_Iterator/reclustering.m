function reclustering(obj)
    import model.phy.SpinCollection.Strategy.FromSpinList
    import model.phy.SpinCollection.SpinCollection            
    import model.phy.SpinCollection.Iterator.ClusterIterator
    import model.phy.SpinCollection.Iterator.ClusterIteratorGen.CCE_Clustering
    
    if obj.Split_Option
        pos_list=find(obj.cluster_size_list > obj.cluster_size_uplimit);
        if ~isempty(pos_list)
            clusters2split=obj.index_list(pos_list);
            %remove the orginal larger clusters and cluster sizes
            obj.index_list(pos_list)=[];
            obj.cluster_size_list(pos_list)=[];
            ncluster2split=length(pos_list);
            for kk=1:ncluster2split
                indices=clusters2split{kk};
                %split the large clusters and generate new small clusters and size list
                [small_clusters,small_cluster_size]=obj.split_single_cluster(indices);                
                obj.index_list=[obj.index_list,small_clusters];
                obj.cluster_size_list=[obj.cluster_size_list,small_cluster_size];             
            end
            obj.sort_clusters();
            obj.CutOff=obj.CutOff-10;
            obj.reclustering();
            
         else
            obj.Split_Option=0;
            obj.cluster_size_list=obj.cluster_size_list';
            obj.index_list=obj.index_list';
            disp(['The size of largest disconnect spin clusters is ' num2str( max(obj.cluster_size_list) ) ' .']);
            disp(['The number of disconnect spin clusters is ' num2str( length(obj.index_list) ) ' .']);
        end
    end

end