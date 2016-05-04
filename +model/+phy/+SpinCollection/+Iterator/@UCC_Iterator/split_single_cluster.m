function [clust_cell,size_list]=split_single_cluster(obj,index)
    import model.phy.SpinCollection.Strategy.FromSpinList
    import model.phy.SpinCollection.SpinCollection            
    import model.phy.SpinCollection.Iterator.ClusterIterator
    import model.phy.SpinCollection.Iterator.ClusterIteratorGen.CCE_Clustering
    
    clst_spin_list=obj.spin_collection.spin_list(index);
    clst_sc=SpinCollection( FromSpinList(clst_spin_list));    
    nspin=clst_sc.getLength;
    
    clu_para.cutoff=obj.CutOff;
    clu_para.max_order=2;
    cce1=CCE_Clustering(clst_sc, clu_para);
    clst_iter=ClusterIterator(clst_sc,cce1);
    connection_mat=clst_iter.index_generator.connection_matrix;
    
    idx=clst_iter.index_list{end,1};
    init_clst={idx};
    init_size=2;
    rowsum=sum(connection_mat,2);
    disconnected_spins=find(rowsum<1)';
    ndisc_spin=length( disconnected_spins);
    
    [clust_cell,size_list]=generate_cluster_cell(connection_mat,init_clst,init_size,idx,disconnected_spins);
    clust_cell(1)=[];
    size_list=size_list(2:end);
    disconnect_spin_cell=cell(1, ndisc_spin);

    if ~isempty(disconnected_spins)
        for kk=1:ndisc_spin
            disconnect_spin_cell{kk}=disconnected_spins(kk);   
        end 
        clust_cell=[clust_cell,disconnect_spin_cell];
        size_list=[size_list,ones(1,ndisc_spin)];
    end 
    
    if sum(size_list)~= nspin
        error('The total numbers of the spin in the cluster_cell are not the same as the total spin number!');
    end
    
    %transform to the original index
    ncell=length(clust_cell);
    for jj=1:ncell
        idx2trans=clust_cell{jj};
        transfered_idx=index(idx2trans);
        clust_cell{jj}=transfered_idx;        
    end

end

function [cluster_cell,cluster_size_list]=generate_cluster_cell(connection_mat,clst_cell,clst_size_list,idx4next_cluster,disconnected_spins)    
    nspin=size(connection_mat,1);
    spin_idx=add_elements(idx4next_cluster,connection_mat);
    cluster_size=length(spin_idx);
    
    cluster_cell=[clst_cell,{sort(spin_idx)}];
    cluster_size_list=[clst_size_list,cluster_size];
    
    
    SubSets=[];
    for kk=1:length(cluster_cell)
        SubSets=[SubSets,cluster_cell{kk}];
    end
    SubSets=unique(SubSets);
    superset=1:nspin;
    superset=setdiff(superset,disconnected_spins);% remove the spin which is connected to no one
    comple_set=sort( setdiff(superset,SubSets) );
    
    if ~isempty(comple_set)
        row_idx=comple_set(1);
        idx4next_cluster=[row_idx,find( connection_mat(row_idx,:) )];
        if isempty(idx4next_cluster)
           error('Something wrong the connection matrix.') 
        end
       [cluster_cell,cluster_size_list]= generate_cluster_cell(connection_mat,cluster_cell,cluster_size_list,idx4next_cluster,disconnected_spins);
    end
        
    
    
end

function new_idx=add_elements(idx,connection_mat)
    len0=length(idx);
    new_idx=idx;
    for kk=1:len0
       row_idx=idx(kk);
       row2add=connection_mat(row_idx,:);
       idx2add=find(row2add); 
       new_idx=[new_idx, idx2add];
       new_idx=unique(new_idx);
    end
    new_idx=unique(new_idx);
    if length(new_idx)>len0
        new_idx=add_elements(new_idx,connection_mat);
    end
    
end