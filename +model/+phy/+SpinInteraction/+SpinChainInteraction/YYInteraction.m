classdef YYInteraction < model.phy.SpinInteraction.SpinChainInteraction.AbstractSpinChainInteraction
    %XYINTERACTION Summary of this class goes here
    % The input parameter can contain the following iterms:
    % para.interaction:: To input the interaction strength matrix.
    % para.AddChain= 0 or 1:: To generate a chain iterator
    % para.AddRing=0 or 1:: To generate a ring iterator
    % para.AddIndexList=0 or 1:: To generate a iterator with the given index_list.
    
    properties
    end
    
    methods
        function obj=YYInteraction(spin_collection, para)
            nspin=spin_collection.getLength;
            has_iter=0;
            
            try
                parameter.interaction=para.interaction;
            catch
                parameter.interaction=ones(nspin);
            end
            
            try
                AddChain=para.AddChain;
            catch
                AddChain=0;
            end
            if AddChain
                iter=model.phy.SpinCollection.Iterator.ChainNeighbouringIterator(spin_collection);
                has_iter=1;
            end
            
            try
                AddRing=para.AddRing;
            catch
                AddRing=0;
            end
            if AddRing
                iter=model.phy.SpinCollection.Iterator.ChainNeighbouringIterator(spin_collection,'ring_index');
                has_iter=1;
            end
            try 
                AddIndexList=para.AddIndexList;
            catch
                AddIndexList=0;
            end
            if AddIndexList
                iter=model.phy.SpinCollection.Iterator.ChainNeighbouringIterator(spin_collection);
                iter.index_list=para.IndexList;
                has_iter=1;
            end
            
            if ~has_iter
                iter=model.phy.SpinCollection.Iterator.ChainNeighbouringIterator(spin_collection);
            end
                
            obj@model.phy.SpinInteraction.SpinChainInteraction.AbstractSpinChainInteraction(parameter, iter);
            obj.nbody=2;
        end
        function skp=single_skp_term(obj)
            spins=obj.iter.currentItem();
            idx=obj.iter.currentIndex();
            spin1=spins{1}; spin2=spins{2};
            coeff=obj.calculate_coeff(idx);
            
            mat1=spin1.sy; mat2=spin2.sy;
            yyTerm=obj.kron_prod(coeff, idx, {mat1, mat2});
            skp=yyTerm;
        end
        
        function mat=calculate_matrix(obj)
            spins=obj.iter.currentItem();
            spin1=spins{1}; spin2=spins{2};
            idx=obj.iter.currentIndex();
            coeff=obj.calculate_coeff(idx);
            mat= coeff*kron(spin1.sy,spin2.sy);
        end
        
        function dataCell=data_cell(obj)
            nPair=obj.iter.getLength;
            nTerms=nPair;
            dataCell=cell(1, nTerms);
            
            obj.iter.setCursor(1);
            for ii=1:nPair
                spins=obj.iter.currentItem();
                spin_idx=obj.iter.currentIndex();
                spin1=spins{1}; spin2=spins{2};
                index1=spin_idx(1); index2=spin_idx(2);
                
                coeff=obj.calculate_coeff(spins);
                
                dataCell{ii}={coeff, obj.nbody, ...
                    index1, spin1.dim, reshape(spin1.sy, [], 1), ...
                    index2, spin2.dim, reshape(spin2.sy, [], 1)
                    };
                obj.iter.nextItem;
            end
        end
        
    end
    
end

