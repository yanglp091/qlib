classdef dipInteraction < model.phy.SpinInteraction.SpinChainInteraction.AbstractSpinChainInteraction
    %DQTINTERACTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj=dipInteraction(spin_collection, para)
            nspin=spin_collection.getLength;
            n_interaction=size(para.interaction,1);
            if n_interaction==nspin-1
                iter=model.phy.SpinCollection.Iterator.ChainNeighbouringIterator(spin_collection);
            elseif n_interaction==nspin
                iter=model.phy.SpinCollection.Iterator.ChainNeighbouringIterator(spin_collection,'ring_index');
            else
                erro('wrong size of the interaction parameter!');
            end
            obj@model.phy.SpinInteraction.SpinChainInteraction.AbstractSpinChainInteraction(para, iter);
            obj.nbody=2;
        end
        
        function skp=single_skp_term(obj)
            spins=obj.iter.currentItem();
            idx=obj.iter.currentIndex();
            spin1=spins{1}; spin2=spins{2};
            coeff=obj.calculate_coeff(idx);
            
            mat1=spin1.sx; mat2=spin2.sx;
            xxTerm=obj.kron_prod(-coeff, idx, {mat1, mat2});
            
            mat1=spin1.sy; mat2=spin2.sy;
            yyTerm=obj.kron_prod(-coeff, idx, {mat1, mat2});
            
            
            mat1=spin1.sz; mat2=spin2.sz;
            zzTerm=obj.kron_prod(2*coeff, idx, {mat1, mat2});
            
            skp=xxTerm+yyTerm+zzTerm;
        end
        
        function mat=calculate_matrix(obj)
            spins=obj.iter.currentItem();
            spin1=spins{1}; spin2=spins{2};
            idx=obj.iter.currentIndex();
            coeff=obj.calculate_coeff(idx);
            mat= -coeff*kron(spin1.sx,spin2.sx)...
                -coeff*kron(spin1.sy,spin2.sy)...
                +2*coeff*kron(spin1.sz,spin2.sz);
        end
        
        function dataCell=data_cell(obj)
            nPair=obj.iter.getLength;
            nTerms=nPair*3;
            dataCell=cell(1, nTerms);
            
            obj.iter.setCursor(1);
            for ii=1:nPair
                spins=obj.iter.currentItem();
                spin_idx=obj.iter.currentIndex();
                spin1=spins{1}; spin2=spins{2};
                index1=spin_idx(1); index2=spin_idx(2);
                
                coeff=obj.calculate_coeff(spin_idx);
                
                dataCell{(ii-1)*2+1}={-coeff, obj.nbody, ...
                    index1, spin1.dim, reshape(spin1.sx, [], 1), ...
                    index2, spin2.dim, reshape(spin2.sx, [], 1)
                    };
                dataCell{(ii-1)*2+2}={-coeff, obj.nbody, ...
                    index1, spin1.dim, reshape(spin1.sy, [], 1), ...
                    index2, spin2.dim, reshape(spin2.sy, [], 1)
                    };
                dataCell{(ii-1)*2+3}={2*coeff, obj.nbody, ...
                    index1, spin1.dim, reshape(spin1.sz, [], 1), ...
                    index2, spin2.dim, reshape(spin2.sz, [], 1)
                    };
                obj.iter.nextItem;
            end
        end
    end
    
end

