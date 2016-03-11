classdef AbstractCCESolution < model.phy.Solution.AbstractSolution
    %ABSTRACTCCESOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods
        function obj=AbstractCCESolution(xml_file)
            obj@model.phy.Solution.AbstractSolution(xml_file);
        end
        function get_parameters(obj, p)
            get_parameters@model.phy.Solution.AbstractSolution(obj, p);
            obj.parameters.LoadSpinCollection=p.get_parameter('SpinCollection','LoadSpinCollection');
            if obj.parameters.LoadSpinCollection
                obj.parameters.SpinCollectionStrategy=p.get_parameter('SpinCollection', 'Source');
                obj.parameters.InputFile=p.get_parameter('SpinCollection', 'FileName');
                obj.BathSpinParameters(p);
                disp('spin collection loaded.');
            else
                obj.BathSpinParameters(p);
            end
                        
            obj.parameters.load_cluster_iter=p.get_parameter('Clustering','LoadCluterIterator');
            if obj.parameters.load_cluster_iter
                CluterIteratorName=p.get_parameter('Clustering','CluterIteratorName');
                data=load([OUTPUT_FILE_PATH,CluterIteratorName]);
                obj.keyVariables('cluster_iterator')=data.cluster_iterator;
                disp('cluster_iterator loaded.');
            else
                obj.parameters.CutOff=p.get_parameter('Clustering', 'CutOff');
                obj.parameters.MaxOrder=p.get_parameter('Clustering', 'MaxOrder');
            end

            obj.CentralSpinParameters(p);                      
            obj.parameters.MagneticField=p.get_parameter('Condition', 'MagneticField');
            obj.parameters.IsSecularApproximation=p.get_parameter('Interaction', 'IsSecular');
 
            obj.parameters.CCEStrategy=p.get_parameter('Dynamics','CCEStrategy');
            NTime=p.get_parameter('Dynamics', 'NTime');
            TMax=p.get_parameter('Dynamics', 'TMax');
            dt=TMax/(NTime-1);
            obj.parameters.TimeList=0:dt:TMax;
            obj.parameters.NTime=NTime;
            obj.parameters.TMax=TMax;
            obj.parameters.NPulse=p.get_parameter('Dynamics', 'NPulse');
            obj.parameters.seed=p.get_parameter('Dynamics','Seed');
        end
        function BathSpinParameters(obj, p)
             %set bath spin
             % get the parameters for setting the bath spins e.g. ZFS, eta, principle_axis, coordinate
                obj.parameters.SetBathSpins.SetSpin=p.get_parameter('SetBathSpins','SetSpin');
                if obj.parameters.SetBathSpins.SetSpin
                    SpeciesNumber=p.get_parameter('SetBathSpins','SpeciesNumber');
                    BathSpinsSettingCell=cell(1,4);
                    NameList=p.get_parameter('SetBathSpins','Name');
                    ZFSList=p.get_parameter('SetBathSpins','ZFS');
                    etaList=p.get_parameter('SetBathSpins','eta');
                    paxisList=p.get_parameter('SetBathSpins','principle_axis');
                    for k=1:SpeciesNumber
                        BathSpinsSettingCell{k}.name=NameList(k,:);
                        BathSpinsSettingCell{k}.ZFS=ZFSList(k,:);
                        if ~strcmp(etaList,'default')                   
                            BathSpinsSettingCell{k}.eta=etaList(k,:);
                        end
                        if ~strcmp(paxisList,'default')                   
                            BathSpinsSettingCell{k}.principle_axis=paxisList(k,:);
                        end                  
                    end
                    obj.parameters.SetBathSpins.BathSpinsSettingCell=BathSpinsSettingCell;               
                end
        end
        
         function CentralSpinParameters(obj,p)
             %set central spin
             % get the parameters for setting the central spin e.g. ZFS, eta, principle_axis, coordinate
                obj.parameters.SetCentralSpin.name=p.get_parameter('SetCentralSpin', 'Name');
                obj.parameters.SetCentralSpin.CentralSpinStates=p.get_parameter('SetCentralSpin', 'CentralSpinStates');
                obj.parameters.SetCentralSpin.SetSpin=p.get_parameter('SetCentralSpin', 'SetSpin');
                Orientation=p.get_parameter('SetCentralSpin', 'Orientation');
                Isotope=p.get_parameter('SetCentralSpin', 'Isotope');
                Coordinate=p.get_parameter('SetCentralSpin', 'Coordinate');
                if obj.parameters.SetCentralSpin.SetSpin
                    if ~strcmp(Orientation,'default') 
                        obj.parameters.SetCentralSpin.orientation=Orientation;
                    end
                    if ~strcmp(Isotope,'default') 
                        obj.parameters.SetCentralSpin.isotope=Isotope;
                    end
                    if ~strcmp(Coordinate,'default') 
                        obj.parameters.SetCentralSpin.coordinate=Coordinate;
                    end
                end
         end
         
         function cluster_iterator=generate_cluster_iterator(obj)
            import model.phy.SpinCollection.Strategy.FromFile
            import model.phy.SpinCollection.Strategy.FromSpinList
            import model.phy.SpinCollection.SpinCollection            
            import model.phy.SpinCollection.Iterator.ClusterIterator
            import model.phy.SpinCollection.Iterator.ClusterIteratorGen.CCE_Clustering
            
            if obj.parameters.load_cluster_iter             
                 cluster_iterator=obj.keyVariables('cluster_iterator');
                 if obj.parameters.SetBathSpins.SetSpin;
                     paraCell=obj.parameters.SetBathSpins.BathSpinsSettingCell;
                     cluster_iterator.spin_collection.set_spin(paraCell);
                 else
                     cluster_iterator.spin_collection.set_spin();
                 end
            else             
                  switch obj.parameters.SpinCollectionStrategy
                       case 'File'
                           spin_collection=SpinCollection( FromFile(...
                               [INPUT_FILE_PATH, obj.parameters.InputFile]));
                           if obj.parameters.SetBathSpins.SetSpin;
                               paraCell=obj.parameters.SetBathSpins.BathSpinsSettingCell;
                               spin_collection.set_spin(paraCell);
                           else
                               spin_collection.set_spin();
                           end
                       case 'SpinList'
                           error('not surported so far.');
                   end
                   obj.keyVariables('spin_collection')=spin_collection;

                   clu_para.cutoff=obj.parameters.CutOff;
                   clu_para.max_order=obj.parameters.MaxOrder;
                   disp('clustering begins...')
                   tic
                   % the cce strategy can be change
                   cce=CCE_Clustering(spin_collection, clu_para);

                   cluster_iterator=ClusterIterator(spin_collection,cce);
                   [~]=cluster_iterator.cross_relation_gen();
                   obj.keyVariables('cluster_iterator')=cluster_iterator;
                   save([OUTPUT_FILE_PATH, 'cluster_iterator', obj.timeTag, '.mat'],'cluster_iterator');
                   toc
             end
         end
         
         function calculate_total_coherence(obj, evolu_para,clst_para,cluster_iter) 
%%           evolution parameters  
%            evolu_para.npulse:: the pulse number
%            evolu_para.is_secular:: determine whether the secular approximation is taken or not
%            evolu_para.MagneticField:: the magnetic field vector             
%            evolu_para.center_spin_states:: the states of the center spin involved in this calculatioin
%            evolu_para.strategy_name:: CCE strategy name
%            evolu_para.timelist:: the evolution time list

%         cluster parameters
%           clst_para.bath_spin_collection:: a SpinCollection, including the  all bath spins
%           clst_para.center_spin:: a Spin
%           clst_para.bath_spin_state:: a list gives out the specific configuration of the bath
%                                       spin state for SingleSampleCCE,
%                                       e.g. [1,1,1,1,1,1,1,...,1] for the ground state.                                      


           ncluster=cluster_iter.getLength;
           ntime=length(evolu_para.timelist);           
           cluster_index_list=cluster_iter.index_list;
           MagneticField=evolu_para.MagneticField;
           strategy_name=evolu_para.strategy_name;

           CoherenceMatrix=zeros(ncluster,ntime);
           disp('calculate the cluster-coherence matrix ...');
           tic
           
          %In order to eliminate the parfor warning, I have to arrange the
          %swith...case... in this form. Beside, I try to add a field "clst_index" in the structure data "clst_para" in every loop first. 
          % But this action is forbidden in parfor circulation. So I have to change the way to construct 
          % the AbstractClusterCoherence class. This is pretty ulgy, but I have to do this. 
           switch strategy_name
              case 'EnsembleCCE'
                      parpool();
                      parfor n=1:ncluster 
                          Condition=model.phy.LabCondition.getCondition;              
                          Condition.setValue('magnetic_field',MagneticField);

                          %calculate cluster coherence              
                          clst_index=cluster_index_list{n,1};
                          import model.phy.Solution.CCESolution.CCECoherenceStrategy.ECCEClusterCoherence 
                          clst_coh=ECCEClusterCoherence(clst_index,clst_para);
                          CoherenceMatrix(n,:)=clst_coh.calculate_cluster_coherence(evolu_para);
                          delete(clst_coh);
                      end
                      delete(gcp('nocreate'));
                                       
              case 'SingleSampleCCE'
%                       parpool();
                      for n=1:ncluster 
%                           disp(['calculating the ' num2str(n) 'th cluster coherence....']);
                          Condition=model.phy.LabCondition.getCondition;              
                          Condition.setValue('magnetic_field',MagneticField);

                          %calculate cluster coherence              
                          clst_index=cluster_index_list{n,1};
                          import model.phy.Solution.CCESolution.CCECoherenceStrategy.SSCCEClusterCoherence
                          clst_coh=SSCCEClusterCoherence(clst_index,clst_para);
                          CoherenceMatrix(n,:)=clst_coh.calculate_cluster_coherence(evolu_para);
                          delete(clst_coh);
                      end
%                       delete(gcp('nocreate'));
               otherwise
                       error('No such CCE method.......');              
            end


           toc
           disp('calculation of the cluster-coherence matrix finished.');          

           obj.calculate_tilde_coh_matrix(CoherenceMatrix,cluster_iter);
            
           if ncluster<20000
                obj.keyVariables('coherence_matrix')=CoherenceMatrix;
           else
                timeTag=datestr(clock,'yyyymmdd_HHMMSS');
                save([OUTPUT_FILE_PATH, 'coherence_matrix', timeTag, '.mat'],'CoherenceMatrix');
                clear CoherenceMatrix;
           end
         end
        function calculate_tilde_coh_matrix(obj,cohmat,cluster_iter)
            subcluster_list=cluster_iter.cluster_info.subcluster_list;
            cluster_number_list=[cluster_iter.cluster_info.cluster_number_list{:,2}];
            
            ncluster=length(subcluster_list);
            ntime=length(cohmat(1,:));
            coh_tilde_mat=zeros(ncluster,ntime);
            coh_total=ones(1,ntime);
            coh=struct();
            
            cceorder=1;
            endpoints=cumsum(cluster_number_list);
            for m=1:ncluster
                subcluster=subcluster_list{m};
                nsubcluster=length(subcluster);
                coh_tilde=cohmat(m,:);
                if nsubcluster==0
                    coh_tilde_mat(m,:)= coh_tilde;
                elseif nsubcluster>0
                    for n=1:nsubcluster;
                        coh_tilde_sub=coh_tilde_mat(subcluster(n),:);
                        coh_tilde=coh_tilde./coh_tilde_sub;
                    end
                    coh_tilde_mat(m,:)=coh_tilde;
                end

                coh_total=coh_total.*coh_tilde;
                if m==endpoints(1,cceorder)
                    field_name=strcat('coherence_cce_',num2str(cceorder));
                    coh.(field_name)=coh_total;
                    cceorder=cceorder+1;
                end
            end
            coh.('coherence')= coh_total;            
           if ncluster<20000          
               obj.keyVariables('coherence_tilde_matrix')=coh_tilde_mat;
           else
               timeTag=datestr(clock,'yyyymmdd_HHMMSS');
               save([OUTPUT_FILE_PATH, 'coherence_tilde_matrix', timeTag, '.mat'],'coh_tilde_mat');
               clear coh_tilde_mat;
           end
           coh.('timelist')=obj.parameters.TimeList;
            obj.keyVariables('coherence')=coh;
        end
    end
    
end

