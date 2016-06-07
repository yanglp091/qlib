classdef ProbeMagR < model.phy.Solution.AbstractSolution
    %MAGRSOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj=ProbeMagR(xml_file)
            obj@model.phy.Solution.AbstractSolution(xml_file);
        end
        
        function get_parameters(obj, p)
            get_parameters@model.phy.Solution.AbstractSolution(obj, p);
            
            obj.parameters.MagneticField=p.get_parameter('Condition', 'MagneticField');
            Condition=model.phy.LabCondition.getCondition;              
            Condition.setValue('magnetic_field',obj.parameters.MagneticField);
            
            obj.BathSpinParameters(p);
            obj.CentralSpinParameters(p);
            
            obj.parameters.AddContact = p.get_parameter('Interaction', 'AddContact');
            obj.parameters.contactInt = p.get_parameter('Interaction', 'ContactInt');
            obj.parameters.AddDipInt = p.get_parameter('Interaction', 'AddDipInt');  
            obj.parameters.IsSecularApproximation=p.get_parameter('Interaction', 'IsSecular');

            NTime=p.get_parameter('Dynamics', 'NTime');
            TMax=p.get_parameter('Dynamics', 'TMax');
            dt=TMax/(NTime-1);
            obj.parameters.TimeList=0:dt:TMax;
            obj.parameters.NTime=NTime;
            obj.parameters.TMax=TMax;
            obj.parameters.NPulse=p.get_parameter('Dynamics', 'NPulse');
            
            obj.parameters.AddDecay = p.get_parameter('Dynamics', 'AddDecay');
            obj.parameters.VerticalDecayRate = p.get_parameter('Dynamics', 'VerticalDecayRate');
            obj.parameters.ParallelDecayRate = p.get_parameter('Dynamics', 'ParallelDecayRate');                      
        end
               
        function perform(obj)
            nseed=obj.parameters.nseed;
            spin_list_tot=obj.keyVariables('SpinListTotal');
            spin_collection_cell=cell(1,nseed);
            
            time=obj.parameters.TimeList;
            ntime=length(time);
            coherence_matrix=zeros(nseed,ntime);
            dip_info=cell(nseed,1);
            pos=1;
            for kk=1:nseed
               idx=pos:1:pos+3;
               spin_list=spin_list_tot(idx);
               [evolution_para,cluster_para,clst_index,spin_collection]=obj.pre_calculation(spin_list);
               spin_collection_cell{kk}=spin_collection;               
               clst_coh=model.phy.Solution.CCESolution.CCECoherenceStrategy.UnCorrelatedClusterCoherence(clst_index,cluster_para);
               coh=clst_coh.calculate_cluster_coherence(evolution_para);
               coherence_matrix(kk,:)=coh;
               
               [obs_tf,obs_pro]=obj.get_cluster_transition_spectrum(clst_coh,evolution_para.IntPara);
               dip_info{kk}.frequency=obs_tf;
               dip_info{kk}.amplitude=obs_pro;
               pos=pos+4;
            end
            if nseed>1
               coherence=sum(coherence_matrix)/nseed;
            else
               coherence=coherence_matrix;
            end
           obj.keyVariables('coherence_matrix')=coherence_matrix;
           obj.keyVariables('coherence')=coherence;
           obj.keyVariables('timelist')=time;
           obj.keyVariables('spin_collection_cell')=spin_collection_cell;
          
           disp('Calculation of this solution finishes.');
           obj.keyVariables('dip_info')=dip_info;
        end
        
       function [evolution_parameter,cluster_parameter,clst_index,spin_collection]=pre_calculation(obj,spin_list)
           import model.phy.PhysicalObject.NV
           import model.phy.SpinCollection.SpinCollection
           import model.phy.SpinCollection.Strategy.FromSpinList
           
           spin_collection=SpinCollection(FromSpinList(spin_list));
           evolution_parameter.center_spin_states=obj.parameters.SetCentralSpin.CentralSpinStates;
           evolution_parameter.timelist=obj.parameters.TimeList;
           evolution_parameter.npulse=obj.parameters.NPulse;
           evolution_parameter.is_secular=obj.parameters.IsSecularApproximation;

           IntPara.AddDipInt=obj.parameters.AddDipInt;          
           IntPara.AddContact=obj.parameters.AddContact;
           IntPara.interaction=obj.parameters.contactInt;
           IntPara.AddIndexList=1;
           idx=cell(6,1);idx{1}=[1,2];idx{2}=[2,3];idx{3}=[3,4];idx{4}=[4,1];idx{5}=[1,3];idx{6}=[2,4];
           IntPara.IndexList=idx;
           evolution_parameter.IntPara=IntPara;
           
           

           evolution_parameter.AddDecay=obj.parameters.AddDecay;           
           evolution_parameter.vertical_decay_rates=obj.parameters.VerticalDecayRate;
           evolution_parameter.parallel_decay_rates=obj.parameters.ParallelDecayRate;           
                      
           center_spin_name=obj.parameters.SetCentralSpin.name;
           para_central_spin=obj.parameters.SetCentralSpin; 
           center_spin=eval(strcat(center_spin_name,'(','para_central_spin',')'));
           obj.keyVariables('center_spin')=center_spin;
           
           cluster_parameter.center_spin=center_spin.espin;
           cluster_parameter.bath_spin_collection=spin_collection;
           nspin=spin_collection.getLength;
           clst_index=1:nspin;
        end
    end
    
end

