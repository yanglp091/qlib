classdef MagRSolution < model.phy.Solution.AbstractSolution
    %MAGRSOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj=MagRSolution(xml_file)
            obj@model.phy.Solution.AbstractSolution(xml_file);
        end
        
        function get_parameters(obj, p)
            get_parameters@model.phy.Solution.AbstractSolution(obj, p);
            
            obj.parameters.MagneticField=p.get_parameter('Condition', 'MagneticField');
            
            obj.parameters.SpinCollectionStrategy=p.get_parameter('SpinCollection', 'Source');
            obj.parameters.InputFile=p.get_parameter('SpinCollection', 'FileName');
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
            Condition=model.phy.LabCondition.getCondition;              
            Condition.setValue('magnetic_field',obj.parameters.MagneticField);
            
           [evolution_para,cluster_para,clst_index]=obj.pre_calculation();
           
           clst_coh=model.phy.Solution.CCESolution.CCECoherenceStrategy.UnCorrelatedClusterCoherence(clst_index,cluster_para);
           coh=clst_coh.calculate_cluster_coherence(evolution_para);
           obj.keyVariables('coherence')=coh;
           obj.keyVariables('timelist')=obj.parameters.TimeList;
           disp('Calculation of this solution finishes.');
           
        end
        
       function [evolution_parameter,cluster_parameter,clst_index]=pre_calculation(obj)
           import model.phy.PhysicalObject.NV

           evolution_parameter.center_spin_states=obj.parameters.SetCentralSpin.CentralSpinStates;
           evolution_parameter.timelist=obj.parameters.TimeList;
           evolution_parameter.npulse=obj.parameters.NPulse;
           evolution_parameter.is_secular=obj.parameters.IsSecularApproximation;
           
           evolution_parameter.AddContact=obj.parameters.AddContact;
           evolution_parameter.contactInt=obj.parameters.contactInt;
           evolution_parameter.AddDipInt=obj.parameters.AddDipInt;
           
           

           evolution_parameter.AddDecay=obj.parameters.AddDecay;           
           evolution_parameter.vertical_decay_rates=obj.parameters.VerticalDecayRate;
           evolution_parameter.parallel_decay_rates=obj.parameters.ParallelDecayRate;           
                      
           center_spin_name=obj.parameters.SetCentralSpin.name;
           para_central_spin=obj.parameters.SetCentralSpin; 
           center_spin=eval(strcat(center_spin_name,'(','para_central_spin',')'));
           obj.keyVariables('center_spin')=center_spin;
           
           cluster_parameter.center_spin=center_spin.espin;
           cluster_parameter.bath_spin_collection=obj.keyVariables('spin_collection');
           nspin=obj.keyVariables('spin_collection').getLength;
           clst_index=1:nspin;
        end
    end
    
end

