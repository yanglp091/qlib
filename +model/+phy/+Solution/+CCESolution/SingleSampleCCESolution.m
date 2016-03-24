classdef SingleSampleCCESolution < model.phy.Solution.CCESolution.AbstractCCESolution
    %ENSEMBLECCESOLUTION Summary of this class goes here
    %   EnsembleCCESolution needs the following input paramters:
    %   1. parameters.SpinCollectionStrategy
    %   2. parameters.InputFile
    %   3. parameters.SetBathSpins
    %   4. parameters.SetCentralSpin 
    %   5. parameters.MagneticField
    %   6.parameters.CutOff
    %   7.parameters.MaxOrder    
    %   9. parameters.IsSecularApproximation
    %   6. parameters.NPulse
    %   10. parameters.NTime
    %   11. parameters.TMax
    %   12. parameters.TimeList

   
    properties
    end
    
    methods
        function obj=SingleSampleCCESolution(xml_file)
            obj@model.phy.Solution.CCESolution.AbstractCCESolution(xml_file);
        end
        
         function perform(obj)
           Condition=model.phy.LabCondition.getCondition;              
           Condition.setValue('magnetic_field',obj.parameters.MagneticField);  
             
           %%  Generate Spin Collection FromSpinList and generate clusters 
           cluster_iterator=obj.generate_cluster_iterator();
           [evolution_para,cluster_para]=obj.pre_calculation(cluster_iterator);
           obj.calculate_total_coherence(evolution_para,cluster_para,cluster_iterator); 
           disp('Calculation of this solution finishes.');
        end
        
        function [evolution_parameter,cluster_parameter]=pre_calculation(obj,cluster_iterator)
           import model.phy.PhysicalObject.NV

           evolution_parameter.center_spin_states=obj.parameters.SetCentralSpin.CentralSpinStates;
           evolution_parameter.timelist=obj.parameters.TimeList;
           evolution_parameter.npulse=obj.parameters.NPulse;
           evolution_parameter.is_secular=obj.parameters.IsSecularApproximation;
           evolution_parameter.MagneticField=obj.parameters.MagneticField; 
           evolution_parameter.strategy_name=obj.parameters.CCEStrategy;           
                      
           center_spin_name=obj.parameters.SetCentralSpin.name;
           para_central_spin=obj.parameters.SetCentralSpin; 
           center_spin=eval(strcat(center_spin_name,'(','para_central_spin',')'));
           obj.keyVariables('center_spin')=center_spin;
           
           cluster_parameter.center_spin=center_spin.espin;
           cluster_parameter.bath_spin_collection=cluster_iterator.spin_collection;
           cluster_parameter.bath_spin_state=obj.generate_bath_spin_state(cluster_iterator);
        end
        function bs_state=generate_bath_spin_state(obj,cluster_iterator)
            seed=obj.parameters.seed;
            nspin=cluster_iterator.spin_collection.getLength;
            dim_list=cluster_iterator.spin_collection.getDimList;
            rng(seed);rand_numbers=randi([1,100],1,nspin);
            
            bs_state=zeros(1,nspin);
            for kk=1:nspin
               dim=dim_list(kk); 
               bs_state(1,kk)=mod(rand_numbers(kk),dim)+1; 
            end
        end
        
    end
    
end

