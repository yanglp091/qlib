classdef BioCompassSimulation  < model.phy.Solution.AbstractSolution
    %BIOCOMPASSSIMULATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj=BioCompassSimulation(xml_file)
            obj@model.phy.Solution.AbstractSolution(xml_file);
        end
        
        function get_parameters(obj, p)
            import model.phy.SpinCollection.Strategy.FromFile
            import model.phy.SpinCollection.Strategy.FromSpinList
            import model.phy.SpinCollection.SpinCollection  

            get_parameters@model.phy.Solution.AbstractSolution(obj, p);

%             obj.parameters.MagneticField=p.get_parameter('Condition', 'MagneticField');
%             Condition=model.phy.LabCondition.getCondition;              
%             Condition.setValue('magnetic_field',obj.parameters.MagneticField);

            obj.SetCentralSpins(p);

            obj.parameters.InputFile=p.get_parameter('SetMagRSpin', 'FileName');    
            MagR_sc=SpinCollection( FromFile([INPUT_FILE_PATH, obj.parameters.InputFile]));
            MagR_sc.set_spin();
            obj.parameters.MagRSpinCollection=MagR_sc;
            obj.parameters.MagRSpinNumber=MagR_sc.getLength;
            obj.parameters.InitialStateSeed=p.get_parameter('SetMagRSpin', 'InitialStateSeed');

            obj.parameters.AddContact = p.get_parameter('Interaction', 'AddContact');
            obj.parameters.contactInt = p.get_parameter('Interaction', 'ContactInt');
            obj.parameters.AddDipInt = p.get_parameter('Interaction', 'AddDipInt');  

            NTime=p.get_parameter('Dynamics', 'NTime');
            TMax=p.get_parameter('Dynamics', 'TMax');
            dt=TMax/(NTime-1);
            obj.parameters.TimeList=0:dt:TMax;
            obj.parameters.NTime=NTime;
            obj.parameters.TMax=TMax;
        %     obj.parameters.NPulse=p.get_parameter('Dynamics', 'NPulse');
        % 
        %     obj.parameters.AddDecay = p.get_parameter('Dynamics', 'AddDecay');
        %     obj.parameters.VerticalDecayRate = p.get_parameter('Dynamics', 'VerticalDecayRate');
        %     obj.parameters.ParallelDecayRate = p.get_parameter('Dynamics', 'ParallelDecayRate');                      
        end
        function mean_values=perform(obj)
            [initial_state,spin_collection] = obj.GetInitialState;
            [hamiltonian_list,NVDipIntStrength] = obj.GetHamiltonian(spin_collection);
%             observables = obj.GetObservables(spin_collection);
%             mean_values = obj.GetSingletProbability(hamiltonian_list,initial_state);
            mean_values=obj.GetObservableValue(hamiltonian_list,NVDipIntStrength,initial_state);
            
            obj.StoreKeyVariables(mean_values);
            
            
        end

    end
    
end

