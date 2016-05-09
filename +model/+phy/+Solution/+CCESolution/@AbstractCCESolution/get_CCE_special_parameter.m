function get_CCE_special_parameter(obj,p)
    switch obj.parameters.CCEStrategy
        case 'EnsembleCCE'

        case 'SingleSampleCCE'
            obj.parameters.seed=p.get_parameter('Dynamics','Seed');
        case 'DissipativeEnsembleCCE'
            obj.parameters.temperature=p.get_parameter('Condition', 'Temperature');
            obj.parameters.vertical_decay_rates=p.get_parameter('Dynamics','VerticalDecayRate');
            obj.parameters.parallel_decay_rates=p.get_parameter('Dynamics','ParallelDecayRate');
        case 'GFNCCESolution'
            obj.parameters.correlation_time=p.get_parameter('Dynamics','CorrelationTime');
        otherwise
            error('No such kind of CCE strategy.');
    end
end