classdef WavefunctionEvolution < model.phy.Dynamics.AbstractEvolutionKernel
    %WAVEFUNCTIONEVOLUTION Summary of this class goes here
    % This class is used to calculate the time evolutuion of the
    % wavefuncton in the Hilbert Spase
    
    properties
        timelist
        initial_state_type
        matrix_prefactor
        matrixList
        result
    end
    
    methods
        function obj=WavefunctionEvolution(qoperatorList,initial_state_type, prefactor)
              noperator=length(qoperatorList);               
              obj.matrixList=cell(1,noperator);
              for n=1:noperator
                   obj.matrixList{1,n}=full(qoperatorList{n}.getMatrix() );
              end
              
              obj.initial_state_type=initial_state_type;
              
              if nargin>2
                   obj.matrix_prefactor=prefactor;
              else
                   obj.matrix_prefactor=-1i;
              end             
              obj.result=[];
        end
          
        function state_out = calculate_evolution(obj, state_in, timelist)
          % This method is used to handle the piecewise evolution case, i.e.,
            % state_out(Nt)=exp(-i*H_N*t)*...*exp(-i*H_2*t)*exp(-i*H_1*t).
            % In this case, the final state of the current step can not be set
            % as the initial state of the next step. The evolution must
            % arranged as like this:
            % state_out(Nt)=exp(-i*H_N*t)*...*exp(-i*H_2t)*exp(-i*H_1*t).
               obj.timelist=timelist;
               ntime=length(timelist);
               dt=timelist(2)-timelist(1);
               obj.result=[obj.result, state_in];

               noperator=length(obj.matrixList);
               core_mat_list=cell(1,noperator);
               evolution_mat_list=cell(1,noperator);
               for m=1:noperator %initial the core matrice                  
                   core_mat_list{m}=expm(obj.matrix_prefactor(m)*dt*obj.matrixList{1,m});
                   evolution_mat_list{m}=1;
               end

               for n=2:ntime                   
                   for m=1:noperator                      
                     evolution_mat_list{1,m}=evolution_mat_list{1,m}*core_mat_list{1,m};
                   end
                   wave_vec=state_in;
                  for m=1:noperator
                       wave_vec=evolution_mat_list{1,m}*wave_vec;
                  end
                  state_out=wave_vec;
                  obj.result=[obj.result,state_out];
               end
        end
                   

        function  mean_val=mean_value(obj, obs_list,varargin)
              ntime=length(obj.timelist);
              len_obs=length(obs_list);
              mean_val=zeros(len_obs,ntime);
              
              for n=1:len_obs
                  mat=obs_list{n}.getMatrix;
                  mat_on_state=mat*obj.result;
                  if nargin>2 
                      % for single sample CCE, we need to calcualte <J|exp(-1i*H^{-}t)*exp(-1i*H^{+}t)|J> .
                      left_vec=varargin{1}{1};
                      mean_val(n,:)=left_vec*mat_on_state;
                  else

                    mean_val(1,:)=diag( (obj.result)' *mat_on_state );
                  end
              end
          end
          
    end
    
end

