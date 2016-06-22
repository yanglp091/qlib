function Sx_val=GetObservableValue(obj,hamiltonian_list,NVDipIntStrength,initial_state)
    TimeList=obj.parameters.TimeList;
    ntime=length(TimeList);
    hami_num=length(hamiltonian_list);
    Vmat_list=cell(1,hami_num);
    Dmat_list=cell(1,hami_num);
    for kk=1:hami_num
       [V,D]=eig( full(hamiltonian_list{kk}) );
       Vmat_list{kk}=V;
       Dmat_list{kk}=D;
    end
    
    % The hamilonian_list is ordered like this: |00>, |11>, |10>, |01>
    % We need to calculate the trace of two list: 
    % exp(i*H_00 * t/2)*exp(i*H_11 * t/2)*exp(-i*H_01 * t/2)*exp(-i*H_10 * t/2) ...
    % exp(i*H_11 * t/2)*exp(i*H_00 * t/2)*exp(-i*H_10 * t/2)*exp(i*H_01 * t/2) ...
    idx_list1=[1,2,4,3];    idx_list2=[2,1,3,4]; 
    dim=size(initial_state,1);
    Sx_val=zeros(1,ntime);
    Sx_val(1)=2;
    for kk=2:ntime
        time=TimeList(kk)/2;
        Umat1=exp(1i*NVDipIntStrength*time);
        Umat2=exp(1i*NVDipIntStrength*time);
        for jj=1:hami_num
            V1=Vmat_list{ idx_list1(jj) };
            D1=Dmat_list{ idx_list1(jj) };            
            d_vec1=diag(D1);
            exp_d_vec_1=exp(-1i*time*d_vec1);
            Umat1=Umat1*V1*diag(exp_d_vec_1)*V1';
            
            V2=Vmat_list{ idx_list2(jj) };
            D2=Dmat_list{ idx_list2(jj) };            
            d_vec2=diag(D2);
            exp_d_vec_2=exp(-1i*time*d_vec2);
            Umat2=Umat2*V2*diag(exp_d_vec_2)*V2';
        end
        Sx_val(kk)=real( trace(Umat1) )/dim + real( trace(Umat2) )/dim;
    end      
end