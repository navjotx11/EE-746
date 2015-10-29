close all

spike_count=0;
iterations = 0;
while spike_count==0,
    Iapp = IappCalculation(we,S,t,N); %% input current.
    
    %% calculation of output voltage by solving diff. eqn.
    
    V(1)= -70e-3;
    U(1)= 0;
    
    spike_count=0;
    for i=1:length(t)-1;
        if V(i) >= 0,
            V(i)=0;
            V(i+1)= Vgamma;
            U(i+1) = U(i)+ b;
            spike_count=spike_count+1;
        else
            k = (1/C)*(-gL*(V(i)-E_L)+gL*DelT*exp((V(i)-V_T)/DelT)-U(i)+Iapp(i));
            l = (1/tau_t)*(a*(V(i)-E_L)-U(i));
            V(i+1)= V(i) + dt*k;
            U(i+1)= U(i) + dt*l;
        end
    end
    
    %% calculation of delta_t
    [y,x] = max(V);
    delta_t=zeros(N,1);
    for i=1:N
        for j=1:length(S{i})
            if t(x)-S{i}(j)<0, break;end
            delta_t(i)=t(x)-S{i}(j);
        end
    end
    %% calculation of new weight.
    
    gamma = 1;
    for i=1:N
        if spike_count~=0, break; end
        delta_we(i) = we(i)*gamma*(exp(-delta_t(i)/tau)-exp(-delta_t(i)/tau_s));
        we(i)= we(i)+ delta_we(i);
        if we(i)>=500, we(i)= 500;end
    end
    
    
    
    iterations = iterations + 1;
end

%% Plots
figure

subplot(2,1,1)
plot(t,Iapp)

subplot(2,1,2)
plot(t,V)

disp('iterations');
disp(iterations);
disp('weights')
disp(we)

