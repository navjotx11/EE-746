close all




delta_t = {};
spike_locations = 1;
iterations = 0;
while ~isempty(spike_locations),
    %% calculation of output voltage by solving diff. eqn.
    spike_locations = [];
    V(1)= -70e-3;
    U(1)= 0;
    S=S2;
    Iapp = IappCalculation(we,S,t,N); %% input current.
    spike_count=0;
    for i=1:length(t)-1;
        if V(i) >= 0,
            V(i)=0;
            V(i+1)= Vgamma;
            U(i+1) = U(i)+ b;
            spike_count=spike_count+1;
            spike_locations(spike_count)=t(i);
            
        else
            k = (1/C)*(-gL*(V(i)-E_L)+gL*DelT*exp((V(i)-V_T)/DelT)-U(i)+Iapp(i));
            l = (1/tau_t)*(a*(V(i)-E_L)-U(i));
            V(i+1)= V(i) + dt*k;
            U(i+1)= U(i) + dt*l;
        end
    end
    
    %% calculation of delta_t
    for i=1:N
        delta_t{i,1} = zeros(1,length(spike_locations));
    end
    
    for i=1:N
        for kk=1:length(spike_locations)
            for j=1:length(S{i})
                if spike_locations(kk)-S{i}(j)<0, break;end
                delta_t{i,1}(kk)=spike_locations(kk)-S{i}(j);
            end
        end
    end
    
    %% calculation of new weight.
    delta_we = zeros(N,1);
    gamma = 1;
    for i=1:N
        
        if isempty(spike_locations), break; end
        for j=1:length(spike_locations)
            delta_we(i) = delta_we(i) - we(i)*gamma*(exp(-delta_t{i}(j)/tau)-exp(-delta_t{i}(j)/tau_s));
        end
        we(i)= we(i)+ delta_we(i);
        if we(i)<=10, we(i)= 10;end
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

