%Changing strengths of synapses(weights) depending on time difference to produce output spikes.

clear all
close all


T = 500e-3;
N=100;
lambda = 1;
dt = 0.1e-3;
t=0:dt:T;

for i=1:N
    time_spike{i,1}(1)=exprnd(1/lambda);
end

for j=1:N
    i=1;
    while time_spike{j,1}(i)<=T
        time_spike{j,1}(i+1)=time_spike{j,1}(i)+exprnd(1/lambda);
        i=i+1;
    end
    time_spike{j,1}=time_spike{j,1}(1:i-1);
    clear i;
end
%AEF model physical parameters
C = 200e-12;
gL= 10e-9;
E_L= -70e-3;
V_T= -50e-3;
DelT = 2e-3;
a = 2e-9;
tau_t = 30e-3;
b = 0;
Vgamma= -58e-3;

%the state variables.
V = zeros(1,length(t));
U = zeros(1,length(t));

Io = 1e-12;
tau = 15e-3;
tau_s = tau/4;
wo =50;
sigma_w = 5;
we = sigma_w*randn(N,1) + wo;
delta_we = zeros(N,1);

Iapp = zeros(1,length(t));

spike_count=0;
iterations = 0;
while spike_count==0,    
   %% calculation of total input current.
   Iapp = zeros(1,length(t));
   for ii=1:N
        for i=1:length(t)
            for j=1:length(time_spike{ii})
                if t(i)<time_spike{ii}(j), break; end
                Iapp(i) = Iapp(i) + Io*we(ii)*(exp(-(t(i)-time_spike{ii}(j))/tau)-exp(-(t(i)-time_spike{ii}(j))/tau_s));
            end
        end
    end

    
    
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

    %% calculation of delta_t (Delta_tk)
    [y,x] = max(V);
    delta_t=zeros(N,1);
    for i=1:N
        for j=1:length(time_spike{i})
            if t(x)-time_spike{i}(j)<0, break;end
            delta_t(i)=t(x)-time_spike{i}(j);
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

    %% Plots
    figure
    
    subplot(3,1,1)
    plot(t,Iapp)
    xlabel('Time')
    ylabel('Input Current')

    subplot(3,1,2)
    plot(t,V)
    xlabel('Time')
    ylabel('Membrane Potential')
    
    subplot(3,1,3)
    plot(delta_t,delta_we)
    xlabel('Delta_t')
    ylabel('Delta_We')
    iterations = iterations + 1;
end
disp('iterations');
disp(iterations);
disp('weights')
disp(we)

