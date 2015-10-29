%Response to a poisson stimulus by a AEF RS neuron.
clear all
close all
T = 500e-3;
lambda = 10;
dt = 0.1e-3;
t=0:dt:T;


time_spike(1)=exprnd(1/lambda);
i=1;
while time_spike(i)<=T
    time_spike(i+1)=time_spike(i)+exprnd(1/lambda);
    i=i+1;
end
time_spike=time_spike(1:i-1);
clear i;

Io = 1e-12;
we = 500;
tau = 15e-3;
tau_s = tau/4;

Iapp = zeros(1,length(t));
for i=1:length(t)
    for j=1:length(time_spike)
        if t(i)<time_spike(j), break; end
        Iapp(i) = Iapp(i) + Io*we*(exp(-(t(i)-time_spike(j))/tau)-exp(-(t(i)-time_spike(j))/tau_s));
    end
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
V(1)= -70e-3; 
U(1)= 0;


for i=1:length(t)-1;
    if V(i) >= 0,
        V(i)=0;
        V(i+1)= Vgamma;
        U(i+1) = U(i)+ b;
    else
        k = (1/C)*(-gL*(V(i)-E_L)+gL*DelT*exp((V(i)-V_T)/DelT)-U(i)+Iapp(i));
        l = (1/tau_t)*(a*(V(i)-E_L)-U(i));
        V(i+1)= V(i) + dt*k;
        U(i+1)= U(i) + dt*l;
    end   
end

ts=zeros(1,length(time_spike))
for v=2:length(time_spike)
    ts(v)=time_spike(v)-time_spike(v-1);
end

subplot(3,1,2)
plot(t,Iapp)
time_spike
xlabel('Time')
ylabel('Input Current')

subplot(3,1,3)
plot(t,V)
xlabel('Time')
ylabel('Membrane Potential')

subplot(3,1,1)
plot(ts)
xlabel('Index')
ylabel('Spike Activity')

length(ts)