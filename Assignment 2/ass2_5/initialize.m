%poisson stimuli S1 and S2 and random weights.

clear all
close all

T = 500e-3;
N=100;
lambda = 1;
dt = 0.1e-3;
t=0:dt:T;


% input spike times allocations (poisson stimuli) (the next 4 for loops are used for this.)
for i=1:N
    S1{i,1}(1)=exprnd(1/lambda);
end

for j=1:N
    i=1;
    while S1{j,1}(i)<=T
        S1{j,1}(i+1)=S1{j,1}(i)+exprnd(1/lambda);
        i=i+1;
    end
    S1{j,1}=S1{j,1}(1:i-1);
    clear i;
end

for i=1:N
    S2{i,1}(1)=exprnd(1/lambda);
end

for j=1:N
    i=1;
    while S2{j,1}(i)<=T
        S2{j,1}(i+1)=S2{j,1}(i)+exprnd(1/lambda);
        i=i+1;
    end
    S2{j,1}=S2{j,1}(1:i-1);
    clear i;
end



wo =200;
sigma_w = 20;
we = sigma_w*randn(N,1) + wo;
tau = 15e-3;
tau_s = tau/4;

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

