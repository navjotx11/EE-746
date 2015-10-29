
%Adaptive Exponential Integrate and Fire (AEF)
%Euler Method

%differential equation parameters.
clear
while true,
    type = input('Enter type of neuron(RS,IB or CH):  ','s');
 
    if strcmp(type,'RS'),
        C = 200e-12;
        gL= 10e-9;
        E_L= -70e-3;
        V_T= -50e-3;
        DelT = 2e-3;
        a = 2e-9;
        tau_t = 30e-3;
        b = 0;
        Vgamma= -58e-3;
        break;
    elseif strcmp(type,'IB'),
        C = 130e-12;
        gL= 18e-9;
        E_L= -58e-3;
        V_T= -50e-3;
        DelT = 2e-3;
        a = 4e-9;
        tau_t = 150e-3;
        b = 120e-12;
        Vgamma= -50e-3;
        break;
    elseif strcmp(type,'CH'),
        C = 200e-12;
        gL= 10e-9;
        E_L= -58e-3;
        V_T= -50e-3;
        DelT = 2e-3;
        a = 2e-9;
        tau_t = 120e-3;
        b = 100e-12;
        Vgamma= -46e-3;
        break;
    else disp('Invalid input!');
    end
end

h=0.1e-3; %time step
tf=500e-3; %the final time(stop time)
% ti=0:h/2:tf; %time for current I.
t =0:h:tf; %time for state variables V and U.

V = zeros(1,length(t));
U = zeros(1,length(t));
V(1)= -70e-3; %the state variables.
U(1)= 0;
I = 1e-12*input('Enter input current(in pA) - 250/350/450: '); %applied current.

for i=1:length(t)-1;
    if V(i) >= 0,
        V(i)=0;
        V(i+1)= Vgamma;
        U(i+1) = U(i)+ b;
    else
        k = (1/C)*(-gL*(V(i)-E_L)+gL*DelT*exp((V(i)-V_T)/DelT)-U(i)+I);
        l = (1/tau_t)*(a*(V(i)-E_L)-U(i));
        V(i+1)= V(i) + h*k;
        U(i+1)= U(i) + h*l;
    end
  
    
end

figure;
plot(t,V);
    