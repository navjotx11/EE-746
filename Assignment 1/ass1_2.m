
%Izhikevich Model
%4th Order Runga-Kutta Method

%differential equation parameters.
clear
while true,
    type = input('Enter type of neuron(RS,IB or CH):  ','s');
   
    if strcmp(type,'RS'),
        C = 100e-12;
        kz= 0.7e-6;
        Er= -60e-3;
        Et= -40e-3;
        a = 0.03e3;
        b = -2e-9;
        c = -50e-3;
        d = 100e-12;
        Vpeak= 35e-3;
        break;
    elseif strcmp(type,'IB'),
        C = 150e-12;
        kz= 1.2e-6;
        Er= -75e-3;
        Et= -45e-3;
        a = 0.01e3;
        b = 5e-9;
        c = -56e-3;
        d = 130e-12;
        Vpeak= 50e-3;
        break;
    elseif strcmp(type,'CH'),
        C = 50e-12;
        kz= 1.5e-6;
        Er= -60e-3;
        Et= -40e-3;
        a = 0.03e3;
        b = 1e-9;
        c = -40e-3;
        d = 150e-12;
        Vpeak= 25e-3;
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
V(1)= c; %the state variables.
U(1)= d;
I = 1e-12*input('Enter input current(in pA) - 400/500/600: '); %applied current.

for i=1:length(t)-1;
   if V(i)< Vpeak, 
        k1 = ((V(i)-Er)*(V(i)-Et)*kz - U(i) + I)/C;
        l1 = a*b*(V(i)-Er) - a*U(i);
        k2 = ((V(i)+ 0.5*h*k1-Er)*(V(i)+0.5*h*k1-Et)*kz - (U(i)+0.5*h*l1) + I)/C;
        l2 = a*b*(V(i)+0.5*h*k1-Er) - a*(U(i)+0.5*h*l1);
        k3 = ((V(i)+ 0.5*h*k2-Er)*(V(i)+0.5*h*k2-Et)*kz - (U(i)+0.5*h*l2) + I)/C;
        l3 = a*b*(V(i)+0.5*h*k2-Er) - a*(U(i)+0.5*h*l2);
        k4 = ((V(i)+ h*k3-Er)*(V(i)+h*k3-Et)*kz - (U(i)+ h*l3) + I)/C;
        l4 = a*b*(V(i)+ h*k3-Er) - a*(U(i)+ h*l3);

        k = (k1+2*k2+2*k3+k4)/6;
        l = (l1+2*l2+2*l3+l4)/6;

        V(i+1)= V(i) + h*k;
        U(i+1)= U(i) + h*l;
        
   else
        V(i+1)= c;
        U(i+1) = U(i)+d;
    end
end

figure;
plot(t,V);
    