
%Hodgkin-Huxley neuron model
%Euler method

clear
close all
%The physical parameters.
milli= 1; 
C=1*milli^2;
ENa=50*milli;
EK=-77*milli;
El=-55*milli;
gNa=120;
gK=36;
gl=0.3;


%anonymous functions for an,am,ah, etc.
an = @(x) 0.01*(x/milli+55)/(1-exp(-0.1*(x/milli+55)));
am = @(x) 0.1*(x/milli+40)/(1-exp(-0.1*(x/milli+40)));
ah = @(x) 0.07*exp(-0.05*(x/milli+65));
bn = @(x) 0.125*exp(-0.0125*(x/milli+65));
bm = @(x) 4*exp(-0.0556*(x/milli+65));
bh = @(x) 1/(1+exp(-0.1*(x/milli+35)));
fx = @(a,b,x) a*(1-x) - b*x; %function for the RHS of the diff. eqn.(11) 

%time variables;
dt=0.01*milli; %time step;
tf = 150*milli; %final time(stop time)
t = 0:dt:tf;

%state variables
V = zeros(1,length(t));
n = zeros(1,length(t));
m = zeros(1,length(t));
h = zeros(1,length(t));

iNa = zeros(1,length(t));
iK = zeros(1,length(t));
il = zeros(1,length(t));

%initializations
V(1) = -65*milli;
n(1) = 0.3159;
m(1) = 0.0515;
h(1) = 0.6017;
iNa(1)= gNa*m(1)^3*h(1)*(V(1)-ENa);
iK(1) = gK*n(1)^4*(V(1)-EK);
il(1) = gl*(V(1)-El);

I = zeros(1,length(t));
Io = 15;
for i=1:length(t)
    if i>=2*30*milli/dt && i<3*30*milli/dt, I(i) = Io;%applied current.
    else I(i) = 0;
    end
end

for i=1:length(t)-1
    V(i+1) = V(i) + dt*(I(i)-iNa(i)-iK(i)-il(i))/C;
    
    n(i+1) = n(i) + dt*fx(an(V(i)),bn(V(i)),n(i));
    m(i+1) = m(i) + dt*fx(am(V(i)),bm(V(i)),m(i));
    h(i+1) = h(i) + dt*fx(ah(V(i)),bh(V(i)),h(i));
    
    iNa(i+1)= gNa*m(i+1)^3*h(i+1)*(V(i+1)-ENa);
    iK(i+1) = gK*n(i+1)^4*(V(i+1)-EK);
    il(i+1) = gl*(V(i+1)-El);
end
%instantaneous power.
PNa = iNa.*(V-ENa);
PK = iK.*(V-EK);
Pl = il.*(V-El);
Pc = V.*(I-iNa-iK-il);

Ec=trapz(t(find(t==60):find(t==90)), Pc(find(t==60):find(t==90)));
% ENa=trapz(t(find(t==60):find(t==90)), PNa(find(t==60):find(t==90)));
% EK=trapz(t(find(t==60):find(t==90)), PK(find(t==60):find(t==90)));
% El=trapz(t(find(t==60):find(t==90)), Pl(find(t==60):find(t==90)));
Ec/3

plot(t,V);
figure
plot(t,m,t,n,t,h);
figure
plot(t,PNa,t,PK,t,Pl,t,Pc);