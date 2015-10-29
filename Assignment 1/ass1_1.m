%Differencial equation is solved using Runge-Kutta 2nd order Heuns method.

clear;
%parameters of differential equation.
C=300e-12;
gL=30e-9;
V_T=20e-3;
E_L=-70e-3;

N=10; % number of neurons.
Vo = zeros(N,1);
for i=1:N       %initializing the membrane voltages for the N neurons.
    Vo(i,1) = E_L;
end
tf = 500e-3;   % final time(stop time).
h=0.1e-3;     % time step(delta t).
t=0:h:tf;
V = zeros(N,length(t));

%current parameters
I = zeros(N,length(t));
alfa=0.1;
Ic=3e-9;
for i=1:N
    I(i,:)= (1+i*alfa)*Ic; %input current.
end


V(:,1)=Vo;
k1 = zeros(N,1);
k2 = zeros(N,1);

index=zeros(N,1);

count=zeros(N,1);
for i=1:length(t)-1
    %V(:,i+1)=V(:,i)+0.5*h*(k1+k2);
    for j=1:N
        if V(j,i)>= V_T, 
            V(j,i+1)=E_L;
            if count(j)==0, 
                index(j)=i;
                count(j)=count(j)+1;
            end
        else
            k1=-gL*V(:,i)/C + gL*E_L/C + I(:,i)/C;
            k2=-gL*(V(:,i)+ k1*h)/C + gL*E_L/C + I(:,i+1)/C;
            V(j,i+1) = V(j,i)+0.5*h*(k1(j)+k2(j));
        end
    end
end
clear count;

% V(1,length(t))
% V(10,length(t))



%spike time vs input current.

t_spike = zeros(N,1);


for i=1:N
    t_spike(i,1)=t(index(i));
end
% disp(I(:,1));
% disp(t_spike);


subplot(2,1,1)
plot(t,V(2,:),t,V(4,:),t,V(6,:),t,V(8,:))
title('V vs Time')
ylabel('Voltage')
xlabel('Time')

subplot(2,1,2)
plot(I(:,1),t_spike)
title('Spike Time vs Input Current')
ylabel('Spike Time')
xlabel('Input Current')