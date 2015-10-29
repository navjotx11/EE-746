%Author%Author: Dharmashloka Debashis
%learning weights to differentiate between two different input stimuli.

S=S1;
Iapp = IappCalculation(we,S,t,N);

% spike_locations = 1;
% iterations = 0;

%% calculation of output voltage by solving diff. eqn.

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




%% Plots
figure

subplot(2,1,1)
plot(t,Iapp)

subplot(2,1,2)
plot(t,V)