%Synaptic Current Exchange in Neurons
%Leaky Integrate and Fire with Refractory Period

N=5;

for i=1:N
spike_time{i}=[ ];
arrival_time{i}=[ ];
strength{i}=[ ];
pre_neuron{i}=[ ];
end

C=300e-12;
gL=30e-9;
V_T=20e-3;
E_L=-70e-3;
Rp= 2e-3;
Ip=50e-9;
Io=1e-12;
tau=15e-3;
tau_s=tau/4;
w=3000;

tau_d=zeros(3,2);   %Axonal Delays  1->2,2->3,3->4;1->1,2->5;
tau_d(1,1)=1e-3;
tau_d(1,2)=8e-3;
tau_d(2,1)=5e-3;
tau_d(2,2)=5e-3;
tau_d(3,1)=9e-3;
tau_d(3,2)=1e-3;

%Specifying Netwok Configuration
for i=1:N
    Fanout{i}=[];
    Weights{i}=[];
    Delay{i}=[];
end
for i=2:4
  Fanout{i}=[1 , 5];
  Weights{i}=[w,w];
end

Delay{2}=[tau_d(1,1) , tau_d(1,2)];
Delay{3}=[tau_d(2,1) , tau_d(2,2)];
Delay{4}=[tau_d(3,1) , tau_d(3,2)];

% Timing Variables
tf = 500e-3;   % final time(stop time).
h=0.1e-3;     % time step(delta t).
t=0:h:tf;

Iapp=zeros(N,length(t));
Isyn=zeros(N,length(t));
V=zeros(N,length(t));

for i=1:length(t)
    if t(i)>=0e-3 && t(i)<1e-3,
        Iapp(2,i)=Ip;
    end
    if t(i)>=4e-3 && t(i)<5e-3,
        Iapp(3,i)=Ip;
    end
    if t(i)>=8e-3 && t(i)<9e-3,
        Iapp(4,i)=Ip;
    end
end
count=0;
V(:,1)=E_L;

for i=1:N
fanout{i}=[];
delay{i}=[];
weight{i}=[];
end


%temp_c=0;
k1=zeros(3,1);
k2=zeros(3,1);

for j=2:4
    for k=1:length(t)-1
        if V(j,k)>= V_T, 
            V(j,k+1)=E_L;
            
            spike_time{j}=[spike_time{j},t(k)];count=20;
            temp=Fanout{j};
            for c=1:length(temp)
                arrival_time{temp(c)}=[arrival_time{temp(c)},t(k)+tau_d(j-1,rem(temp(c),3))];
                pre_neuron{temp(c)}=[pre_neuron{temp(c)},j];
                strength{temp(c)}=[strength{temp(c)}, w];
            end
        else
            if count~=0
                V(j,k+1)=E_L;
                count=count-1;
            else
            k1=-gL*V(j,k)/C + gL*E_L/C + Iapp(j,k)/C + Isyn(j,k)/C;
            k2=-gL*(V(j,k)+ k1*h)/C + gL*E_L/C + Iapp(j,k+1)/C +Isyn(j,k+1)/C;
            V(j,k+1) = V(j,k)+0.5*h*(k1+k2);
            end
        end
    end
end

for i=1:2
    I_temp=zeros(1,length(t));
    for j=2:4
        temp_sn=1;
        while temp_sn < length(spike_time{j})
                           
            for k=round(spike_time{j}(temp_sn)/h)+ round(tau_d(j-1,i)/h):length(t)%round(spike_time{j}(temp_sn+1)/h)
                I_temp(k)=I_temp(k)+(Io*w*(exp(-(t(k)-spike_time{j}(temp_sn))/tau) - exp(-(t(k)-spike_time{j}(temp_sn))/tau_s)));
            end
               
            temp_sn=temp_sn+1;
            
        end
        
         for k=round(spike_time{j}(temp_sn)/h) + round(tau_d(j-1,i)/h):length(t)
             I_temp(k)=I_temp(k)+(Io*w*(exp(-(t(k)-spike_time{j}(temp_sn))/tau) - exp(-(t(k)-spike_time{j}(temp_sn))/tau_s)));
         end
        
    end
   Isyn((4*i)-3,:)=Isyn((4*i)-3,:)+I_temp;
end
       
for j=1:4:5
    for k=1:length(t)-1
        if V(j,k)>= V_T, 
            V(j,k+1)=E_L;
            spike_time{j}=[spike_time{j},t(k)];count=20;
            temp=Fanout{j};
            for c=1:length(temp)
                arrival_time{temp(c)}=[arrival_time{temp(c)},t(k)+tau_d(j-1,rem(temp(c),3))];
                pre_neuron{temp(c)}=[pre_neuron{temp(c)},j];
                strength{temp(c)}=[strength{temp(c)}, w];
            end
        else
            if count~=0
                V(j,k+1)=E_L;
                count=count-1;
            else
            k1=-gL*V(j,k)/C + gL*E_L/C + Iapp(j,k)/C + Isyn(j,k)/C;
            k2=-gL*(V(j,k)+ k1*h)/C + gL*E_L/C + Iapp(j,k+1)/C +Isyn(j,k+1)/C;
            V(j,k+1) = V(j,k)+0.5*h*(k1+k2);
            end
        end
    end
end
subplot(2,1,1);
plot(t,Isyn(1,:),t,Isyn(5,:));
ylabel('Synaptic Currents');
xlabel('Time');

subplot(2,1,2);

plot(t,V(1,:),t,V(2,:),t,V(3,:),t,V(4,:),t,V(5,:));
ylabel('Membrane Potential');
xlabel('Time');