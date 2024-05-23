function [time,speed,speed_cycle,time_cycle,time_cycle_size]=...
    speed_cycle_elab_v01(input_file_speedname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function read the recorded speed and extract an equivalent speed 
% cycle to
% be used for FCGR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input_file_speedname='speed_recorded.txt';
recorded_data=import_speedrecorded_v00(input_file_speedname);

time=recorded_data(:,1);
speed=recorded_data(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Speed cycle arrangement - only open cycle crack are good
% Speed first derivative computation
Dspeed_Dtime=zeros(length(time)-1,1);
sign_der=zeros(length(time)-1,1);
for i=1:length(time)-1
    Dspeed_Dtime(i,1)=(speed(i+1)-speed(i))/(time(i+1)-time(i));
end

temp=1;
Cycle_data(temp,1)=time(1);
Cycle_data(temp,2)=speed(1);

for i=2:length(time)-1
    sign_der(i,1)=Dspeed_Dtime(i-1,1)*Dspeed_Dtime(i,1); % searching for a derivative sign change
    if sign_der(i,1)<0
        temp=temp+1;
        Cycle_data(temp,1)=time(i);
        Cycle_data(temp,2)=speed(i);
    end
end

if speed(length(time))>speed(length(time)-1)
    temp=temp+1;
    Cycle_data(temp,1)=time(length(time));
    Cycle_data(temp,2)=speed(length(time));
end

% Effective detected cycle
[time_cycle_size,~]=size(Cycle_data);
time_cycle=Cycle_data(:,1);
speed_cycle=Cycle_data(:,2);

% figure(1)
% plot(time,speed,'k-','LineWidth',2);
% hold on
% plot(time_cycle,speed_cycle,'b-','LineWidth',2);
% grid on
% hold on 
% [size1_cycle,~]=size(cycle_stress);
% for i=1:size1_cycle
% plot([cycle_stress(i,1),cycle_stress(i,2)],[cycle_stress(i,3),cycle_stress(i,4)],'m-+','LineWidth',2);
% hold on
% end
% xlim([3900,4000]);
end
