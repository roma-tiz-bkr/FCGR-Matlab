%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Postprocessing macro
% - INTRODUCED IN v05: output directory myfolder2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
subplot(1,2,1)
txt1='FAD failure locus';
plot(Lr_fail,kr_fail,'k-','DisplayName',txt1,'LineWidth',1);
hold on
for rrr=1:n_tot_scn
if isempty(fracture_res_ALL{rrr,1})==0
txt=strcat('NODE=',num2str(NODE_start(rrr,1)));
txt2=strcat('Crack starting point-',txt);
plot(fracture_res_ALL{rrr,1}(1,4),fracture_res_ALL{rrr,1}(1,1),'bO','DisplayName',txt2,'LineWidth',1);
hold on
if n_cyc_tot(rrr,1)>1
    txt3=strcat('Crack actual propagation-',txt);
    plot(fracture_res_ALL{rrr,1}(1:ncycles_grow,4),fracture_res_ALL{rrr,1}(1:ncycles_grow,1),'m-','DisplayName',txt3,'LineWidth',1);
    hold on
    txt4=strcat('Crack prevision propagation-',txt);
    plot(fracture_res_ALL{rrr,1}(1+ncycles_grow:end,4),fracture_res_ALL{rrr,1}(1+ncycles_grow:end,1),'b--','DisplayName',txt4,'LineWidth',1);
else
    txt3=strcat('Crack actual propagation-',txt);
    plot(fracture_res_ALL{rrr,1}(1:n_max_cycles(rrr),4),fracture_res_ALL{rrr,1}(1:n_max_cycles(rrr),1),'m-','DisplayName',txt3,'LineWidth',1);
end
end
end
title(strcat('Worst node -->',num2str(NODE_start(index_NODE_MIN_LIFE,1))));
grid on
xlabel('L_r');
ylabel('k_r=K_{I}^{max}/K_{IC}');

subplot(1,2,2)
txt1=strcat('Failure locus',n_mat_par_string);
plot(Lr_fail,kr_fail,'k-','DisplayName',txt1,'LineWidth',1);
hold on
rrr=index_NODE_MIN_LIFE;
txt=strcat('NODE=',num2str(NODE_start(rrr,1)));
txt2=strcat('Crack starting point-',txt);
plot(fracture_res_ALL{rrr,1}(1,4),fracture_res_ALL{rrr,1}(1,1),'bO','DisplayName',txt2,'LineWidth',1);
hold on
if n_cyc_tot(rrr,1)>1
    txt3=strcat('Crack actual propagation-',txt);
    plot(fracture_res_ALL{rrr,1}(1:ncycles_grow,4),fracture_res_ALL{rrr,1}(1:ncycles_grow,1),'m-','DisplayName',txt3,'LineWidth',1);
    hold on
    txt4=strcat('Crack prevision propagation-',txt);
    plot(fracture_res_ALL{rrr,1}(1+ncycles_grow:end,4),fracture_res_ALL{rrr,1}(1+ncycles_grow:end,1),'b--','DisplayName',txt4,'LineWidth',1);

    if min_LIFE==ncycles_grow*n_cyc_tot(index_NODE_MIN_LIFE,1)
        subplot(1,2,2)
        title(strcat('Actual life=',num2str(min_LIFE),'Cycles','--> Counter saturated.'));   
    else
        subplot(1,2,2)
        title(strcat('Actual life=',num2str(min_LIFE),'Cycles','- Estimated next maintenance in=',num2str(min_LIFE-ncycles_grow),'Cycles. Node min life=',num2str(NODE_MIN_LIFE)));
    end

else
    txt3=strcat('Crack actual propagation-',txt);    
    plot(fracture_res_ALL{rrr,1}(1:n_max_cycles(rrr),4),fracture_res_ALL{rrr,1}(1:n_max_cycles(rrr),1),'m-','DisplayName',txt3,'LineWidth',1);
    if min_LIFE==ncycles_grow*n_cyc_tot(index_NODE_MIN_LIFE,1)
        subplot(1,2,2)
        title(strcat('Actual life=',num2str(min_LIFE),'Cycles','--> Counter saturated.'));   
    else
        subplot(1,2,2)
        title(strcat('Actual life=',num2str(min_LIFE),'Cycles','--> Rotor inspection is required.'));     
    end
end

grid on
xlabel('L_r');
ylabel('k_r=K_{I}^{max}/K_{IC}');
legend show

% figure(1)
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf,'BS7910_2019_60speed_MCS_cycle','png');
figure(2)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,strcat(myfolder2,'BS7910_2019_FAD_',title_pic_end),'png');
saveas(gcf,strcat(myfolder2,'BS7910_2019_FAD_',title_pic_end),'fig');

fileID = fopen(strcat(myfolder2,'res_min_LIFE_',title_pic_end,'.txt'),'w');
fprintf(fileID, 'NODE min_life\n\n');
for rrr=1:n_tot_scn
tempami9999(rrr,1)=n_max_cycles(rrr);
fprintf(fileID,'%f %f\n',[NODE_start(rrr,1), tempami9999(rrr,1)]);
end
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Worst node path plot and data print - FEA data
% cycle_to_display=1;
% Worst cycle point
rrr=index_NODE_MIN_LIFE;

% % First cycle point linearization LOW SPEED
% path_FEA(:,1)=ALL_RES_SEC{rrr,1}(:,1).*SCALING_FCT; % The path is the same for all LS
% 
% ww2=LS_eval_start-n_start_postproc_LC+1;
% pos_temp=1;
% display_S1_FEA_zero(:,pos_temp)=S1_int_min_zero_store(:,cycle_to_display,rrr);
% display_S1_FEA(:,pos_temp)=S1_int_min_store(:,cycle_to_display,rrr);
% display_Pm_seqv_FEA(1,pos_temp)=ALL_RES_SEC{rrr,1}(1,(ww2-1)*3+3);
% display_Pb_seqv_FEA(1,pos_temp)=ALL_RES_SEC{rrr,1}(2,(ww2-1)*3+3);
% 
% % Second cycle point linearization HIGH SPEED
% ww2=LS_eval_end-n_start_postproc_LC+1;
% pos_temp=2;
% display_S1_FEA_zero(:,pos_temp)=S1_int_max_zero_store(:,cycle_to_display,rrr);
% display_S1_FEA(:,pos_temp)=S1_int_max_store(:,cycle_to_display,rrr);
% display_Pm_seqv_FEA(1,pos_temp)=ALL_RES_SEC{rrr,1}(1,(ww2-1)*3+3);
% display_Pb_seqv_FEA(1,pos_temp)=ALL_RES_SEC{rrr,1}(2,(ww2-1)*3+3);
% 
% % Delta cycle curve and points linearization
% TOT_DS1_zero=DS1_int_zero_store(:,cycle_to_display,rrr);
% TOT_DS1=DS1_int_store(:,cycle_to_display,rrr);

Pm_max_S1_trend=Pm_max_S1_store(:,rrr);
Pb_max_S1_trend=Pb_max_S1_store(:,rrr);
DPm_max_S1_trend=Delta_Pm_S1_max_store(:,rrr);
DPb_max_S1_trend=Delta_Pb_S1_max_store(:,rrr);
a_trend=a_store(:,rrr);

figure(3)
% subplot(1,4,1)
% txt=strcat('S_1 @',num2str(FEA_ref_speed(1,LS_eval_end-n_start_postproc_LC+1)),'[m/s] FEA');
% plot(path_FEA(:,1),display_S1_FEA_zero(:,2),'b-','DisplayName',txt,'LineWidth',1);
% hold on
% txt=strcat('S_1 @',num2str(FEA_ref_speed(1,LS_eval_start-n_start_postproc_LC+1)),'[m/s] FEA');
% plot(path_FEA(:,1),display_S1_FEA_zero(:,1),'m-','DisplayName',txt,'LineWidth',1);
% hold on
% txt=strcat('\Delta_{S_1}');
% plot(path_FEA(:,1),TOT_DS1_zero(:,1),'r-','DisplayName',txt,'LineWidth',1);
% grid on
% xlabel('Path length [mm]');
% ylabel('First principal stress [MPa]');
% legend show
% % title(strcat('Seqv: Pm=',num2str(round(display_Pm_seqv_FEA(1,2))),' [MPa] ',' Pb=',num2str(round(display_Pb_seqv_FEA(1,2))),' [MPa]'));

% subplot(1,4,2)
% txt=strcat('S_1 @',num2str(FEA_ref_speed(1,LS_eval_end-n_start_postproc_LC+1)),'[m/s] FEA not zero');
% plot(path_FEA(:,1),display_S1_FEA(:,2),'b-','DisplayName',txt,'LineWidth',1);
% hold on
% txt=strcat('S_1 @',num2str(FEA_ref_speed(1,LS_eval_start-n_start_postproc_LC+1)),'[m/s] FEA not zero');
% plot(path_FEA(:,1),display_S1_FEA(:,1),'m-','DisplayName',txt,'LineWidth',1);
% hold on
% txt=strcat('\Delta_{S_1} not zero');
% plot(path_FEA(:,1),TOT_DS1(:,1),'r-','DisplayName',txt,'LineWidth',1);
% grid on
% xlabel('Path length [mm]');
% ylabel('First principal stress [MPa]');
% legend show
% title(strcat('Seqv: Pm=',num2str(round(display_Pm_seqv_FEA(1,2))),' [MPa] ',' Pb=',num2str(round(display_Pb_seqv_FEA(1,2))),' [MPa]'));

subplot(1,3,1)
if min_LIFE>=ncycles_grow
    txt=strcat('Pm_{S_1}');
    plot(time_store(1:min_LIFE,rrr),Pm_max_S1_trend(1:min_LIFE,1),'b-','DisplayName',txt,'LineWidth',1);
    hold on
    txt=strcat('Pb_{S_1}');
    plot(time_store(1:min_LIFE,rrr),Pb_max_S1_trend(1:min_LIFE,1),'k-','DisplayName',txt,'LineWidth',1);
    hold on
    txt=strcat('\Delta Pm_{S_1}');
    plot(time_store(1:min_LIFE,rrr),DPm_max_S1_trend(1:min_LIFE,1),'b--','DisplayName',txt,'LineWidth',1);
    hold on
    txt=strcat('\Delta Pb_{S_1}');
    plot(time_store(1:min_LIFE,rrr),DPb_max_S1_trend(1:min_LIFE,1),'k--','DisplayName',txt,'LineWidth',1);
else
    txt=strcat('Pm_{S_1}');
    plot(time_store(1:min_LIFE+1,rrr),Pm_max_S1_trend(1:min_LIFE+1,1),'b-','DisplayName',txt,'LineWidth',1);
    hold on
    txt=strcat('Pb_{S_1}');
    plot(time_store(1:min_LIFE+1,rrr),Pb_max_S1_trend(1:min_LIFE+1,1),'k-','DisplayName',txt,'LineWidth',1);
    hold on
    txt=strcat('\Delta Pm_{S_1}');
    plot(time_store(1:min_LIFE+1,rrr),DPm_max_S1_trend(1:min_LIFE+1,1),'b--','DisplayName',txt,'LineWidth',1);
    hold on
    txt=strcat('\Delta Pb_{S_1}');
    plot(time_store(1:min_LIFE+1,rrr),DPb_max_S1_trend(1:min_LIFE+1,1),'k--','DisplayName',txt,'LineWidth',1);
end
grid on
xlabel('Time');
ylabel('Pm , Pb [MPa]');
legend show

subplot(1,3,2)
if min_LIFE>=ncycles_grow
txt='Crack length';
plot(time_store(1:min_LIFE,rrr) ,a_trend(1:min_LIFE,1),'b-','DisplayName',txt,'LineWidth',1);
else
txt='Crack length';
plot(time_store(1:min_LIFE+1,rrr) ,a_trend(1:min_LIFE+1,1),'b-','DisplayName',txt,'LineWidth',1);
end
grid on
xlabel('Time');
ylabel('Crack length [mm]');
legend show

subplot(1,3,3)
plot(time_ref,speed_ref,'k-','LineWidth',2);
hold on
plot(time_cycle,speed_cycle,'b-','LineWidth',2);
grid on
hold on 
[size1_cycle,~]=size(cycle_stress);
for i=1:size1_cycle
plot([cycle_stress(i,1),cycle_stress(i,2)],[cycle_stress(i,3),cycle_stress(i,4)],'m-+','LineWidth',2);
hold on
end
% xlim([0,50]);

figure(3)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,strcat(myfolder2,'Path_worst_data_',title_pic_end),'png');
saveas(gcf,strcat(myfolder2,'Path_worst_data_',title_pic_end),'fig');