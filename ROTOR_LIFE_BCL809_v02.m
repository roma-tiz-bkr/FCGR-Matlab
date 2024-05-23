diary off
delete diary.log
diary diary.log
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - FCGR ROTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ROTOR_NAME='NEOM_BCL809';

ref_critical_phi=224;
ref_critical_imp_name=strcat('T52_',num2str(ref_critical_phi),'_030_D2_935_F22HT8');

ref_folder1='\T52_207_030_D2_935_F22HT8\';
ref_folder2='\T52_300_030_D2_935_F22HT8\';

input_file_speedname_ref='speed_recorded_const_809_days_speed.txt';
speed_file_loc=pwd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do you have a new rotor?
new_rotor='No'; % Answers: 'Yes' 'No'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(new_rotor,'Yes')==1
    %% - ROTOR GLOBAL SUMMARY REPORT
    txt1='*********************************';
    txt2='+++++++++ ROTOR REPORT ++++++++++';
    txt3='*********************************';
    txt4=strcat('Rotor name:',ROTOR_NAME);

    txt5=strcat('Rotor status: 1 - Green flag');

    txt6=strcat('Rotor working-days:',num2str(0));
    txt7=strcat('Rotor expected life [days]: NA');

    txt8=strcat('Rotor working-cycles: NA');
    txt9=strcat('Rotor expected cycles: NA');

    mylog_summary=fopen(strcat(speed_file_loc,'\ROTOR_SUMMARY_',ROTOR_NAME,'.txt'),'w');
    formatSpec3='%s\n';
    fprintf(mylog_summary,formatSpec3,txt1);
    fprintf(mylog_summary,formatSpec3,txt2);
    fprintf(mylog_summary,formatSpec3,txt3);
    fprintf(mylog_summary,formatSpec3,txt4);
    fprintf(mylog_summary,formatSpec3,txt5);
    fprintf(mylog_summary,formatSpec3,txt6);
    fprintf(mylog_summary,formatSpec3,txt7);
    fprintf(mylog_summary,formatSpec3,txt8);
    fprintf(mylog_summary,formatSpec3,txt9);

    fclose('all');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation of ROTOR condition: can we go?
filename_rotor_summary=strcat(speed_file_loc,'\ROTOR_SUMMARY_',ROTOR_NAME,'.txt');
ROTOR_STATUS = importfile_ROTOR_STATUS(filename_rotor_summary);
if ROTOR_STATUS(1,1)==1 % Can we switch on the rotor?
    input_file_speedname=strcat(speed_file_loc,'\',input_file_speedname_ref);

    cd(strcat(speed_file_loc,ref_folder1));
    NC_FCGR_launcher_T52_207_v03

    cd(strcat(speed_file_loc,ref_folder2));
    NC_FCGR_launcher_T52_300_v03

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ADD SECTION FOR GLOBAL RESULTS ELABORATION
    cd(strcat(speed_file_loc,ref_folder1));
    addpath(myfolder)
    for i=1:n_comp
        cd(strcat(speed_file_loc,ref_folder1,component_names{1,i},'\04_MATLAB_FCGR\'));
        filename=strcat('res_min_LIFE_T52_207_030_935_',component_names{1,i},'.txt');
        resminLIFE1{i,1} = importfile_resminlife_v00(filename);
        minlife1(i,1)=min(resminLIFE1{i,1}(:,2));
    end
    cd(strcat(speed_file_loc,ref_folder1));
    rmpath(myfolder)
    cd(strcat(speed_file_loc,ref_folder2));
    addpath(myfolder)
    for i=1:n_comp
        cd(strcat(speed_file_loc,ref_folder2,component_names{1,i},'\04_MATLAB_FCGR\'));
        filename=strcat('res_min_LIFE_T52_300_030_935_',component_names{1,i},'.txt');
        resminLIFE2{i,1} = importfile_resminlife_v00(filename);
        minlife2(i,1)=min(resminLIFE2{i,1}(:,2));
    end
    cd(strcat(speed_file_loc,ref_folder2));
    rmpath(myfolder)

    minlife_imp1=min(minlife1);
    minlife_imp2=min(minlife2);

    % Actually the rotor life is calculated with simulation of cycles. If the
    % rotor cycle is too short a fake life is appearing.
    LIFE_ROTOR=interp1([207 300],[minlife_imp1 minlife_imp2],ref_critical_phi);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% - ROTOR GLOBAL SUMMARY REPORT
    txt1='*********************************';
    txt2='+++++++++ ROTOR REPORT ++++++++++';
    txt3='*********************************';
    txt4=strcat('Rotor name:',ROTOR_NAME);

    if LIFE_ROTOR>ncycles_grow
        txt5=strcat('Rotor status: 1 - Green flag');
    else
        txt5=strcat('Rotor status: 0 - Red flag - maintenance required');
    end

    txt6=strcat('Rotor working-days:',num2str(time_ref(end,1)-time_ref(1,1)+ROTOR_STATUS(2,1)));
    txt7=strcat('Rotor expected life [days]:',num2str(round(LIFE_ROTOR/ncycles_grow*(time_ref(end,1)-time_ref(1,1)))));

    txt8=strcat('Rotor working-cycles:',num2str(ncycles_grow));
    txt9=strcat('Rotor expected cycles:',num2str(round(LIFE_ROTOR)));

    mylog_summary=fopen(strcat(speed_file_loc,'\ROTOR_SUMMARY_',ROTOR_NAME,'.txt'),'w');
    formatSpec3='%s\n';
    fprintf(mylog_summary,formatSpec3,txt1);
    fprintf(mylog_summary,formatSpec3,txt2);
    fprintf(mylog_summary,formatSpec3,txt3);
    fprintf(mylog_summary,formatSpec3,txt4);
    fprintf(mylog_summary,formatSpec3,txt5);
    fprintf(mylog_summary,formatSpec3,txt6);
    fprintf(mylog_summary,formatSpec3,txt7);
    fprintf(mylog_summary,formatSpec3,txt8);
    fprintf(mylog_summary,formatSpec3,txt9);

    fclose('all');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ROTOR_STATUS(1,1)==0% Unstable rotor identified
    disp('ROTOR NOT READY TO GO: unstable defects may occur: INSPECTION REQUIRED');
else 
    disp('Calculation cannot be performed!');
end

cd(speed_file_loc);
diary off