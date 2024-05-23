%% - CRACK PROPAGATION WITH NOT CONSTANT CYCLE
% BS 7910:2019 
% Secondary stress is not considered

%-----------------------
% - Actual macro version
%----------------------- Implemented functions version
% - Material properties evaluation: mat_prop_FCGR_v02.m
% - LIFING: LIFE_computation_v10.m
% - FAD material curve definition: failure_locus.m
% - Linearization: linearized_s1_v01.m
% - Postproc: Postproc_macro_v04.m
% - Speed cycle definition: speed_cycle_elab_v00.m
% - Import Ansys files: import_ALL_RES_SEC.m
% - Import Ansys files: importfile_Section_data.m
% - Ansys_macro1: ALL_n_path_data_extract_v11.dat
% - Ansys_macro2: wp_for_crack_grad_ver03.mac
% - INPUT: speed_recorded.txt
% -------------------------------------------------
% >>>>>>>>>>>>> INTRODUCED IN v17: <<<<<<<<<<<<<<<<
% -> RUN BY FCGR_LAUNCHER -------------------------
% -> input and output directories -----------------
% -> for loop on LS_start and LS_end --------------
% -> integrated .json file generation -------------
% -> added ALL_LS_min_life.txt file generation ----
% -> v16 works with Envelop_data_S1_delta_S1_v10 --
% -> v16 works with Postproc_macro_v05 ------------
% -> v16 works with mat_prop_FCGR_v03 ------------
% -------------------------------------------------

% Output controls
error_flag_ALL=0; % Check parametes for errors

%% Data import
% Input directories
myfolder1=strcat(comp_name,myfolder_path);
myfolder2=strcat(comp_name,myfolder_res);
addpath(myfolder1);

% - Import material properties
[Sy_mat,UTS_mat,A,m,DKth,KiC,error_flag_ALL,n_mat_par_string]=mat_prop_FCGR_v06(n_mat_par);

% - Import file names
% -> Section general information 
filename_Section_data='Section_data.txt';
Section_data = importfile_Section_data(filename_Section_data);
n_tot_scn=Section_data(end,1);      % Number of total evaluated section

% - Sections data extraction
n_start_postproc_LC=Section_data(1,end-2);
n_end_postproc_LC=Section_data(1,end-1);
% n_dir_postproc_LC=Section_data(1,end);

% -> Path data import
ALL_RES_SEC=cell(n_tot_scn,1);
for ww1=1:n_tot_scn
    ALL_RES_SEC{ww1,1} = import_ALL_RES_SEC(strcat('ALL_RES_SEC',num2str(ww1),'.txt'),n_end_postproc_LC-n_start_postproc_LC+1);
end

% -> Title picture name
title_pic_end=strcat(imp_name,component_names{1,xxxxx});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Speed cycle arrangement - only open cycle crack are good
local_loc=pwd;
cd(speed_file_loc);
[time_ref,speed_ref,speed_cycle,time_cycle,time_cycle_size]=...
    speed_cycle_elab_v01(input_file_speedname_ref);

% This cycle will be repeated to simulate the total CC train life
TOT_TIME_CYCLE=time_cycle(end)-time_cycle(1); %[days]
TARGET_TOTAL_LIFE=30*365; %[days]

cd(local_loc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - FCGR elaboration
NODE_start=zeros(n_tot_scn,1);
NODE_end=zeros(n_tot_scn,1);

fracture_res_ALL=cell(n_tot_scn,1);

% if lr_comp_ref~=1
%     load(strcat('n_disk2',myfolder_res,'Lr_new_max_disk2.mat'));
% end

% Import defect dimension
if strcmp(new_rotor,'Yes')==1
    def_dim_iter=[];
    def_dim_iter=ones(n_tot_scn,2);
    def_dim_iter(:,1)=def_dim_iter(:,1).*a0;
    def_dim_iter(:,2)=def_dim_iter(:,2).*c0;
    mylog_def=fopen(strcat(myfolder2,'\',imp_name,component_names{1,xxxxx},'_def_summary','.txt'),'w');
    formatSpec2='%15.8f %15.8f\n';
    fprintf(mylog_def,formatSpec2,def_dim_iter);
    fclose('all');
elseif strcmp(new_rotor,'No')==1
    def_dim_iter=[];
    filename_def_dim_iter=strcat(myfolder2,'\',imp_name,component_names{1,xxxxx},'_def_summary','.txt');
    def_dim_iter=importfile_def_dim_start(filename_def_dim_iter);
end
for rrr=1:n_tot_scn 
    clc
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    display(strcat('Process elaboration...',num2str(round(rrr/n_tot_scn*100)),'%'));
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

    % - Analyzed node
    NODE_start(rrr,1)=Section_data(rrr,2);
    NODE_end(rrr,1)=Section_data(rrr,3);

    if NODE_end(rrr,1)==0
        skippami=1;
    else
        skippami=0;
    end

    if skippami==0 

    % Plate geometry parameters
    B=Section_data(rrr,4)*SCALING_FCT;    %[mm]    Plate Thickness
    W=Section_data(rrr,5)*SCALING_FCT; %[mm]    Plate Width

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % - FEA Pm Pb computation - Linearization with first point and a0 point.

    path_FEA=[];
    path_FEA(:,1)=ALL_RES_SEC{rrr,1}(:,1).*SCALING_FCT; % The path is the same for all LS - its length is scaled with impeller sizes
    S1_FEA=zeros(length(path_FEA(:,1)),n_end_postproc_LC-n_start_postproc_LC+1);

    Pm_seqv_FEA=zeros(n_end_postproc_LC-n_start_postproc_LC+1,1);
    Pb_seqv_FEA=zeros(n_end_postproc_LC-n_start_postproc_LC+1,1);

    for ww2=1:n_end_postproc_LC-n_start_postproc_LC+1
        S1_FEA(:,ww2)=ALL_RES_SEC{rrr,1}(:,(ww2-1)*3+2);
        Pm_seqv_FEA(ww2,1)=ALL_RES_SEC{rrr,1}(1,(ww2-1)*3+3);
        Pb_seqv_FEA(ww2,1)=ALL_RES_SEC{rrr,1}(2,(ww2-1)*3+3);
    end
    %-----------$ Speed        $ Pm S1   $ Pb S1   $ Pm Seqv   $ Pb Seqv   $ S1-1st   $ S1-a0

    ALL_FEA_res=[FEA_ref_speed',Pm_seqv_FEA,Pb_seqv_FEA];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % - FEA interpolation for different speed -> linear interpolation
    n_pts_path=21;
    for ww4=1:n_pts_path
        S1_int(ww4,:)=interp1(FEA_ref_speed',S1_FEA(ww4,:)',speed_cycle);
    end
    Pm_seqv=interp1(ALL_FEA_res(:,1),ALL_FEA_res(:,2),speed_cycle,'linear');
    Pb_seqv=interp1(ALL_FEA_res(:,1),ALL_FEA_res(:,3),speed_cycle,'linear');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% - K calculation
    % cycle delta determination
    Dspeed_Dtime2=zeros(time_cycle_size-1,1);
    temp=0;
    for ww5=1:time_cycle_size-1
        Dspeed_Dtime2(ww5,1)=(speed_cycle(ww5+1)-speed_cycle(ww5))/(time_cycle(ww5+1)-time_cycle(ww5));
        if Dspeed_Dtime2(ww5,1)>0
            temp=temp+1;
            cycle_stress(temp,1)=time_cycle(ww5,1);    % min time cycle
            cycle_stress(temp,2)=time_cycle(ww5+1,1);  % max time cycle
            cycle_stress(temp,3)=speed_cycle(ww5,1);   % min speed cycle
            cycle_stress(temp,4)=speed_cycle(ww5+1,1); % max speed cycle
            cycle_stress(temp,5)=Pm_seqv(ww5+1,1);     % max Pm seqv cycle 
            cycle_stress(temp,6)=Pb_seqv(ww5+1,1);     % max Pb seqv cycle 
            S1_int_max(:,temp)=S1_int(:,ww5+1);        % max S1 path interpolated cycle
            S1_int_min(:,temp)=S1_int(:,ww5);          % min S1 path interpolated cycle

            DS1_FEA_int(:,temp)=S1_int_max(:,temp)-S1_int_min(:,temp);  % Delta S1 path interpolated cycle

            for ww7=1:length(S1_int_max(:,temp))
                if S1_int_max(ww7,temp)<0
                    S1_int_max_zero(ww7,temp)=0;
                else
                    S1_int_max_zero(ww7,temp)=S1_int_max(ww7,temp);
                end
            end

            for ww7=1:length(S1_int_min(:,temp))
                if S1_int_min(ww7,temp)<0
                    S1_int_min_zero(ww7,temp)=0;
                else
                    S1_int_min_zero(ww7,temp)=S1_int_min(ww7,temp);
                end
            end

    %         % Delta S1 path interpolated cycle
    %         for ww7=1:length(DS1_FEA_int(:,temp))
    %             if DS1_FEA_int(ww7,temp)<0
    %                 DS1_FEA_int_zero(ww7,temp)=0;
    %             else
    %                 DS1_FEA_int_zero(ww7,temp)=DS1_FEA_int(ww7,temp);
    %             end
    %         end        
            % The cycle is open close, then I assume to use cautiosly the abs
            % delta value
            DS1_FEA_int_zero(:,temp)=abs(S1_int_max_zero(:,temp)-S1_int_min_zero(:,temp));  % Delta S1 path interpolated cycle
        end
    end    

    [ncycles_grow,~]=size(cycle_stress); % Recorded cycles  it is equal for all nodes

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% - Failure reference locus
    [Lr_fail, kr_fail]=failure_locus(Sy_mat,UTS_mat); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% - Crack propagazion computation
    a_new(1,1)=def_dim_iter(rrr,1); % starting crack depth - a0 - size
    c_new(1,1)=def_dim_iter(rrr,2); % starting crack latera - c0 - size
    a_new(1,2)=cycle_stress(1,1); % starting time in which I have this crack [hours]

    flag=0;
    % Cycle repetition
    n_cyc_tot(rrr,1)=1;
    Delta_time=0;
    ncycles_start(rrr,1)=1;
	NC_LIFE_computation_v14

    % Cycles to be calculated
    maxxxxx_ITER=round(TARGET_TOTAL_LIFE/TOT_TIME_CYCLE); % number of cycle iterations to be performed

    while flag==0
        if  n_cyc_tot(rrr,1)>=maxxxxx_ITER
            flag=1;
            pause(0.25);
        else
        n_cyc_tot(rrr,1)=n_cyc_tot(rrr,1)+1;
        ncycles_start(rrr,1)=ncycles_grow*(n_cyc_tot(rrr,1)-1)+1;
        Delta_time=cycle_stress(end,1)*(n_cyc_tot(rrr,1)-1);
        NC_LIFE_computation_v14
        end
    end    

    fracture_res_ALL{rrr,1}=fracture_res;
    fracture_res=[];
    else
        fracture_res_ALL{rrr,1}=[];
        n_max_cycles(rrr)=1e9; % these nodes are not calculated!
    end 
end

%% - MIN LIFE CALCULATION - skip not performed nodes
[min_LIFE,index_NODE_MIN_LIFE]=min(n_max_cycles);
NODE_MIN_LIFE=NODE_start(index_NODE_MIN_LIFE,1);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(strcat('Impeller MIN LIFE ->',num2str(min_LIFE)));
disp(strcat('Node min life ->',num2str(NODE_MIN_LIFE)));
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% Worst Section?
WRST_SCN=find(Section_data(:,2)==NODE_MIN_LIFE,1); 

% Writing report
mylog=fopen(strcat(myfolder2,'\',imp_name,component_names{1,xxxxx},'_report','.txt'),'w');
formatSpec='%s %d\n';

fprintf(mylog,'Impeller min life -> %d\n',min_LIFE);
fprintf(mylog,'Node min life -> %d\n',NODE_MIN_LIFE);
fprintf(mylog,'Section min life -> %d\n',WRST_SCN);

if min_LIFE==ncycles_grow*n_cyc_tot(index_NODE_MIN_LIFE,1)
    msg=['Counter saturated: the rotor life is grather then the computed life.'];
    disp(msg);
    fprintf(mylog,'%s\n',msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Postproc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plottami_or_not==1
    NC_Postproc_macro_v05
else
    fileID = fopen(strcat(myfolder2,'res_min_LIFE_',title_pic_end,'.txt'),'w');
    fprintf(fileID, 'NODE min_life\n\n');
    for rrr=1:n_tot_scn
    tempami9999(rrr,1)=n_max_cycles(rrr);
    fprintf(fileID,'%f %f\n',[NODE_start(rrr,1), tempami9999(rrr,1)]);
    end
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Check global error procedure
if error_flag_ALL==1
    clear all
    close all
    clc
    msg='Errors are found!!!! Check the INPUT parameters definition.';
    disp(msg);
    fprintf(mylog,'%s\n',msg);
    delete(strcat('res_min_LIFE_',title_pic_end,'.txt'));
    delete(strcat('BS7910_2019_FAD_',title_pic_end,'.png'));
end
if error_int==1
    msg='Errors into the linearization are found!!!! Check the linearization step.';
    disp(msg);
    fprintf(mylog,'%s\n',msg);
end

if min(NODE_end)==0
    msg='The following sections are not extracted from Ansys:';
    disp(Section_data(NODE_end==0,1));
    disp(msg);
    fprintf(mylog,'%s %d\n',msg,Section_data(NODE_end==0,1));
else
    msg='All sections are well extracted from Ansys.';
    disp(msg);
    fprintf(mylog,'%s\n',msg);
end

close all

fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Defect dimension archive at required life\minlife

iter_extr=min(ncycles_grow,min_LIFE);
for rrr=1:n_tot_scn 
    if isempty(fracture_res_ALL{rrr,1})==0
        def_dim_iter(rrr,1)=fracture_res_ALL{rrr,1}(iter_extr,2); %a0 updated
        def_dim_iter(rrr,2)=fracture_res_ALL{rrr,1}(iter_extr,6); %c0 updated
    else
        % section stress data are not available
        def_dim_iter(rrr,1)=0; %a0 updated
        def_dim_iter(rrr,2)=0; %c0 updated
    end
end

%- Print defect sizes
[def_dim_a,def_dim_b]=size(def_dim_iter);

% save('def_dim_summary.mat','def_dim_a','def_dim_b');

mylog_def=fopen(strcat(myfolder2,'\',imp_name,component_names{1,xxxxx},'_def_summary','.txt'),'w');
formatSpec2='%15.8f %15.8f\n';
fprintf(mylog_def,formatSpec2,def_dim_iter);
fclose('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Printing and exporting results
% final_table=array2table([mycycles life_node_scn],'VariableNames',{'LS_end','LS_start','Min_Life','Worst_Node','Worst_SCN'});
% disp(final_table)
% 
% fileID = fopen(strcat(myfolder2,'ALL_LS_min_life_',comp_name,'.txt'),'w');
% formatSpec1='%9s %9s %9s %9s %9s\n';
% formatSpec2='%9d %9d %9d %9d %9d\n';
% txt1='LS end';
% txt2='LS start';
% txt3='min life';
% txt4='node';
% txt5='SCN';
% fprintf(fileID,formatSpec1,txt1,txt2,txt3,txt4,txt5);
% for i=1:size(mycycles,1)
% fprintf(fileID,formatSpec2,mycycles(i,1),mycycles(i,2),life_node_scn(i,1),life_node_scn(i,2),life_node_scn(i,3));
% end
% fclose(fileID);

%% .json file generation
if json_GT==1
    GT_tool_ALL_data_v04
end

% if lr_comp_ref==1
%     if min_LIFE==ncycles_grow*n_cyc_tot(index_NODE_MIN_LIFE,1)
%         [~,index_NODE_MIN_LIFE_disk2]=max(Lr_new_max_disk2(1,:));
%     else
%         index_NODE_MIN_LIFE_disk2=index_NODE_MIN_LIFE;
%     end
%     save(strcat(myfolder2,'Lr_new_max_disk2.mat'),'Lr_new_max_disk2','index_NODE_MIN_LIFE_disk2');
% end
%% Directories restoring
rmpath(myfolder1);