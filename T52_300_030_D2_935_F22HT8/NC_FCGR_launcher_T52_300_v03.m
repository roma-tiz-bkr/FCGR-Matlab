clearvars -except speed_file ROTOR_NAME ref_critical_phi ref_critical_imp_name ref_folder1 ref_folder2 speed_file_loc input_file_speedname input_file_speedname_ref new_rotor ROTOR_STATUS
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - FCGR tool - lr defined by disk2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% ---------------------------- GENERAL INPUTS -----------------------------
% -------------------------------------------------------------------------
% Impeller name
imp_name='T52_300_030_935_'; % impeller_phi_mach_diameter_

% Number of components to be calculated
n_comp=5; 

% Impeller component names
component_names=cell(1,n_comp);

component_names{1,1}='n_disk2'; % This component must be always the first for lr computation
component_names{1,2}='n_blade1';
component_names{1,3}='n_disk1';
component_names{1,4}='n_flowpath1';
component_names{1,5}='n_flowpath2';

% Folder for matlab function
myfolder='01_MATLAB_functions/';
myfolder_path='/03_ANSYS_path/';
myfolder_res='/04_MATLAB_FCGR/';

% Speed recorded input file name
% input_file_speedname='speed_recorded.txt';
% input_file_speedname='speed_recorded_short.txt';
% input_file_speedname='speed_recorded_const.txt';

% -------------------------------------------------------------------------
% ----------------------- INPUTS FOR FCGR ANALYSIS ------------------------
% -------------------------------------------------------------------------

% Output controls
plottami_or_not=0; % 1 for plots / 0 for no plots

% Initial semi-elliptical crack dimensions if rotor is new
a0=0.7;   %[mm]
c0=a0;    %[mm]

% Scaling geometry factor
SCALING_FCT=935/390;

% Chose material:
n_mat_par=1;    % type 1 for F22HT8-ITN06710+2sigma-R0.1-H2-70bar
                %      2 for OptiMo-ITN4200013+2sigma-R0.1-H2-70bar
                %      3 for OptiMo-ITN4200013-air

% Reference Load Step FEA speed [m/s]
FEA_ref_speed=[...
    0,...
40.0 ,...
80.0 ,...
120.0,...
160.0,...
200.0,...
240.0,...
265.0,...
305.0,...
345.0,...
385.0,...
441.7,...
450.5,...
    ]; %[m/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - EXECUTION PHASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

starting_dir=pwd;

fracture_res_ALL_comp=cell(n_comp,1);

for xxxxx=1:n_comp
    clearvars -except json_GT component_names xxxxx n_comp comp_name starting_dir imp_name myfolder myfolder_path myfolder_res starting_dir plottami_or_not  a0 c0 SCALING_FCT n_mat_par MCS_burst MCS_low LS_dir maxxxxx_ITER FEA_ref_speed speed_file ROTOR_NAME ref_critical_phi ref_critical_imp_name ref_folder1 ref_folder2 speed_file_loc input_file_speedname input_file_speedname_ref new_rotor ROTOR_STATUS fracture_res_ALL_comp
    comp_name=component_names{1,xxxxx}; % component
    lr_comp_ref=strcmp(component_names{1,xxxxx},'n_disk2');
    close all
    clc
    cd(starting_dir);
    json_GT=0;  % 0 = NOT print .json %   json_GT=1; % 1 = print .json for GT FCGR tool

    % Input macro matlab directory
    addpath(myfolder);

    % ---------------------------------------------------------------------
    % -------------------- RUN FOR FCGR ANALYSIS! -------------------------
    % ---------------------------------------------------------------------

    NC_BS7910_2019_ncloud_v21
    fracture_res_ALL_comp{xxxxx,1}=fracture_res_ALL;
    % Path restoring
    rmpath(myfolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Defect dimension archive at required life\minlife
cd(starting_dir);
addpath(myfolder);
for xxxxx=1:n_comp
    cd(strcat(starting_dir,'\',component_names{1,xxxxx},'\04_MATLAB_FCGR\'));
    filename=strcat('res_min_LIFE_T52_300_030_935_',component_names{1,xxxxx},'.txt');
    resminLIFE_loc{xxxxx,1} = importfile_resminlife_v00(filename);
    minlife_loc(xxxxx,1)=min(resminLIFE_loc{xxxxx,1}(:,2));
end

% Impeller min LIFE, compared to speed cycle life
iter_extr=min(min(minlife_loc),ncycles_grow);

for xxxxx=1:n_comp
    cd(strcat(starting_dir,'\',component_names{1,xxxxx},'\04_MATLAB_FCGR\'));
    n_tot_scn=[];
    def_dim_iter=[];
    [n_tot_scn,~]=size(fracture_res_ALL_comp{xxxxx,1});
    for rrr=1:n_tot_scn 
        if isempty(fracture_res_ALL_comp{xxxxx,1}{rrr,1})==0
            def_dim_iter(rrr,1)=fracture_res_ALL_comp{xxxxx,1}{rrr,1}(iter_extr,2); %a0 updated at required cycle
            def_dim_iter(rrr,2)=fracture_res_ALL_comp{xxxxx,1}{rrr,1}(iter_extr,6); %c0 updated at required cycle
        else
            % section stress data are not available
            def_dim_iter(rrr,1)=0; %a0 updated
            def_dim_iter(rrr,2)=0; %c0 updated
        end
    end

    %- Print defect sizes
    [def_dim_a,def_dim_b]=size(def_dim_iter);

    % save('def_dim_summary.mat','def_dim_a','def_dim_b');
    mylog_def=[];
    mylog_def=fopen(strcat(starting_dir,'\',component_names{1,xxxxx},'\04_MATLAB_FCGR\',imp_name,component_names{1,xxxxx},'_def_summary','.txt'),'w');
    formatSpec2='%15.8f %15.8f\r\n';
    fprintf(mylog_def,formatSpec2,def_dim_iter');
    fclose('all');
end