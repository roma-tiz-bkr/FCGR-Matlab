function [Sy_mat,UTS_mat,A,m,DKth,KiC,error_flag_ALL,n_mat_par_string]=mat_prop_FCGR_v05(n_mat_par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% - Material properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v03 --> implemented switch rather than if

switch n_mat_par
    case 1 
        n_mat_par_string='F22HT8-ITN06710-mean-R0.5-H2-10bar';
        % Paris' data
%         A=0.0000000630678;    % [mm/cycle]
%         A=1.58786E-7;    % [mm/cycle] updated 06 Jul 23
% 	    A=8.31921E-8;    % [mm/cycle] updated 03 August 23 by Granta
%         A=8.31921E-8;    % [mm/cycle] updated 03 August 23 by Granta
        A=1.71744E-07;   % 30 bar!
        % [mm/cycle] updated 20 November 23 @50 bar!!!
%         m=3.389518198;
%         m=3.102739838; % updated 06 Jul 23
        m=3.1902720199103; % updated 03 August 23 by Granta
%         m=3.20; % updated 03 August 23 by Granta
%         DKth=	5.33;     % [MPa sqrt(m)]

        DKth=   5.47;     % [MPa sqrt(m)] updated 03 August 23 ESTIMATED by Granta
        KiC=65;       % [MPa sqrt(m)]
        % Tensile properties
        Sy_mat=734.79780047;   % [MPa]
        UTS_mat=944.84027256;  % [MPa]
    case 2 
        n_mat_par_string='OptiMo-ITN4200013+2sigma-R0.1-H2-70bar';
        % Paris' data
        A=8.859e-11;    % [mm/cycle]
        m=5.42;
        DKth=4.2;     % [MPa sqrt(m)]
        KiC=55;       % [MPa sqrt(m)]
        % Tensile properties
        Sy_mat=1210;   % [MPa]
        UTS_mat=1853;  % [MPa]
    case 3 
        n_mat_par_string='OptiMo-ITN4200013-air'; % estimated!!!!
        % Paris' data
        A=3.70051326305486e-08;    % [mm/cycle]
        m=2.58866501727698;
        DKth=7.42;     % [MPa sqrt(m)]
        KiC=136;       % [MPa sqrt(m)]
        % Tensile properties
        Sy_mat=1210;   % [MPa]
        UTS_mat=1853;  % [MPa]
    case 4 
        n_mat_par_string='API6ACRA-718-H2'; % estimated!!!!
        % Paris' data
        A=4e-17;    % [mm/cyc/(MPa sqrt(m))^(1/m)]
        m=8.3;
        DKth=4;     % [MPa sqrt(m)] NPD ASSUMPTION
        KiC=60;       % [MPa sqrt(m)]
        % Tensile properties
        Sy_mat=1034;   % [MPa]
        UTS_mat=1676;  % [MPa]
    otherwise
	    disp('Error! The requested material FCGR material properties are not found!');
        error_flag_ALL=1;
        
end
error_flag_ALL=0;
end