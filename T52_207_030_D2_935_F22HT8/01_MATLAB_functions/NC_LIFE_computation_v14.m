%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - LIFE calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pigreco=pi;
for i=ncycles_start(rrr,1):ncycles_grow*n_cyc_tot(rrr,1)
    if flag==0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Linearization with defect size increment
        % if there are S1<0 they are considered as zero
        
        [Pm_max,Pb_max,error_int]=linearized_s1_v01(path_FEA(:,1),S1_int_max_zero(:,i-ncycles_grow*(n_cyc_tot(rrr,1)-1)),a_new(i,1));
        [Delta_Pm_S1_LIN,Delta_Pb_S1_LIN,error_int]=linearized_s1_v01(path_FEA(:,1),DS1_FEA_int_zero(:,i-ncycles_grow*(n_cyc_tot(rrr,1)-1)),a_new(i,1));
        
        time_store(i,rrr)=cycle_stress(i-ncycles_grow*(n_cyc_tot(rrr,1)-1),1)+cycle_stress(end,1)*(n_cyc_tot(rrr,1)-1); 
        Pm_max_S1_store(i,rrr)=Pm_max;
        Pb_max_S1_store(i,rrr)=Pb_max;
        Delta_Pm_S1_max_store(i,rrr)=Delta_Pm_S1_LIN;
        Delta_Pb_S1_max_store(i,rrr)=Delta_Pb_S1_LIN;
        a_store(i,rrr)=a_new(i,1);

        S1_int_max_store(:,i,rrr)=S1_int_max(:,i-ncycles_grow*(n_cyc_tot(rrr,1)-1));
        S1_int_min_store(:,i,rrr)=S1_int_min(:,i-ncycles_grow*(n_cyc_tot(rrr,1)-1));
        DS1_int_store(:,i,rrr)=DS1_FEA_int(:,i-ncycles_grow*(n_cyc_tot(rrr,1)-1));
        
        S1_int_max_zero_store(:,i,rrr)=S1_int_max_zero(:,i-ncycles_grow*(n_cyc_tot(rrr,1)-1));
        S1_int_min_zero_store(:,i,rrr)=S1_int_min_zero(:,i-ncycles_grow*(n_cyc_tot(rrr,1)-1));
        DS1_int_zero_store(:,i,rrr)=DS1_FEA_int_zero(:,i-ncycles_grow*(n_cyc_tot(rrr,1)-1));
        
        % Seqv stress cycle extraction
        Pm_max_seqv_cyc=cycle_stress(i-ncycles_grow*(n_cyc_tot(rrr,1)-1),5);        
        Pb_max_seqv_cyc=cycle_stress(i-ncycles_grow*(n_cyc_tot(rrr,1)-1),6);        
        
        if W>=2*(c_new(i,1)+B)
            alpha=(a_new(i,1)/B)/(1+(B/c_new(i,1)));
        else
            alpha=(2*a_new(i,1)/B)*(c_new(i,1)/W);            
        end
       
        seqv_max_ref(i,1)=(Pb_max_seqv_cyc+(Pb_max_seqv_cyc^2+9*(Pm_max_seqv_cyc^2)*(1-alpha)^2)^0.5)/(3*(1-alpha)^2);
        
%         if lr_comp_ref==1
            Lr_new_max(i,1)=seqv_max_ref(i,1)/Sy_mat;
%             Lr_new_max_disk2(i,rrr)=Lr_new_max(i,1);
%         else
%             Lr_new_max(i,1)=Lr_new_max_disk2(i,index_NODE_MIN_LIFE_disk2); % The disk2 lr is used for FCGR of other components
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ARRIVATO QUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - coefficients calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEMBRANE loading calculation
        % Conditions check
        if a_new(i,1)/(2*c_new(i,1)) > 1
            disp('You are violating membrane and bending crack conditions');
        end
        
        if a_new(i,1)/(2*c_new(i,1)) <= 0.1
            if (a_new(i,1)/B) >= 1.25*((a_new(i,1)/c_new(i,1))+0.6)
                disp('You are violating membrane and bending crack conditions');
            end
        elseif a_new(i,1)/(2*c_new(i,1)) > 0.1 || a_new(i,1)/(2*c_new(i,1)) <= 1
            if (a_new(i,1)/B) >= 1
                disp('You are violating membrane and bending crack conditions');
            end
        end
        fw=(sec((pigreco*c_new(i,1)/W)*(a_new(i,1)/B)^0.5))^0.5;
        if a_new(i,1)/(2*c_new(i,1)) <=0.5 || a_new(i,1)/(2*c_new(i,1)) >0       
            M_1=1.13-0.09*a_new(i,1)/c_new(i,1);
            M_2=(0.89/(0.2+a_new(i,1)/c_new(i,1)))-0.54;
            M_3=0.5-1/(0.65+a_new(i,1)/c_new(i,1))+14*(1-a_new(i,1)/c_new(i,1))^24;
        elseif a_new(i,1)/(2*c_new(i,1)) >0.5   
            M_1=(c_new(i,1)/a_new(i,1))^0.5*(1+0.04*(c_new(i,1)/a_new(i,1)));
            M_2=0.2*(c_new(i,1)/a_new(i,1))^4;
            M_3=-0.11*(c_new(i,1)/a_new(i,1))^4;
        end
            
        % MEMBRANE semplification theta=90 --> crack front (deepest point)
        g_90=1;
        if a_new(i,1)/(2*c_new(i,1)) <=0.5 || a_new(i,1)/(2*c_new(i,1)) >0
            f_theta_90=1;
        elseif a_new(i,1)/(2*c_new(i,1)) >0.5 % || a_new(i,1)/(2*c_new(i,1)) <=1 
            f_theta_90=(c_new(i,1)/a_new(i,1))^0.5;
        end
        
        % MEMBRANE semplification theta=0 --> crack end
        if a_new(i,1)/(2*c_new(i,1)) <=0.5 || a_new(i,1)/(2*c_new(i,1)) >0
            g_0=1.1+0.35*(a_new(i,1)/B)^2;
            f_theta_0=(a_new(i,1)/c_new(i,1))^0.5;
        elseif a_new(i,1)/(2*c_new(i,1)) >0.5 % || a_new(i,1)/(2*c_new(i,1)) <=1 
            g_0=1.1+0.35*(c_new(i,1)/a_new(i,1))*(a_new(i,1)/B)^2;
            f_theta_0=1;
        end
        
        % Complete elliptic integral of the second kind
        if a_new(i,1)/(2*c_new(i,1)) <=0.5 || a_new(i,1)/(2*c_new(i,1)) >0
            PHI=(1+1.464*(a_new(i,1)/c_new(i,1))^1.65)^0.5;
        elseif a_new(i,1)/(2*c_new(i,1)) >0.5 % || a_new(i,1)/(2*c_new(i,1)) <=1 
            PHI=(1+1.464*(c_new(i,1)/a_new(i,1))^1.65)^0.5;
        end
        Mm_0=(M_1+M_2*(a_new(i,1)/B)^2+M_3*(a_new(i,1)/B)^4)*g_0*f_theta_0*fw/PHI;
        Mm_90=(M_1+M_2*(a_new(i,1)/B)^2+M_3*(a_new(i,1)/B)^4)*g_90*f_theta_90*fw/PHI;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BENDING loading calculation
        if a_new(i,1)/(2*c_new(i,1)) <=0.5 || a_new(i,1)/(2*c_new(i,1)) >0
            q=0.2+(a_new(i,1)/c_new(i,1))+0.6*(a_new(i,1)/B);
            H1=1-0.34*(a_new(i,1)/B)-0.11*(a_new(i,1)/c_new(i,1))*(a_new(i,1)/B);
            G1=-1.22-0.12*(a_new(i,1)/c_new(i,1));
            G2=0.55-1.05*(a_new(i,1)/c_new(i,1))^0.75+0.47*(a_new(i,1)/c_new(i,1))^1.5;
        elseif a_new(i,1)/(2*c_new(i,1)) >0.5 % || a_new(i,1)/(2*c_new(i,1)) <=1 
            q=0.2+(c_new(i,1)/a_new(i,1))+0.6*(a_new(i,1)/B);
            H1=1-(0.04+0.41*(c_new(i,1)/a_new(i,1)))*(a_new(i,1)/B)+...
                (0.55-1.93*(c_new(i,1)/a_new(i,1))^0.75+1.38*(c_new(i,1)/a_new(i,1))^1.5)*(a_new(i,1)/B)^2;
            G1=-2.11+0.77*(c_new(i,1)/a_new(i,1));
            G2=0.55-0.72*(c_new(i,1)/a_new(i,1))^0.75+0.14*(c_new(i,1)/a_new(i,1))^1.5;            
        end
        H2=1+G1*(a_new(i,1)/B)+G2*(a_new(i,1)/B)^2;
        H_0=H1+(H2-H1)*sin(deg2rad(0))^q;
        H_90=H1+(H2-H1)*sin(deg2rad(90))^q;
        
        Mb_0=H_0*Mm_0;
        Mb_90=H_90*Mm_90;

        ktm=1; % not set because FEA results (to be confirmed)
        ktb=1; % not set because FEA results (to be confirmed)
        
        Mkm=1; % welding not present
        Mkb=1; % welding not present
        
        km=1;  % not set because FEA results (to be confirmed)
        
        M=1;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - kr calculation        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 0 deg
%         YSp_min_0=M*(ktm*Mkm*Mm_0*Pm_min+ktb*Mkb*Mb_0*(Pb_min+(km-1)*Pm_min)); % Primary stress
%         YSp_min_90=M*(ktm*Mkm*Mm_90*Pm_min+ktb*Mkb*Mb_90*(Pb_min+(km-1)*Pm_min)); % Primary stress
%         YSs_min_0=0;                                                          % Secondary stress ---> YSs=Mm*Qm+Mb*Qb; NOT CONSIDERED
%         YSs_min_90=0;                                                          % Secondary stress ---> YSs=Mm*Qm+Mb*Qb; NOT CONSIDERED
        
        YSp_max_0=M*(ktm*Mkm*Mm_0*Pm_max+ktb*Mkb*Mb_0*(Pb_max+(km-1)*Pm_max)); % Primary stress
        YSp_max_90=M*(ktm*Mkm*Mm_90*Pm_max+ktb*Mkb*Mb_90*(Pb_max+(km-1)*Pm_max)); % Primary stress
        YSs_max_0=0;                                                          % Secondary stress ---> YSs=Mm*Qm+Mb*Qb;
        YSs_max_90=0;                                                          % Secondary stress ---> YSs=Mm*Qm+Mb*Qb;
        
%         YS_min_0=YSp_min_0+YSs_min_0;
        YS_max_0=YSp_max_0+YSs_max_0;        
%         YS_min_90=YSp_min_90+YSs_min_90;
        YS_max_90=YSp_max_90+YSs_max_90;
        
%         Ki_cycle_min_0=YS_min_0*sqrt(pigreco*a_new(i,1)/1000);
%         Ki_cycle_min_90=YS_min_90*sqrt(pigreco*a_new(i,1)/1000);
        Ki_cycle_max_0=YS_max_0*sqrt(pigreco*c_new(i,1)/1000);
        Ki_cycle_max_90=YS_max_90*sqrt(pigreco*a_new(i,1)/1000);
        
%         kr_min_0(i,1)=Ki_cycle_min_0/KiC; % kr min in that cycle
%         kr_min_90(i,1)=Ki_cycle_min_90/KiC; % kr min in that cycle
        kr_max_0(i,1)=Ki_cycle_max_0/KiC; % kr max in that cycle
        kr_max_90(i,1)=Ki_cycle_max_90/KiC; % kr max in that cycle

        YDeltaS_p_0=M*(ktm*Mkm*Mm_0*(Delta_Pm_S1_LIN)+ktb*Mkb*Mb_0*((Delta_Pb_S1_LIN)+(km-1)*(Delta_Pm_S1_LIN)));
        YDeltaS_p_90=M*(ktm*Mkm*Mm_90*(Delta_Pm_S1_LIN)+ktb*Mkb*Mb_90*((Delta_Pb_S1_LIN)+(km-1)*(Delta_Pm_S1_LIN)));
        YDeltaS_s_0=0; % Secondary stress ---> YDeltaS_s=Mm*(Qm_max-Qm_min)+Mb*(Qb_max-Qb_min);        
        YDeltaS_s_90=0; % Secondary stress ---> YDeltaS_s=Mm*(Qm_max-Qm_min)+Mb*(Qb_max-Qb_min);        
        
        YDeltaS_0=YDeltaS_p_0+YDeltaS_s_0;
        YDeltaS_90=YDeltaS_p_90+YDeltaS_s_90;
        Delta_Ki_0=YDeltaS_0*sqrt(pigreco*c_new(i,1)/1000);
        Delta_Ki_90=YDeltaS_90*sqrt(pigreco*a_new(i,1)/1000);
        Delta_kri_0=Delta_Ki_0/KiC;
        Delta_kri_90=Delta_Ki_90/KiC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% - Crack propagation computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Lr_new_max(i,1)>=max(Lr_fail)
            KIC_fail_LOCAL=1e-9; % artificial KIC because Lr max is reached
        else
            KIC_fail_LOCAL=interp1(Lr_fail,kr_fail,Lr_new_max(i,1))*KiC;
        end
        Ki_cycle_max=max(Ki_cycle_max_0,Ki_cycle_max_90);
        kr_cycle_max=max(kr_max_0(i,1),kr_max_90(i,1));
        Delta_Ki_max=max(Delta_Ki_0,Delta_Ki_90);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fracture_res definition summary
        fracture_res(i,1)=kr_cycle_max; %kr max
        fracture_res(i,2)=a_new(i,1); % a0 crack size [mm]
        fracture_res(i,3)=a_new(i,2); % time at which I have this new crack size
        fracture_res(i,4)=Lr_new_max(i,1); % lr
        fracture_res(i,5)=Delta_kri_90; % Delta kri 90  
        fracture_res(i,6)=c_new(i,1); % c0 crack size [mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if Delta_Ki_max>DKth && Ki_cycle_max < KIC_fail_LOCAL
            a_new(i+1,1)=a_new(i,1)+A*Delta_Ki_90^m; % [mm]
            c_new(i+1,1)=c_new(i,1)+A*Delta_Ki_0^m;  % [mm]
            a_new(i+1,2)=cycle_stress(i-ncycles_grow*(n_cyc_tot(rrr,1)-1),2)+Delta_time; %  [days] time at which I have this new crack size
            
%             check_me(i,1)=a_new(i,1);
%             check_me(i,2)=2*c_new(i,1);
%             check_me(i,3)=Pm_max;
%             check_me(i,4)=Delta_Pm_S1_LIN;
%             check_me(i,5)=Pb_max;
%             check_me(i,6)=Delta_Pb_S1_LIN;
%             check_me(i,7)=0;
%             check_me(i,8)=0;
%             check_me(i,9)=0;
%             check_me(i,10)=0;
%             check_me(i,11)=PHI;
%             check_me(i,12)=M_1;
%             check_me(i,13)=M_2;
%             check_me(i,14)=M_3;
%             check_me(i,15)=g_0;
%             check_me(i,16)=g_90;
%             check_me(i,17)=f_theta_0;
%             check_me(i,18)=f_theta_90;
%             check_me(i,19)=fw;
%             check_me(i,20)=Mm_0;
%             check_me(i,21)=Mm_90;
%             check_me(i,22)=H1;
%             check_me(i,23)=H2;
%             check_me(i,24)=H2;
%             check_me(i,25)=H2;
%             check_me(i,26)=Mb_0;
%             check_me(i,27)=Mb_90;
%             check_me(i,28)=kr_max_0(i,1);
%             check_me(i,29)=Delta_kri_0;
%             check_me(i,30)=kr_max_90(i,1);
%             check_me(i,31)=Delta_kri_90;
%             check_me(i,32)=0;
%             check_me(i,33)=0;
%             check_me(i,34)=0;
%             check_me(i,35)=0;
%             check_me(i,36)=0;
%             check_me(i,37)=0;
%             check_me(i,38)=0;
%             check_me(i,39)=0;
%             check_me(i,40)=0;
%             check_me(i,41)=0;
%             check_me(i,42)=kr_max_0(i,1);
%             check_me(i,43)=Delta_kri_0;
%             check_me(i,44)=kr_max_90(i,1);
%             check_me(i,45)=Delta_kri_90;
%             check_me(i,46)=Pm_max_seqv_cyc;
%             check_me(i,47)=Pb_max_seqv_cyc;
%             check_me(i,48)=alpha;
%             check_me(i,49)=alpha;
%             check_me(i,50)=alpha;
%             check_me(i,51)=seqv_max_ref(i,1);
%             check_me(i,52)=kr_cycle_max;
%             check_me(i,53)=Lr_new_max(i,1);
%             check_me(i,54)=KIC_fail_LOCAL/KiC;
%             check_me(i,55)=Delta_Ki_0;
%             check_me(i,56)=0;
%             check_me(i,57)=Delta_Ki_90;
%             check_me(i,58)=0;
%             check_me(i,59)=A;
%             check_me(i,60)=m;
%             check_me(i,61)=a_new(i+1,1);
%             check_me(i,62)=c_new(i+1,1)*2;
%             check_me(i,63)=KIC_fail_LOCAL/KiC-kr_cycle_max;
              n_max_cycles(rrr)=i;
        elseif Ki_cycle_max >= KIC_fail_LOCAL
            flag=1;
            disp('Crack can be unstable!! Please turn off the CC train and check your rotor.');
            n_max_cycles(rrr)=i-1;
            i=ncycles_grow; % to end the for loop
        else
            % crack does not growth
            a_new(i+1,1)=a_new(i,1);
            c_new(i+1,1)=c_new(i,1);
            a_new(i+1,2)=cycle_stress(i-ncycles_grow*(n_cyc_tot(rrr,1)-1),2)+Delta_time; %  [hours] time at which I have this new crack size
            n_max_cycles(rrr)=i;
        end
    end
end