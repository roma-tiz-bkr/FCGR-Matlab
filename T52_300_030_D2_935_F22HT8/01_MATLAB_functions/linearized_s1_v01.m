function [Pm_s1,Pb_s1,error_int]=linearized_s1_v01(path_FEA,S1_FEA_loc,atemp_macro)

% - S1 linearization on the flaw surface
tk_plate=path_FEA(end,1);
S1_pt_a0=interp1(path_FEA,S1_FEA_loc,atemp_macro);

pos_near_all=find(path_FEA<=atemp_macro);
pos_near=pos_near_all(end,1);  % starting point location

min_dfct=zeros(pos_near,1);

for cyc=1:pos_near
    path_S1_pos_near=path_FEA(cyc,1);
    S1_near=S1_FEA_loc(cyc,1);  

    S1_1st=S1_near+(S1_pt_a0-S1_near)/(atemp_macro-path_S1_pos_near)*(0-path_S1_pos_near);
    S1_2nd=S1_near+(S1_pt_a0-S1_near)/(atemp_macro-path_S1_pos_near)*(tk_plate-path_S1_pos_near);

    S1_lin=[];
    S1_lin=interp1([0;tk_plate],[S1_1st;S1_2nd],path_FEA);
    delta_fct=[];
    delta_fct=S1_lin-S1_FEA_loc;
    min_dfct(cyc,1)=min(delta_fct(1:pos_near,1));
end

pto_LIN=find(min_dfct>=-0.01,1);

if isempty(pto_LIN)==0    
    cyc=pto_LIN;
    path_S1_pos_near=path_FEA(cyc,1);
    S1_near=S1_FEA_loc(cyc,1);

    S1_1st=S1_near+(S1_pt_a0-S1_near)/(atemp_macro-path_S1_pos_near)*(0-path_S1_pos_near);
    S1_2nd=S1_near+(S1_pt_a0-S1_near)/(atemp_macro-path_S1_pos_near)*(tk_plate-path_S1_pos_near);

    Pm_s1=(S1_1st+S1_2nd)/2;
    Pb_s1=(S1_1st-S1_2nd)/2;
    error_int=0;
else
    error_int=1;
    S1_1st=1e9;
    S1_2nd=1e9;
    Pm_s1=1e9;
    Pb_s1=1e9;
end

end