function [lr, kr]=failure_locus(Sy_mat,UTS_mat)

% DTS06.26_A_2 eq.4 and eq.5 as reference to compute the FAD

n_div_fail=101;
i=1;
lr(i,1)=0;
kr(i,1)=1;
for i=2:n_div_fail
    lr(i,1)=lr(i-1,1)+(Sy_mat+UTS_mat)/2/Sy_mat/100;
    kr(i,1)=(1-0.14*lr(i,1)^2)*(0.3+0.7*exp(-0.65*lr(i,1)^6));
end
i=n_div_fail+1;
    lr(i,1)=lr(i-1,1)+1e-9;
    kr(i,1)=0;
% figure
% plot(lr,kr)