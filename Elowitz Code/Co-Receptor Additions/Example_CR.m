%Ligand Space [Concentrations: (X,Y)]
XC=logspace(-1.4,2.5,13);
YC=logspace(-0.9,3.5,13);


%Original Parameters from Elowitz (C_h added)
% prms = [
%     0.07,0.0072,...A_j receptors
%     14.5,0.002,...B_k receptors
%     18.6,0.001,...C_h coreceptors
%     0.55,0.45,0.032,0.97,...kD_f_ij
%     0.12,0.22,0.012,0.12,...kD_r_ij
%     0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,...kT_f_ijk
%     0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,...kT_r_ijk
%     0.55,0.45,0.032,0.97,...kCL_f_hj
%     0.12,0.22,0.012,0.12,...kCL_r_hj
%     0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,...kCD_f_hij
%     0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,...kCD_r_hij
%     0.012,0.48,0.097,0.047,0.18,0.044,0.087,0.048...e_ijk
%     ];
% 
prms = [
    3,0.001,...A_j receptors
    14.5,2,...B_k receptors
    10^-12,10^-12,...C_h coreceptors %Using C_h = 0 will cause errors in fmincon
    0.55,0.45,0.032,0.97,...kD_f_ij
    1,1,1,1,...kD_r_ij %Changing the reverse rate to 1 allows for ratio for comparison to original
    0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,...kT_f_ijk
    1,1,1,1,1,1,1,1,...kT_r_ijk %Changing the reverse rate to 1 allows for ratio for comparison to original
    0.55,0.45,0.032,0.97,...kCL_f_hj
    1,1,1,1,...kCL_r_hj
    0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,...kCD_f_hij
    1,1,1,1,1,1,1,1,...kCD_r_hij
    0.012,0.48,0.097,0.047,0.18,0.044,0.087,0.048...e_ijk
    ];
    
% % 
% prms = lognrnd(1,2,[1,62]);
% % 
% prms(5:6) = [10^-10, 10^-10]

% prms(7:54) = prms(7:54)*3

prms(31:34) = prms(31:34)*3;
prms(35:38) = prms(35:38)*2;
prms(39:46) = prms(39:46)*3;
prms(47:54) = prms(47:54)*2;






%Plot Signaling Heatmap
figure(7)
title('Coreceptor')
[MIMB,T] =Model_2d_CR(XC, YC, prms);
imagesc(MIMB)
set(gca,'YDir','normal')

