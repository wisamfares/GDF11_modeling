%This code requires the 'Optimization Toolbox' to use fmincon properly,
%otherwise you'll recieve an error


XC=logspace(-1.4,2.5,13);
YC=logspace(-0.9,3.5,13);

% XC=logspace(-4,4,13);
% YC=logspace(-4,4,13);

% Edited Parameters from Elowitz (4 A receptors)
% prms = [
%     0.07,0.0072,0.021,0.0022,...A_j receptors
%     14.5,0.002,...B_k receptors
%     0.55,0.45,0.032,0.97,0.7150,0.5850,0.0416,1.2610,...KA_ij
%     0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,...KB_ijk
%     0.012,0.48,0.097,0.047,0.18,0.044,0.087,0.048, 0.012,0.48,0.097,0.047,0.18,0.044,0.087,0.048...e_ijk
%     ];


% A = rand(1,2)
% A = [0.7*A 0.3*A]
% B = rand(1,2);
% KA = rand(1,4);
% KA = [KA 1.3*KA];
% KB = rand(1,8);
% KB = [KB KB];
% E = rand(1,8);
% E = [E E];
% 
% prms = [A B KA KB E]
 
prms = [
    0.07,0.0072,...A_j receptors
    14.5*10^-10,0.02*10^-10,...B_k receptors
    0.55,0.45,0.032,0.97,...KA_ij
    0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,...KB_ijk
    0.012,0.48,0.097,0.047,0.18,0.044,0.087,0.048...e_ijk
    ];
% 
% prms = [
%     0.0072,0.07,...A_j receptors
%     0.002,14.5,...B_k receptors
%     0,0, 0.55, 0.45,...KA_ij
%     0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,...KB_ijk
%     0.012,0.48,0.097,0.047,0.18,0.044,0.087,0.048...e_ijk
%     ];


prms = [
    3,0.001,...A_j receptors
    14.5,2,...B_k receptors
    0.55,0.45,0.032,0.97,...KA_ij
    0.10,0.039,0.02,0.14,0.28,0.0067,0.27,0.13,...KB_ijk
    0.012,0.48,0.097,0.047,0.18,0.044,0.087,0.048...e_ijk
    ];

figure(4)
MIMB=Model_2d_WF(XC, YC, prms);
imagesc(MIMB)
set(gca,'YDir','normal')

