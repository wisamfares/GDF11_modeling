%Ligand Space [Concentrations: (X,Y)]
XC=logspace(-2.5,2.5,13);
YC=logspace(-2.5,2.5,13);


%Original Parameters from Elowitz (C_h added)
prms = [
    0.07,0.0072,...A_j receptors
    14.5,0.002,...B_k receptors
    10^-5,10^-5,...C_h coreceptors
    1,1,1,1,...kD_f_ij
    0.55^-1,0.45^-1,0.032^-1,0.97^-1,...kD_r_ij
    1,1,1,1,1,1,1,1,...kT_f_ijk
    0.10^-1,0.039^-1,0.02^-1,0.14^-1,0.28^-1,0.0067^-1,0.27^-1,0.13^-1,...kT_r_ijk
    1,1,1,1,...kCL_f_hj
    5^-1,4^-1,8^-1,2^-1,...kCL_r_hij
    1,1,1,1,1,1,1,1,...kCD_f_hij
    0.08^-1,0.02^-1,0.012^-1,0.018^-1,0.03^-1,0.0014^-1,0.016^-1,0.118^-1,...kCD_r_hij
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
figure(10)
[MIMB,T] =Model_2d_CR(XC, YC, prms);
imagesc(MIMB)
set(gca,'YDir','normal')
set(gca,'dataAspectRatio',[1 1 1]);
expon = linspace(-2.5,2.5,13);
label = cell(1,length(XC));
for i = 1:length(XC)
    label{i} =['10^{' num2str(round(expon(i),2)) '}'];
end
set(gca,'XTick',1:2:13,'XTickLabel',label(1:2:13))
set(gca,'YTick',1:2:13, 'YTickLabel',label(1:2:13))

xlabel('Ligand 1')
ylabel('Ligand 2')


