function [Mtx,T]=Model_2d_CR(Xconc,Yconc,prm,Solver)

if ~exist('Solver','var')
    Solver=@ModelSolver_2step_CR;
end

A=prm(1:2);
B=prm(3:4);
C=prm(5:6);
kD_f=prm(7:10);
kD_r=prm(11:14);
kT_f=prm(15:22);
kT_r=prm(23:30);
kCL_f=prm(31:34);
kCL_r=prm(35:38);
kCD_f=prm(39:46);
kCD_r=prm(47:54);
E=prm(55:62);

Mtx=zeros(length(Xconc),length(Yconc)); %Store total signal per (X,Y)

for Xidx=1:length(Xconc)
    for Yidx=1:length(Yconc)
        [Mtx(Xidx,Yidx),~,T] =Solver([Xconc(Xidx);Yconc(Yidx)], A, B, C, ...Receptor Concentrations
            kD_f, kD_r, kT_f, kT_r, ... Rate constants
            kCL_f, kCL_r, kCD_f, kCD_r, E);
    end
end



end