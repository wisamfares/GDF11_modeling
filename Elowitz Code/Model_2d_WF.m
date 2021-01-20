function [Mtx,T]=Model_2d_WF(Xconc,Yconc,prm,Solver)

if ~exist('Solver','var')
    Solver=@ModelSolver_2step_WF;
end

A=prm(1:4);
B=prm(5:6);
KA=prm(7:14);
KB=prm(15:30);
E=prm(31:end);

Mtx=zeros(length(Xconc),length(Yconc));

for Xidx=1:length(Xconc)
    for Yidx=1:length(Yconc)
        [Mtx(Xidx,Yidx),~,T]=Solver([Xconc(Xidx);Yconc(Yidx)], A, B, KA, KB, E);
    end
end



end