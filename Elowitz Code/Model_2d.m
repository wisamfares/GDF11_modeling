function [Mtx,T]=Model_2d(Xconc,Yconc,prm,Solver)

if ~exist('Solver','var')
    Solver=@ModelSolver_2step;
end

A=prm(1:2);
B=prm(3:4);
KA=prm(5:8);
KB=prm(9:16);
E=prm(17:24);

Mtx=zeros(length(Xconc),length(Yconc));

for Xidx=1:length(Xconc)
    for Yidx=1:length(Yconc)
        [Mtx(Xidx,Yidx),~,T]=Solver([Xconc(Xidx);Yconc(Yidx)], ...
            A, B, KA, KB, E);
    end
end



end