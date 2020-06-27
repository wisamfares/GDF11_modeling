function [Sig,err,T]=ModelSolver_2step(L, A, B, KA, KB, E)

%
% vector order: A1L1B1, A2L1B1, A1L2B1, A2L2B1, and same with B2.
%

err=0;
% initialization
if ~exist('L', 'var')
    L=[1,0];
end
if ~exist('A', 'var')
    A = [1,1];
end
if ~exist('B', 'var')
    B = [1,1];
end

NA = length(A);
NL = length(L);
NB = length(B);
N = NA*NL*NB;

if ~exist('KA', 'var')
    KA = 1 * ones(1,NA*NL);
end
if ~exist('KB', 'var')
    KB = 1 * ones(1,NA*NL*NB);
end
if ~exist('E', 'var')
    E = 1 * ones(1,NA*NL*NB);
end

Ai = reshape(A,NA,1,1);
Lj = reshape(L,1,NL,1);
Bk = reshape(B,1,1,NB);
KAij = reshape(KA,NA,NL,1);
KBijk = reshape(KB,NA,NL,NB);
E = reshape(E,1,NA*NL*NB);


%constraints:
Acon=[1,1,1,1,0,0,0,0;
    0,0,0,0,1,1,1,1];
Bcon=B;
%upper bound:
UB=reshape(repmat(B,4,1),[],1);
LB=zeros(N,1);

% initial guess:
%
T0 = zeros(N,1);
options= optimoptions('lsqnonlin','display','off','MaxFunEvals',8000,'MaxIter',4000);
optionsFMC= optimoptions('fmincon','display','off','MaxFunEvals',8000,'MaxIter',4000);
% try first with no bounds as it is faster
[T,~,~,exitflag] = lsqnonlin(@errfunc, T0, [], [], options);
count=0;
while any(T<0) || any(sum(sum(reshape(T,2,2,2),1),2)>Bk)
    [T,~,exitflag] = fmincon(@(x) sum(errfunc(x).^2), T0/(2^count), Acon, Bcon, [],[],LB,UB,[],optionsFMC);
    T0=T;
    [T,~,~,exitflag] = lsqnonlin(@errfunc, T, [], [], options);    
    if count>5
        break
    end
    count=count+1;
end
if any(T<0) || any(sum(sum(reshape(T,2,2,2),1),2)>Bk)
    warning('Negative receptor Level. L: %d %d, A: %d %d B: %d %d',L, A, B)
    err=1;
end
if exitflag==0
    warning('No minima reached. L: %d %d, A: %d %d B: %d %d',L, A, B)
    err=1;
end
    
Sig=E*T;

    function err = errfunc(x)
        xijk=reshape(x,NA,NL,NB);
        
        %calculate the various elements, and duplicate to a full 3d
        %matrices (i,j,k) if needed
        Cijk = KBijk;
        Ck = Bk-sum(sum(xijk,1),2);
        Ck = repmat(Ck, [NA,NL,1]);
        Cij = KAij;
        Cij = repmat(Cij,[1, 1, NB]);
        Ci = (Ai-squeeze(sum(sum(xijk,2),3))) ./ (1+KAij*Lj');
        Ci = repmat(Ci,[1, NL, NB]);
        Cj = Lj;
        Cj = repmat(Cj,[NA, 1, NB]);
        
        err = reshape(Cijk.*Ck.*Cij.*Ci.*Cj,[],1) - x;
    end
end
