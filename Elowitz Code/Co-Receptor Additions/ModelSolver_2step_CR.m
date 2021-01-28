function [Sig,err,T]=ModelSolver_2step_CR(L, A, B, C, kD_f, kD_r, kT_f, kT_r, ...
                                          kCL_f, kCL_r, kCD_f, kCD_r, E)

%
% vector order: A1L1B1, A2L1B1, A1L2B1, A2L2B1, and same with B2.
%

err=0;
%% Intialize Variables
if ~exist('L', 'var')
    L=[1,0];
end
if ~exist('A', 'var')
    A = [1,1];
end
if ~exist('B', 'var')
    B = [1,1];
end
if ~exist('C', 'var')
    C = [1,1];
end

NA = length(A);
NL = length(L);
NB = length(B);
NC = length(C); %check if NC is required for N calculation
N = NA*NL*NB + NL*NC; %N will be the length of input vector for error function (12 entries for (2,2,2,2))

if ~exist('kD_f', 'var')
    KA = 1 * ones(1,NA*NL);
end
if ~exist('kD_r', 'var')
    KA = 1 * ones(1,NA*NL);
end

if ~exist('kT_f', 'var')
    KB = 1 * ones(1,NA*NL*NB);
end
if ~exist('kT_r', 'var')
    KB = 1 * ones(1,NA*NL*NB);
end

if ~exist('kCL_f', 'var')
    KB = 1 * ones(1,NL*NC);
end
if ~exist('kCL_r', 'var')
    KB = 1 * ones(1,NL*NC);
end

if ~exist('kCD_f', 'var')
    KB = 1 * ones(1,NA*NL*NC);
end
if ~exist('kCD_r', 'var')
    KB = 1 * ones(1,NA*NL*NC);
end

if ~exist('E', 'var')
    E = 1 * ones(1,NA*NL*NB);
end
%% Reshape (CHECK FOR DIMENSION ALIGNMENT ESPECIALLY WITH MATRIX MULTIPLICATION)
%worst case, add 4th dimension to reshapes to reserve dimension for Ch
Ai = reshape(A,1,NA,1,1);
Lj = reshape(L,1,1,NL,1);
Bk = reshape(B,1,1,1,NB);

Ch = reshape(C,NC,1,1,1); %See if this configuration makes sense, add 4th dim if necessary
kD_fij = reshape(kD_f,1,NA,NL,1);
kD_rij = reshape(kD_r,1,NA,NL,1);
kT_fijk = reshape(kT_f,1,NA,NL,NB); %kT_f and kT_r can be RECOMBINED into KT (no direct relation to CL)
kT_rijk = reshape(kT_r,1,NA,NL,NB);
kCL_fhj = reshape(kCL_f,NC,1,NL,1); %Again, check config of reshape
kCL_rhj = reshape(kCL_r,NC,1,NL,1);
kCD_fhij = reshape(kCD_f,NC,NA,NL,1); %This one might get tricky, this will need to be used with B,
kCD_rhij = reshape(kCD_r,NC,NA,NL,1); %NOTE: consider making 1st dimension C_h, to keep order
E = reshape(E,1,NA*NL*NB);

%% Constraints:

%Acon multiples by x (T_ijk) in fmincon to return amount of B_k there is in
%each trimer, note Acon * T_ijk = 2x8 * 8x1 = 2x1 for direct comparison to
%B vector (renamed to Bcon)
Acon=[  1,1,0,0,0,0,0,0,0,0,0,0; %Concentration check for C_1
        0,0,1,1,0,0,0,0,0,0,0,0; %Concentration check for C_2
        0,0,0,0,1,1,1,1,0,0,0,0; %Concentration check for B_1
        0,0,0,0,0,0,0,0,1,1,1,1]; %Concentration cehck for B_2

%Keeping consistent with variable names for fmincon
Bcon=[C'; B'];

%upper & lower bounds:
UB = ones(12,1);
UB(1:4) = reshape(repmat(Ch,2,1),[],1);
UB(5:12)=reshape(repmat(Bk,4,1),[],1); %T_ijk < B_k for all k
LB=zeros(N,1); %Lower bound for T_ijk is obviously 0 

% initial guess:
%
T0 = zeros(N,1); %length NA*NL*NB + NC*NL (12 for (2,2,2,2) case
options= optimoptions('lsqnonlin','display','off','MaxFunEvals',8000,'MaxIter',4000);
optionsFMC= optimoptions('fmincon','display','off','MaxFunEvals',8000,'MaxIter',4000);
% try first with no bounds as it is faster -- NOT SURE IF THIS IS GOOD
[T,~,~,exitflag] = lsqnonlin(@errfunc_CR, T0, [], [], options);
count=0;

T_check = T(5:12); %Need to add a similiar check for CL_hj < Ch
CL_check = T(1:4);
while any(T<0) || any(sum(sum(reshape(T_check,1,2,2,2),2),3)>Bk) || any(sum(reshape(CL_check,2,1,2,1),3)>Ch)
    [T,~,exitflag] = fmincon(@(x) sum(errfunc_CR(x).^2), T0/(2^count), Acon, Bcon, [],[],LB,UB,[],optionsFMC);
    T0=T;
    [T,~,~,exitflag] = lsqnonlin(@errfunc_CR, T, LB, [], options);    
    if count>5
        break
    end
    count=count+1;
    
    T_check = T(5:12); %Need to add a similiar check for CL_hj < Ch
    CL_check = T(1:4);
end

T_check = T(5:12);%Need to add a similiar check for CL_hj < Ch
CL_check = T(1:4); 
if any(T<0) || any(sum(sum(reshape(T_check,1,2,2,2),2),3)>Bk) || any(sum(reshape(CL_check,2,1,2,1),3)>Ch)
    warning('Negative receptor Level. L: %d %d, A: %d %d B: %d %d',Lj(1), Lj(2), Ai(1), Ai(2) , Bk(1), Bk(2))
    err=1;
end
if exitflag==0
    warning('No minima reached. L: %d %d, A: %d %d B: %d %d',Lj, Ai, Bk)
    err=1;
end
    
Sig=E*T(5:12);

    function err = errfunc_CR(x) %x has 12 entries, 4 CL_hj, 8 T_ijk
        y = x(1:4); %temporary CL_hj storage from input (first 4 entries)
        z = x(5:12);%temporary T_ijk storage from input (last 8 entries)
        y_hj = reshape(y,NC,1,NL,1);
        z_ijk=reshape(z,1,NA,NL,NB);
        KT = kT_fijk./kT_rijk;
        
        Vk = Bk -sum(sum(z_ijk,2),3); %Bk accounting
       
        
        Vh = Ch - sum(y_hj,3); %Ch accounting
        
        Den_ij = kD_rij + sum(kCD_rhij.*Vh, 1); %Denominator in most equation components 
        Num_ij = kD_fij.*Lj + sum(kCD_fhij.*y_hj,1); %Numerator in most equation components
        
        Vij = Num_ij ./ Den_ij; %Numerator/Denominator (2x2) 
        
     
        Vi = (Ai - squeeze(sum(sum(z_ijk,3),4))); %Ai accounting
        
        
        Vhj = Vh.*(kCL_fhj.*Lj + sum(kCD_rhij.*Vi.*Vij./(1+sum(Vij,3)),2))./(kCL_rhj + sum(kCD_fhij.*Vi./(1+sum(Vij,3)),2));   %Calculated CL_hj at given parameter stage
        
        Vijk = KT.*(Vi./(1+sum(Vij,3))).*Vij.*Vk; %Calculated T_ijk at given parameter stage
        
        Vhj = reshape(Vhj, NC, 1 , NL, 1);
        Vijk = reshape(Vijk, 1, NA, NL, NB);
        
        err = ones(1,12);
        
        err(1:4) = Vhj - y_hj;
        err(5:12) = Vijk - z_ijk;
    end
end
%%