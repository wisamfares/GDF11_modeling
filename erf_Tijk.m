function [E_Tijk] = erf_Tijk(Ai_0, Lj, Bk_0,K_ij, gamma_j, K_ijk, T_ijk)
%Error Function as described in Elowtiz paper, Eq. 24
%   Note, Eq. 24 includes a minus sign instead of a plus sign, this is
%   most likely a typo, given the derivation in Eqs. 21-23. 
E_Tijk = 0; 


for i = 1:2
   
    W = 0;
    for j = 1:2
    kl = (1/gamma_j(j))*K_ij(i,j)*Lj(j);
    W = W + kl;
    end 
    
    for j = 1:2 
        for k = 1:2
            %X = (Ai_0(i) - sum(sum(T_ijk(i,:,:))))/(1+sum((1/gamma_j(j))*K_ij(i,:)*Lj(j))); %The K_ij*Lj sum is NOT summing correctly
            X = (Ai_0(i) - sum(sum(T_ijk(i,:,:))))/(1+W);
            Y = (Bk_0(k) - sum(sum(T_ijk(:,:,k))));
            Z = K_ijk(i,j,k)*(1/gamma_j(j))*K_ij(i,j)*X*Lj(j)*Y-T_ijk(i,j,k);
            E_Tijk = E_Tijk + Z^2;
        end 
    end
end
            
end

