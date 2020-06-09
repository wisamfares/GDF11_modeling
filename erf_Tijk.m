function [E_Tijk] = erf_Tijk(Ai_0, Lj, Bk_0,K_ij, gamma_j, K_ijk, T_ijk)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
E_Tijk = 0; 
for i = 1:2
    for j = 1:2 
        for k = 1:2
            X = (Ai_0(i) - sum(sum(T_ijk(i,:,:))))/(1+sum(gamma_j(j)*K_ij(i,:)*Lj(j)));
            Y = (Bk_0(k) - sum(sum(T_ijk(:,:,k))));
            Z = K_ijk(i,j,k)*gamma_j(j)*K_ij(i,j)*X*Lj(j)*Y-T_ijk(i,j,k);
            E_Tijk = E_Tijk + Z^2;
        end 
    end
end
            
end

