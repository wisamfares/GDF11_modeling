rng(0)

%epsilon relates to how "active" a 2*Receptor-Ligand complex is
%(phos:dephos), similar to K_D for enzymes.
epsilon_ijk = rand(2,2,2); %phosphorylation rate (element of symbol)/ desphosphorylation (gamma) of T_ijk 

epsilon_ijk(:,:,1) = [0.05 0.60; 0.02 0.15];
epsilon_ijk(:,:,2) = [0.03 0.01; 0.06 0.08];

% alpha = 1/sum(sum(sum(epsilon_ijk))); %alpha used for normalization so that sum of epsilons equals one
% epsilon_ijk = alpha*epsilon_ijk %normalization of epsilon with alpha


%K_ijk is the forward rate/reverse rate of (D_ij + B_k  <==> T_ijk)
K_ijk = rand(2,2,2);

% beta = sum(sum(sum(K_ijk)))
% K_ijk = 1/beta*K_ijk
K_ijk(:,:,1) = [0.13 0.19; 0.29 0.01];
K_ijk(:,:,2) = [0.02 0.06; 0.14 0.16];

Ai_0 = [.9 .44]; 
Bk_0 = [.29 .63];


%K_ij is the forward/reverse rate of (A_i + L_j <==> D_ij)
%K_ij = rand(2,2)
% K_ij = 1/gamma_j(1)*K_ij
% K_ij_2 = 1/gamma_j(2)*K_ij
K_ij = [0.16 0.02; 0.84 0.14];
gamma_j = [sum(K_ij(:,1)) sum(K_ij(:,2))];

Lj = [0 0]
S = zeros(7)
for l = -7:7
    for m = -7:7
        Lj(1) = 10^l;
        Lj(2) = 10^m;
        
        %Error Function
        T_ijk_0 = zeros(2,2,2);


        T_ijk = fmincon(@(T_ijk)erf_Tijk(Ai_0, Lj, Bk_0, K_ij, gamma_j, K_ijk, T_ijk), T_ijk_0);
        %E_Tijk = erf_Tijk(Ai_0, Lj, Bk_0, K_ij, K_ijk, T_ijk);

        S(l+8,m+8) = sum(sum(sum(epsilon_ijk.*T_ijk)));
    end 
end
S = flipud(S); %flip matrix about vertical axis (10^-7, 10^-7 in bottom left) 
heatmap(S/max(max(S)),'Colormap', parula, 'CellLabelColor','none')





