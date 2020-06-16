%This code requires the 'Optimization Toolbox' to use fmincon properly,
%otherwise you'll recieve an error

rng(2)

%epsilon relates to how "active" a 2*Receptor-Ligand complex is
%(phos:dephos), similar to K_D for enzymes.
epsilon_ijk = rand(2,2,2); %phosphorylation rate (element of symbol)/ desphosphorylation (gamma) of T_ijk 
alpha = 1/sum(sum(sum(epsilon_ijk))); %alpha used for normalization so that sum of epsilons equals one
epsilon_ijk = alpha*epsilon_ijk %normalization of epsilon with alpha


%K_ijk is the forward rate/reverse rate of (D_ij + B_k  <==> T_ijk)
K_ijk = rand(2,2,2);
beta = sum(sum(sum(K_ijk)))
K_ijk = 1/beta*K_ijk %normalization of K_ijk so that sum equals 1 


%Log uniform distribution of receptor concentrations, A (type II) and B
%(type I) [10^-3, 10^3]
a = 6*rand(1,2)-3;
b = 6*rand(1,2)-3;
Ai_0 = 10.^a;
Bk_0 = 10.^b;

%K_ij is the forward/reverse rate of (A_i + L_j <==> D_ij)
K_ij = rand(2,2)

gamma_j = [sum(K_ij(:,1)) sum(K_ij(:,2))];

K_ij = 1/gamma_j(1)*K_ij
K_ij_2 = 1/gamma_j(2)*K_ij



    % %Values from Elowitz Paper, Fig. S4
    % K_ijk(:,:,1) = [0.13 0.19; 0.29 0.01]; 
    % K_ijk(:,:,2) = [0.02 0.06; 0.14 0.16]; 
    % Ai_0 = [.90 .44];
    % Bk_0 = [.29 .63]; 
    % K_ij = [0.16 0.02; 0.84 0.14];
    % gamma_j = [1 1]


Lj = [0 0]
S = zeros(7)
for l = -7:7
    for m = -7:7
        Lj(1) = 10^l;
        Lj(2) = 10^m;
        
        %Error Function
        T_ijk_0 = zeros(2,2,2);


        T_ijk = fmincon(@(T_ijk)erf_Tijk(Ai_0, Lj, Bk_0, K_ij, gamma_j, K_ijk, T_ijk), T_ijk_0, [], [], [], [], T_ijk_0, Inf*ones(2,2,2));
        
        S(l+8,m+8) = sum(sum(sum(epsilon_ijk.*T_ijk)));
    end 
end
S = flipud(S); %flip matrix about vertical axis (10^-7, 10^-7 in bottom left) 

% %Lables for Heatmap
xvalues = {'10^-7' '10^-6' '10^-5' '10^-4' '10^-3'  '10^-2' '10^-1' '10^0' ...
     '10^1' '10^2' '10^3' '10^4' '10^5' '10^6' '10^7'}
yvalues = flip(xvalues)

%figure(1)
%subplot(5,5,w)
%h = heatmap(S/max(max(S)), 'Colormap', parula, 'CellLabelColor','none')

h = heatmap(xvalues, yvalues, S/max(max(S)), 'Colormap', parula, 'CellLabelColor','none');

% h.XLabel = 'Ligand 1 (log)'
% h.YLabel = 'Ligand 2 (log)'





