%% Funktion Wellenwiderstand

function [delta_c_WM, delta_Ma] = Transsonischer_W(Ma_unendlich, c_A_F)

load Ergebnisse_Fluegel_Tank_NP.mat;

% PS4 S.5 Formel 18
k_vector = [0.758, 0.1, -0.090, 0, -0.100];

% for n_DD = 1:5
%     M_DD_profil_phi25_vec(1,n_DD) = k_vector(1,n_DD) .* c_A_F.^(n_DD-1);
% end    
M_DD_profil_phi25 = k_vector(1,1) .* c_A_F.^(0) + k_vector(1,2) .* c_A_F.^(1) +...
    k_vector(1,3) .* c_A_F.^(2) + k_vector(1,4) .* c_A_F.^(3) + k_vector(1,5) .* c_A_F.^(4);

% PS4 S.5 Formel 17
delta_Ma_mat = (Ma_unendlich - M_DD_profil_phi25./(sqrt(cos(Ergebnisse_Fluegel.phi_25_max))));
delta_Ma =(delta_Ma_mat);
% PS4 S.4 Formel 16
delta_c_WM = (0.002 * exp(60 * delta_Ma_mat));

end
