%% Funktion zur berechnung von Gleichung 12 PS7 s.3 
% Funktion liefert das Schubzu-Standschubverhaeltnis
function [S_S0_KF_j] = S_S0_KF_j(D, rho_rho0_H, Ma_j, p_p0, bypass)
%     load Projekt_specs.mat;
%     Faktor = specs.Schubfaktor;         % zur erhöhung des Schubes
    S_S0_KF_j = D .* rho_rho0_H .* exp(-0.35 .* Ma_j .* p_p0 .* sqrt(bypass)); %  Faktor .*
end
