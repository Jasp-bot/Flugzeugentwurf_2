clc
clear all
close all
%% Import durch Funktion

Berechnung_PS4_basis_stat_Massen;

%% Laden von Daten

load Projekt_specs.mat;
load Zwischenergebnisse_PS4_basis_stat_Massen.mat;
load Ergebnisse_Basis_stat_m.mat;




%% Berechnung der regressionsanalys

% f_gamma = fit((basis_stat_m.R_OPS_mittel), (basis_stat_m.gamma_DP), 'power1');
% coeffizienten_gamma = coeffvalues(basis_stat_m.f_gamma);
% 
%f_kappa = fit(basis_stat_m.R_OPS_mittel, basis_stat_m.kappa_DP, 'power1');
% coeffizienten_kappa = coeffvalues(basis_stat_m.f_kappa);


%% Plot
hold on
grid on

plot(basis_stat_m.f_gamma, basis_stat_m.R_OPS_mittel, basis_stat_m.gamma_DP, '.r');

plot(basis_stat_m.f_kappa, basis_stat_m.R_OPS_mittel, basis_stat_m.kappa_DP, '.b');



% plot Linie designrange


X(1,1) = specs.max_range_basis_km;
X(2,1) = specs.max_range_basis_km;
Y(1,1) = 0;
Y(2,1) = 0.7;
Z(1,1) = specs.max_range_Shrink_km;
Z(2,1) = specs.max_range_Shrink_km;
plot(X,Y);
plot(Z,Y);

plot(specs.max_range_Shrink_km, Ergebnis_basis_m.kappa_DP_Shrink,'*k');
plot(specs.max_range_Shrink_km, Ergebnis_basis_m.gamma_DP_Shrink,'*k');

xlabel('Reichweite in km','FontSize',16);
ylabel('Gamma oder Kappa in 1','FontSize',16);
title(['Kraftstoff- und Nutzlastmassenfaktor'],'FontSize',18);
legend('Gamma', 'Gammakurve', 'Kappa', 'Kapppakurve','Designreichweite', 'Designreichweite Shrink','Location','northwest','FontSize',18);

