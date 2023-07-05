%% PS7 Flugleistung 1
clc
clear all
close all


load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Basis_stat_m.mat;
load Ergebnisse_Widerstand_FE2.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Start_Landeanforderungen.mat;
addpath('Unterfunktionen Widerstand');

%% Annahmen

% Achtung noch Hardcodet!!!!!!!!!!!!!!!!!!!!!!!!!
 % stand 12uhr 03.07.23: annahme das es sich um eine Laufvariabele handelt

% Annahme das wir offdesign berechnen
c_A_j = Ergebnisse_Widerstand_FE2.c_A_ges_off_D;
c_A_F_j = Ergebnisse_Widerstand_FE2.c_A_F_off_D_vec;
c_W_inkomp_j = Ergebnisse_Widerstand_FE2.c_W_ges_off_D_inkomp;
j = Ergebnisse_Widerstand_FE2.Annahmen.stuetzstellen;

% Annahmen fuer Schub !!!!!!!!!!!!!!!!!!!!!!!! nachfragen nicht sicher
TO_Masse = 250000;
S0_GTO = startschub.S0 / (TO_Masse * specs.g);
Momentane_Masse = 200000;
GTO_G = TO_Masse / Momentane_Masse;



%% Horizontalflugdiagramm

% Flughoehe = specs.flight_level * 10^2 ;                         % in ft
% hoehe =round(unitsratio('m','ft')*Flughoehe);                     % in m

schritte = 20;

hoehe = (round(linspace(1000, specs.flight_level*10^2, schritte)));
hoehe_m = round(unitsratio('m','ft')*hoehe);


for n_datensatz = 1:schritte

[rho_rho0_H, T_H, a_H, p_p0_H] = Atmos_H(hoehe_m(1, n_datensatz));

% PS7 S2 Formel 3
v_EAS_j(n_datensatz, :) = sqrt(( (2) ./ (ISA.rho_0 .* c_A_j) ) .* Ergebnisse_stat_Flaechenbelastung.Fleachenbelastung);

% PS7 S2 Formel 4 Gleitzahl inkompressibel
eps_inkomp_j = c_W_inkomp_j ./ c_A_j;

% PS7 S2 Formel 4 Flugmachzahl Ma_j in TAS
v_TAS_j(n_datensatz, :) = (v_EAS_j(n_datensatz, :))./(sqrt(rho_rho0_H));
Ma_j(n_datensatz, :) = v_TAS_j(n_datensatz, :) ./ a_H;

% Transsonischer Widerstand mit Funktion aus PS4
for zv = 1:j
delta_c_WM(n_datensatz,zv) = Transsonischer_W(Ma_j(n_datensatz,zv), c_A_F_j(1,zv));
end

% Gesamtwiderstansbeiwert PS7 S3 Formel 10
c_W_j(n_datensatz, :) = c_W_inkomp_j + delta_c_WM(n_datensatz, :);


% PS7 S3 Formel 11 Gleitzahl kompressibel
eps_kompr_j(n_datensatz, :) = c_W_j(n_datensatz, :) ./ c_A_j; % eps_kompr_j = (S/G)_erf_j
S_G_erf_j = eps_kompr_j;%(n_datensatz, :);


%% Vorhandenes Schub-Gewichts-Verhaeltnis

% PS7 S3 Formel 12
S_S0_KF_j_1(n_datensatz, :) =  S_S0_KF_j(specs.Drosselgrad(1,1), rho_rho0_H, Ma_j(n_datensatz, :), p_p0_H, specs.bypass);


% Einlaufverluste PS7 S3 Formel 13
dp_p0 = 0.02;
S_S0_E(n_datensatz, :) = 1 - (1.3 + 0.25 * specs.bypass) * dp_p0; 

 
% PS7 S4 Formel 14
S_S0_j(n_datensatz, :) = S_S0_KF_j_1(n_datensatz, :) + S_S0_E(n_datensatz, :);

% PS7 S4 Formel 15

S_G_vorh_j(n_datensatz, :) = S_S0_j(n_datensatz, :) .* S0_GTO .* GTO_G;

%% Optimale Leistungszustände
% Spez Schubueberschuss SET

SET(n_datensatz, :) = S_G_vorh_j(n_datensatz, :) - eps_kompr_j(n_datensatz, :);







%% PLOTS
figure(1) % Horiontalflugdiagramm
hold on
grid on

ylabel('S/G');
xlabel('v_{EAS}');
xlim([0 300]);
ylim([0 1]);

plot(v_EAS_j(n_datensatz, :), S_G_vorh_j(n_datensatz, :), 'r');
plot(v_EAS_j(n_datensatz, :), S_G_erf_j(n_datensatz, :), 'b');
%plot(v_EAS_j(n_datensatz, :), eps_inkomp_j(n_datensatz, :), '--b');
hold off

figure(2) % SET
hold on 
grid on
ylabel('SET');
xlabel('v_{TAS}');
xlim([0 350]);
ylim([0 1]);

plot(v_TAS_j(n_datensatz, :), SET(n_datensatz, :), 'r');

hold off

end







%% Funktion zur Berechnung der Atmospherischen Daten
function [rho_rho0_H, T_H, a_H, p_p0_H] = Atmos_H(hoehe_m)
    load Ergebnisse_ISA_DATA.mat
    rho_rho0_H = ISA.rho(hoehe_m) / ISA.rho_0;
    T_H = ISA.T(hoehe_m);
    a_H = ISA.a(hoehe_m);
    p_p0_H = ISA.p(hoehe_m) / ISA.p0;
end



%% Funktion zur berechnung von Gleichung 12 PS7 s.3 
function [S_S0_KF_j] = S_S0_KF_j(D, rho_rho0_H, Ma_j, p_p0, bypass)
    
    S_S0_KF_j = D .* rho_rho0_H .* exp(-0.35 .* Ma_j .* p_p0 .* sqrt(bypass));
end


