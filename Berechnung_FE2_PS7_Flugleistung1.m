function Berechnung_FE2_PS7_Flugleistung1
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
load Ergebnisse_Massen_FE2.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Hochauftrieb_1.mat;
load Ergebnisse_Hochauftrieb_2.mat;
addpath('Unterfunktionen Widerstand');

%% Annahmen

% Bitte noch verändern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_A_max = Ergebnisse_stat_Flaechenbelastung.C_A_CR; %%%%%%%%%%%%%%% Achtung random wert, bitte von mac geben lassen
c_A_max_LDG = HA2.CA_MAX_LDG;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Achtung noch Hardcodet!!!!!!!!!!!!!!!!!!!!!!!!!
 % stand 12uhr 03.07.23: annahme das es sich um eine Laufvariabele handelt

% Annahme das wir offdesign berechnen
c_A_j = Ergebnisse_Widerstand_FE2.c_A_ges_off_D;
%c_A_j = Ergebnisse_Widerstand_FE2.c_A_F_off_D;
c_A_F_j = Ergebnisse_Widerstand_FE2.c_A_F_off_D;
c_W_inkomp_j = Ergebnisse_Widerstand_FE2.c_W_ges_off_D_inkomp;
j = Ergebnisse_Widerstand_FE2.Annahmen.stuetzstellen;

% Annahmen fuer Schub !!!!!!!!!!!!!!!!!!!!!!!! nachfragen nicht sicher
TO_Masse =  Ergebnisse_Massen_FE2.M_TO;
G_TO = TO_Masse * specs.g;
k_CR = 0.98;
Momentane_Masse_ICA = TO_Masse * 0.98 * specs.g;
Momentane_Masse_DEC = (Ergebnisse_Massen_FE2.M_TO - Ergebnisse_Massen_FE2.M_Z_Tripfuel)*specs.g;
Momentane_Masse_CL = TO_Masse * specs.g;
Momentane_Masse_CL_ALT = (Ergebnisse_Massen_FE2.M_TO - Ergebnisse_Massen_FE2.M_Z_Tripfuel)*specs.g;
Masse_ICA = TO_Masse * 0.98;
Masse_CL = TO_Masse;
Masse_DEC = (Ergebnisse_Massen_FE2.M_TO - Ergebnisse_Massen_FE2.M_Z_Tripfuel);
Masse_CL_ALT = (Ergebnisse_Massen_FE2.M_TO - Ergebnisse_Massen_FE2.M_Z_Tripfuel)
GTO_G_ICA = G_TO / Momentane_Masse_ICA;
GTO_G_CL = 1;
GTO_G_DEC = G_TO / Momentane_Masse_DEC;
GTO_G_CL_ALT = G_TO / Momentane_Masse_CL_ALT;
S0 = k_CR * G_TO * (1/schub_CR.Eta) / (schub_CR.S_S0_CR * schub_CR.S_S0_E);
S0_GTO = S0/G_TO;

G = Momentane_Masse_ICA; %%%%%%%%%%%%%%%%%%% ACHTUNG

% Annahmen fuer Flugbereichsdiagramm



hoehe_CR = round(unitsratio('m','ft')*(specs.flight_level*10^2));
hoehe_ALT = round(unitsratio('m','ft')*(specs.flight_level_ALT*10^2));

%% Horizontalflugdiagramm

% Flughoehe = specs.flight_level * 10^2 ;                         % in ft
% hoehe =round(unitsratio('m','ft')*Flughoehe);                     % in m

schritte = 100;
hoehe_plus = 18;
hoehe = (round(linspace(1500, ((specs.flight_level + hoehe_plus)*10^2), schritte))); % 1000: 3000: specs.flight_level*10^2;
% hoehe = 1000: 1000: specs.flight_level*10^2;
hoehe_m = round(unitsratio('m','ft').*hoehe);

% Plot Luca
% Step_Plot= abs(hoehe_m(1,1)- hoehe_m(1,2) );

% if (Step_Plot ~= 0)
%     figure()
%     grid on 
%     hold on
%     xlim([0 350])
%     ylim([0 0.4])
%     Legend=cell(2*Step_Plot+2,1);%  two positions 
% end
% k=0;
% ende code Luca

%% Erstellungen der Vektoren

[rho_rho0_H, T_H, a_H, p_p0_H, rho_H] = Atmos_H(hoehe_m.');


% PS7 S2 Formel 3
v_EAS_j = sqrt(( (2)./ (ISA.rho_0 .* c_A_j) ) .* (Momentane_Masse_ICA ./Ergebnisse_Fluegel.F));
v_EAS_j_CL = sqrt(( (2) ./ (ISA.rho_0 .* c_A_j) ) .* (Momentane_Masse_CL ./Ergebnisse_Fluegel.F));
v_EAS_j_DEC = sqrt(( (2) ./ (ISA.rho_0 .* c_A_j) ) .* (Momentane_Masse_DEC./Ergebnisse_Fluegel.F));

v_EAS_j_CL_ALT = sqrt(( (2) ./ (ISA.rho_0 .* c_A_j) ) .* (Momentane_Masse_CL_ALT ./Ergebnisse_Fluegel.F));



% PS7 S2 Formel 4 Flugmachzahl Ma_j in TAS
v_TAS_j = (v_EAS_j)./(sqrt(rho_rho0_H));
v_TAS_j_CL = (v_EAS_j_CL)./(sqrt(rho_rho0_H));
v_TAS_j_DEC = (v_EAS_j_DEC)./(sqrt(rho_rho0_H));

v_TAS_j_CL_ALT = (v_EAS_j_CL_ALT)./(sqrt(rho_rho0_H));


Ma_j = v_TAS_j ./ a_H;
Ma_j_CL = v_TAS_j_CL ./ a_H;
Ma_j_DEC = v_TAS_j_DEC ./ a_H;

Ma_j_CL_ALT = v_TAS_j_CL_ALT ./ a_H;


% Vorhandenes Schub-Gewichts-Verhaeltnis

% PS7 S3 Formel 12
S_S0_KF_j_CR = S_S0_KF_j(specs.Drosselgrad(1,3), rho_rho0_H, Ma_j, p_p0_H, specs.bypass);
S_S0_KF_j_CL = S_S0_KF_j(specs.Drosselgrad(1,2), rho_rho0_H, Ma_j_CL, p_p0_H, specs.bypass);
S_S0_KF_j_DEC = S_S0_KF_j(specs.Drosselgrad(1,3), rho_rho0_H, Ma_j_DEC, p_p0_H, specs.bypass);

S_S0_KF_j_CL_ALT = S_S0_KF_j(specs.Drosselgrad(1,2), rho_rho0_H, Ma_j_CL_ALT, p_p0_H, specs.bypass);


% Einlaufverluste PS7 S3 Formel 13
dp_p0 = 0.02;
S_S0_E = 1 - (1.3 + 0.25 * specs.bypass) * dp_p0; 

 
% PS7 S4 Formel 14
S_S0_j = S_S0_KF_j_CR .* S_S0_E;
S_S0_j_CL = S_S0_KF_j_CL .* S_S0_E;
S_S0_j_DEC = S_S0_KF_j_DEC .* S_S0_E;

S_S0_j_CL_ALT = S_S0_KF_j_CL_ALT .* S_S0_E;

% PS7 S4 Formel 15

S_G_vorh_j = S_S0_j .* S0_GTO .* GTO_G_ICA;
S_G_vorh_j_CL = S_S0_j_CL .* S0_GTO .* GTO_G_CL;
S_G_vorh_j_DEC = S_S0_j_DEC .* S0_GTO .* GTO_G_DEC;

S_G_vorh_j_CL_ALT = S_S0_j_CL_ALT .* S0_GTO .* GTO_G_CL_ALT;

% SFC
[sfc_daNh_CR,  sfc_1PERh_CR, sfc_1PERs_Horizontalflug] = SFC_vec(hoehe_m, Ma_j, specs.bypass);
[sfc_daNh_CL,  sfc_1PERh_CL, sfc_1PERs_Horizontalflug_CL] = SFC_vec(hoehe_m, Ma_j_CL, specs.bypass);
[sfc_daNh_DEC, sfc_1PERh_DEC, sfc_1PERs_Horizontalflug_DEC] = SFC_vec(hoehe_m, Ma_j_DEC, specs.bypass);

[sfc_daNh_CL_ALT,  sfc_1PERh_CL_ALT, sfc_1PERs_Horizontalflug_CL_ALT] = SFC_vec(hoehe_m, Ma_j_CL_ALT, specs.bypass);



% % Anzahl der Plots
% numPlots = length(hoehe_m);
% 
% % Farbverlauf definieren
% colorStart = [1, 0.2, 0];   % Startfarbe (RGB)
% colorEnd = [0, 0, 1];     % Endfarbe (RGB)
% 
% % Farbwerte für jeden Plot berechnen
% colors_1 = zeros(numPlots, 3);
% colors_2 = zeros(numPlots, 3);
% for i = 1:numPlots
%     colors_1(i, :) = colorStart + (i-1) * (colorEnd - colorStart) / (numPlots-1);
%     colors_2(i, :) = colorStart + (i-1) * (colorEnd - colorStart) / (numPlots-1);
% end
% 
% 



for n_datensatz = 1:length(hoehe_m)


% Transsonischer Widerstand mit Funktion aus PS4
% for zv = 1:j
% delta_c_WM(n_datensatz,zv)= Transsonischer_W(Ma_j(n_datensatz,zv), c_A_F_j(1,zv));
% delta_c_WM_CL(n_datensatz,zv) = Transsonischer_W(Ma_j_CL(n_datensatz,zv), c_A_F_j(1,zv));
% delta_c_WM_DEC(n_datensatz,zv) = Transsonischer_W(Ma_j_DEC(n_datensatz,zv), c_A_F_j(1,zv));
% end


delta_c_WM(n_datensatz, :)= Transsonischer_W(Ma_j(n_datensatz,:), c_A_F_j);
delta_c_WM_CL(n_datensatz, :) = Transsonischer_W(Ma_j_CL(n_datensatz,:), c_A_F_j);
delta_c_WM_DEC(n_datensatz, :) = Transsonischer_W(Ma_j_DEC(n_datensatz,:), c_A_F_j);

delta_c_WM_CL_ALT(n_datensatz, :) = Transsonischer_W(Ma_j_CL_ALT(n_datensatz,:), c_A_F_j);


% Gesamtwiderstansbeiwert PS7 S3 Formel 10
c_W_j(n_datensatz, :) = c_W_inkomp_j + delta_c_WM(n_datensatz, :);
c_W_j_CL(n_datensatz, :) = c_W_inkomp_j + delta_c_WM_CL(n_datensatz, :);
c_W_j_DEC(n_datensatz, :) = c_W_inkomp_j + delta_c_WM_DEC(n_datensatz, :);

c_W_j_CL_ALT(n_datensatz, :) = c_W_inkomp_j + delta_c_WM_CL_ALT(n_datensatz, :);


% PS7 S3 Formel 11 Gleitzahl kompressibel
eps_kompr_j(n_datensatz, :) = c_W_j(n_datensatz, :) ./ c_A_j; % eps_kompr_j = (S/G)_erf_j
eps_kompr_j_CL(n_datensatz, :) = c_W_j_CL(n_datensatz, :) ./ c_A_j;
eps_kompr_j_DEC(n_datensatz, :) = c_W_j_DEC(n_datensatz, :) ./ c_A_j;

eps_kompr_j_CL_ALT(n_datensatz, :) = c_W_j_CL_ALT(n_datensatz, :) ./ c_A_j;

S_G_erf_j(n_datensatz, :) = c_W_j(n_datensatz, :) ./ c_A_j;
S_G_erf_j_CL(n_datensatz, :) = c_W_j_CL(n_datensatz, :) ./ c_A_j;
S_G_erf_j_DEC(n_datensatz, :) = c_W_j_DEC(n_datensatz, :) ./ c_A_j;

S_G_erf_j_CL_ALT(n_datensatz, :) = c_W_j_CL_ALT(n_datensatz, :) ./ c_A_j;



%% Optimale Leistungszustände
%% Spez Schubueberschuss SET

SET(n_datensatz, :) = S_G_vorh_j(n_datensatz, :) - eps_kompr_j(n_datensatz, :);
[Hochpunkte.SET_y(n_datensatz,:), Hochpunkte.SET_x(n_datensatz,:)] = max(SET(n_datensatz,:));
TAS_SET_H_vec(n_datensatz,:) = [v_TAS_j(n_datensatz, Hochpunkte.SET_x(n_datensatz,1)), Hochpunkte.SET_y(n_datensatz,1), hoehe_m(1,n_datensatz)];


SET_CL(n_datensatz, :) = S_G_vorh_j_CL(n_datensatz, :) - eps_kompr_j_CL(n_datensatz, :);
[Hochpunkte.SET_CL_y(n_datensatz,:), Hochpunkte.SET_CL_x(n_datensatz,:)] = max(SET_CL(n_datensatz,:));
TAS_SET_H_CL_vec(n_datensatz,:) = [v_TAS_j_CL(n_datensatz, Hochpunkte.SET_CL_x(n_datensatz,1)), Hochpunkte.SET_CL_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SET_DEC(n_datensatz, :) = S_G_vorh_j_DEC(n_datensatz, :) - eps_kompr_j_DEC(n_datensatz, :);
[Hochpunkte.SET_DEC_y(n_datensatz,:), Hochpunkte.SET_DEC_x(n_datensatz,:)] = max(SET_DEC(n_datensatz,:));
TAS_SET_H_DEC_vec(n_datensatz,:) = [v_TAS_j_DEC(n_datensatz, Hochpunkte.SET_DEC_x(n_datensatz,1)), Hochpunkte.SET_DEC_y(n_datensatz,1), hoehe_m(1,n_datensatz)];




%% Spezifischer Leistungueberschuss SEP

SEP(n_datensatz, :) = v_TAS_j(n_datensatz, :) .* (S_G_vorh_j(n_datensatz, :) - eps_kompr_j(n_datensatz, :));
[Hochpunkte.SEP_y(n_datensatz,:), Hochpunkte.SEP_x(n_datensatz,:)] = max(SEP(n_datensatz,:));
TAS_SEP_H_vec(n_datensatz,:) = [v_TAS_j(n_datensatz, Hochpunkte.SEP_x(n_datensatz,1)), Hochpunkte.SEP_y(n_datensatz,1), hoehe_m(1,n_datensatz)];


SEP_CL(n_datensatz, :) = v_TAS_j_CL(n_datensatz, :) .* (S_G_vorh_j_CL(n_datensatz, :) - eps_kompr_j_CL(n_datensatz, :));
[Hochpunkte.SEP_CL_y(n_datensatz,:), Hochpunkte.SEP_CL_x(n_datensatz,:)] = max(SEP_CL(n_datensatz,:));
TAS_SEP_H_CL_vec(n_datensatz,:) = [v_TAS_j_CL(n_datensatz, Hochpunkte.SEP_CL_x(n_datensatz,1)), Hochpunkte.SEP_CL_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SEP_DEC(n_datensatz, :) = v_TAS_j_DEC(n_datensatz, :) .* (S_G_vorh_j_DEC(n_datensatz, :) - eps_kompr_j_DEC(n_datensatz, :));
[Hochpunkte.SEP_DEC_y(n_datensatz,:), Hochpunkte.SEP_DEC_x(n_datensatz,:)] = max(SEP_DEC(n_datensatz,:));
TAS_SEP_H_DEC_vec(n_datensatz,:) = [v_TAS_j_DEC(n_datensatz, Hochpunkte.SEP_DEC_x(n_datensatz,1)), Hochpunkte.SEP_DEC_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

% fuer Ausweichflug ALT
SEP_CL_ALT(n_datensatz, :) = v_TAS_j_CL_ALT(n_datensatz, :) .* (S_G_vorh_j_CL_ALT(n_datensatz, :) - eps_kompr_j_CL_ALT(n_datensatz, :));
[Hochpunkte.SEP_CL_ALT_y(n_datensatz,:), Hochpunkte.SEP_CL_ALT_x(n_datensatz,:)] = max(SEP_CL_ALT(n_datensatz,:));
TAS_SEP_H_CL_ALT_vec(n_datensatz,:) = [v_TAS_j_CL_ALT(n_datensatz, Hochpunkte.SEP_CL_ALT_x(n_datensatz,1)), Hochpunkte.SEP_CL_ALT_y(n_datensatz,1), hoehe_m(1,n_datensatz)];



%% Spezifische Reichweite SR

SR(n_datensatz,:) = (v_TAS_j(n_datensatz,:) .* specs.g)./(sfc_1PERs_Horizontalflug(n_datensatz,:) .* eps_kompr_j(n_datensatz,:) .* Momentane_Masse_ICA);
[Hochpunkte.SR_y(n_datensatz,:), Hochpunkte.SR_x(n_datensatz,:)] = max(SR(n_datensatz,:));
TAS_SR_H_vec(n_datensatz,:) = [v_TAS_j(n_datensatz, Hochpunkte.SR_x(n_datensatz,1)), Hochpunkte.SR_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SR_CL(n_datensatz,:) = (v_TAS_j_CL(n_datensatz,:) .* specs.g)./(sfc_1PERs_Horizontalflug_CL(n_datensatz,:) .* eps_kompr_j_CL(n_datensatz,:) .* Momentane_Masse_CL);
[Hochpunkte.SR_CL_y(n_datensatz,:), Hochpunkte.SR_CL_x(n_datensatz,:)] = max(SR_CL(n_datensatz,:));
TAS_SR_H_CL_vec(n_datensatz,:) = [v_TAS_j_CL(n_datensatz, Hochpunkte.SR_CL_x(n_datensatz,1)), Hochpunkte.SR_CL_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SR_DEC(n_datensatz,:) = (v_TAS_j_DEC(n_datensatz,:) .* specs.g)./(sfc_1PERs_Horizontalflug_DEC(n_datensatz,:) .* eps_kompr_j_DEC(n_datensatz,:) .* Momentane_Masse_DEC);
[Hochpunkte.SR_DEC_y(n_datensatz,:), Hochpunkte.SR_DEC_x(n_datensatz,:)] = max(SR_DEC(n_datensatz,:));
TAS_SR_H_DEC_vec(n_datensatz,:) = [v_TAS_j_DEC(n_datensatz, Hochpunkte.SR_DEC_x(n_datensatz,1)), Hochpunkte.SR_DEC_y(n_datensatz,1), hoehe_m(1,n_datensatz)];


%% Spezifische Flugdauer SE

SE(n_datensatz,:) = 1 ./ (sfc_1PERs_Horizontalflug(n_datensatz,:) .* eps_kompr_j(n_datensatz,:) .* Masse_ICA); 
[Hochpunkte.SE_y(n_datensatz,:), Hochpunkte.SE_x(n_datensatz,:)] = max(SE(n_datensatz,:));
TAS_SE_H_vec(n_datensatz,:) = [v_TAS_j(n_datensatz, Hochpunkte.SE_x(n_datensatz,1)), Hochpunkte.SE_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SE_CL(n_datensatz,:) = 1 ./ (sfc_1PERs_Horizontalflug_CL(n_datensatz,:) .* eps_kompr_j_CL(n_datensatz,:) .* Masse_CL); 
[Hochpunkte.SE_CL_y(n_datensatz,:), Hochpunkte.SE_CL_x(n_datensatz,:)] = max(SE_CL(n_datensatz,:));
TAS_SE_H_CL_vec(n_datensatz,:) = [v_TAS_j_CL(n_datensatz, Hochpunkte.SE_CL_x(n_datensatz,1)), Hochpunkte.SE_CL_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

% SE_DEC(n_datensatz,:) = 1 ./ (sfc_1PERs_Horizontalflug_DEC(n_datensatz,:) .* eps_kompr_j_DEC(n_datensatz,:) .* Masse_DEC .* specs.g); 
% [Hochpunkte.SE_DEC_y(n_datensatz,:), Hochpunkte.SE_DEC_x(n_datensatz,:)] = max(SE_DEC(n_datensatz,:));
% TAS_SE_H_DEC_vec(n_datensatz,:) = [v_TAS_j_DEC(n_datensatz, Hochpunkte.SE_DEC_x(n_datensatz,1)), Hochpunkte.SE_DEC_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SE_DEC(n_datensatz,:) = 1 ./ (sfc_1PERs_Horizontalflug_DEC(n_datensatz,:) .* eps_kompr_j_DEC(n_datensatz,:) .* Masse_DEC); 
[Hochpunkte.SE_DEC_y(n_datensatz,:), Hochpunkte.SE_DEC_x(n_datensatz,:)] = max(SE_DEC(n_datensatz,:));
TAS_SE_H_DEC_vec(n_datensatz,:) = [v_TAS_j_DEC(n_datensatz, Hochpunkte.SE_DEC_x(n_datensatz,1)), Hochpunkte.SE_DEC_y(n_datensatz,1), hoehe_m(1,n_datensatz)];





%% Fuer Flugbereichsdiagramm
% berechnung v_min v_max aus Horizontalflugdiagramm

intersection_ICA = InterX([v_EAS_j(1,:); S_G_erf_j(n_datensatz,:)],[v_EAS_j(1,:); S_G_vorh_j(n_datensatz,:)]);
v_TAS_HFD(n_datensatz,:) = intersection_ICA(1,:) ./ sqrt(rho_rho0_H(n_datensatz,:));
S_G_inter_HFD(n_datensatz,:) = intersection_ICA(2,:);

intersection_CL = InterX([v_EAS_j_CL(1,:); S_G_erf_j_CL(n_datensatz,:)],[v_EAS_j_CL(1,:); S_G_vorh_j_CL(n_datensatz,:)]);
v_TAS_HFD_CL(n_datensatz,:) = intersection_CL(1,:) ./ sqrt(rho_rho0_H(n_datensatz,:));
S_G_inter_HFD_CL(n_datensatz,:) = intersection_CL(2,:);

intersection_DEC = InterX([v_EAS_j(1,:); S_G_erf_j(n_datensatz,:)],[v_EAS_j(1,:); S_G_vorh_j(n_datensatz,:)]);
v_TAS_HFD_DEC(n_datensatz,:) = intersection_DEC(1,:) ./ sqrt(rho_rho0_H(n_datensatz,:));
S_G_inter_HFD_DEC(n_datensatz,:) = intersection_DEC(2,:);


end



%% Flugbereichsdiagramm CR


% Physikalische Grenzen

% Formel 27 S10
v_s_1g = sqrt((2./(rho_H .* c_A_max)) .* (G ./ Ergebnisse_Fluegel.F));
v_s_min_DEV = 0.94 .* v_s_1g;



m_BOP = 0.795;
M_BO = m_BOP/sqrt(cos(Ergebnisse_Fluegel.phi_25_max));
v_BO = M_BO .* a_H;
v_MO = Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.M_Rumpf.v_D_EAS .* sqrt(ISA.rho_0./ ISA.rho(hoehe_m));


% test v min
intersection = InterX([v_EAS_j(1,:); S_G_erf_j(1,:)],[v_EAS_j(1,:); S_G_vorh_j(1,:)]);
intersection_2 = InterX([v_EAS_j; S_G_erf_j],[v_EAS_j; S_G_vorh_j]);
v_min_HFD = v_TAS_HFD(:,1);
v_max_HFD = v_TAS_HFD(:,2);


%% Dienstgipfelhoehe DEC

% Physikalische Grenzen

% Formel 27 S10
v_s_1g_DEC = sqrt((2./(rho_H .* c_A_max_LDG)) .* (Momentane_Masse_DEC ./ Ergebnisse_Fluegel.F));
v_s_min_DEC = 0.94 .* v_s_1g_DEC;


% Dienstgipfelhoehe: Bedingung SEP < 0.5 [m/s]
zaehler2 = 1;
for z = 1 : length(TAS_SEP_H_DEC_vec)
    if TAS_SEP_H_DEC_vec(z,2) > 0.5
        zaehler2 = z;
    else
    end
end
H_Dienstgipfel_DEC= hoehe_m(1,zaehler2);
% Maximale Kabinendruckhoehe soll 1000m oder 1000ft ueber H_Dienstgipfel
% liegen
H_Kabienendruck_DEC = H_Dienstgipfel_DEC + (unitsratio('m','ft').* 1500); 



m_BOP = 0.795;
M_BO_DEC = m_BOP/sqrt(cos(Ergebnisse_Fluegel.phi_25_max));
v_BO_DEC = M_BO .* a_H;
v_MO_DEC = Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.M_Rumpf.v_D_EAS .* sqrt(ISA.rho_0./ ISA.rho(hoehe_m));


% test v min
intersection_DEC = InterX([v_EAS_j_DEC(1,:); S_G_erf_j_DEC(1,:)],[v_EAS_j_DEC(1,:); S_G_vorh_j_DEC(1,:)]);
intersection_2_DEC = InterX([v_EAS_j_DEC; S_G_erf_j_DEC],[v_EAS_j_DEC; S_G_vorh_j_DEC]);
v_min_HFD_DEC = v_TAS_HFD_DEC(:,1);
v_max_HFD_DEC = v_TAS_HFD_DEC(:,2);



%nachtrag Flugbereichsdiagramm CR
% Dienstgipfelhoehe: Bedingung SEP < 0.5 [m/s]
zaehler1 = 1;
for z = 1 : length(TAS_SEP_H_vec)
    if TAS_SEP_H_vec(z,2) > 0.5
        zaehler1 = zaehler1 +1;
    end
    if zaehler1 >= length(TAS_SEP_H_vec);
       zaehler1 = length(TAS_SEP_H_vec);
    end
end
H_Dienstgipfel_CR = hoehe_m(1,zaehler1);
% Maximale Kabinendruckhoehe soll 1000m oder 1000ft ueber H_Dienstgipfel
% liegen
H_Kabienendruck_CR = H_Dienstgipfel_DEC + (unitsratio('m','ft').* 1500); 





% figure(7)
% hold on
% grid on
% 
% title('Flugbereichsdiagramm', FontSize=25)
% xlabel('v_{TAS} in m/s', FontSize=20);
% ylabel('H in m', FontSize=20);
% 
% pl(1) = plot(v_s_1g, hoehe_m, '--b');
% pl(2) = plot(v_s_min, hoehe_m, 'k');
% pl(3) = plot(v_BO, hoehe_m, 'r');
% 
% pl(4) = plot(TAS_SR_H_vec(:,1), TAS_SR_H_vec(:,3),'.-k');
% pl(5) = plot(TAS_SEP_H_vec(:,1), TAS_SEP_H_vec(:,3),'.-b');
% pl(6) = plot(TAS_SET_H_vec(:,1), TAS_SET_H_vec(:,3),'.-r');
% pl(7) = plot(TAS_SE_H_vec(:,1), TAS_SE_H_vec(:,3),'.-m');
% 
% pl(8) = plot(v_MO, hoehe_m, 'g');
% 
% % 
% pl(9) = plot(v_min_HFD, hoehe_m.', '.-g');
% pl(10) = plot(v_max_HFD, hoehe_m, '.-k');
% 
% pl(11) = plot(specs.Ma_CR * ISA.a(hoehe_CR), hoehe_CR, '*r'); % Design point
% 
% legend('v_s_{1g}', 'v_s_{min}', 'Ma_{BO}', 'SR_{max}', 'SEP_{max}', 'SET_{max}', 'SE_{max}', 'v_{MO}', 'DP', 'FontSize',15 ,Location='eastoutside');

% 
% 


%% Speichern von Ergebnissen


Ergebnisse_Flugleistung_1.hoehe_m = hoehe_m;
Ergebnisse_Flugleistung_1.c_A_j = c_A_j;
Ergebnisse_Flugleistung_1.c_A_F_j = c_A_F_j;
Ergebnisse_Flugleistung_1.c_W_inkomp_j = c_W_inkomp_j;
Ergebnisse_Flugleistung_1.j = j;
Ergebnisse_Flugleistung_1.v_EAS_j = v_EAS_j;
Ergebnisse_Flugleistung_1.v_EAS_j_CL = v_EAS_j_CL;
Ergebnisse_Flugleistung_1.v_EAS_j_DEC = v_EAS_j_DEC;
Ergebnisse_Flugleistung_1.v_TAS_j = v_TAS_j;
Ergebnisse_Flugleistung_1.v_TAS_j_CL = v_TAS_j_CL;
Ergebnisse_Flugleistung_1.v_TAS_j_DEC = v_TAS_j_DEC;
Ergebnisse_Flugleistung_1.v_TAS_j_CL_ALT = v_TAS_j_CL_ALT;
Ergebnisse_Flugleistung_1.Ma_j = Ma_j;
Ergebnisse_Flugleistung_1.Ma_j_CL = Ma_j_CL;
Ergebnisse_Flugleistung_1.Ma_J_DEC = Ma_j_DEC;
Ergebnisse_Flugleistung_1.Ma_j_CL_ALT = Ma_j_CL_ALT;
Ergebnisse_Flugleistung_1.S_S0_KF_j_CR = S_S0_KF_j_CR;
Ergebnisse_Flugleistung_1.S_S0_KF_j_CL = S_S0_KF_j_CL;
Ergebnisse_Flugleistung_1.S_S0_KF_j_DEC = S_S0_KF_j_DEC;
Ergebnisse_Flugleistung_1.S_S0_KF_j_CL_ALT = S_S0_KF_j_CL_ALT;
Ergebnisse_Flugleistung_1.S_S0_E = S_S0_E;
Ergebnisse_Flugleistung_1.S_S0_j_CR = S_S0_j;
Ergebnisse_Flugleistung_1.S_S0_j_CL = S_S0_j_CL;
Ergebnisse_Flugleistung_1.S_S0_j_DEC = S_S0_j_DEC;
Ergebnisse_Flugleistung_1.S_S0_j_CL_ALT = S_S0_j_CL_ALT;
Ergebnisse_Flugleistung_1.S_G_vorh_j = S_G_vorh_j;
Ergebnisse_Flugleistung_1.S_G_vorh_j_CL = S_G_vorh_j_CL;
Ergebnisse_Flugleistung_1.S_G_vorh_j_DEC = S_G_vorh_j_DEC;
Ergebnisse_Flugleistung_1.sfc_1PERs_Horizontalflug = sfc_1PERs_Horizontalflug;
Ergebnisse_Flugleistung_1.sfc_1PERs_Horizontalflug_CL = sfc_1PERs_Horizontalflug_CL;
Ergebnisse_Flugleistung_1.sfc_1PERs_Horizontalflug_DEC = sfc_1PERs_Horizontalflug_DEC;
Ergebnisse_Flugleistung_1.sfc_daNh_CR = sfc_daNh_CR;
Ergebnisse_Flugleistung_1.sfc_daNh_CL = sfc_daNh_CL;
Ergebnisse_Flugleistung_1.sfc_daNh_DEC = sfc_daNh_DEC;
Ergebnisse_Flugleistung_1.c_W_j = c_W_j;
Ergebnisse_Flugleistung_1.c_W_j_CL = c_W_j_CL;
Ergebnisse_Flugleistung_1.c_W_j_DEC = c_W_j_DEC;
Ergebnisse_Flugleistung_1.S_G_erf_j = S_G_erf_j; 
Ergebnisse_Flugleistung_1.S_G_erf_j_CL = S_G_erf_j_CL;
Ergebnisse_Flugleistung_1.S_G_erf_j_DEC = S_G_erf_j_DEC;
Ergebnisse_Flugleistung_1.SET = SET;
Ergebnisse_Flugleistung_1.TAS_SET_H_vec = TAS_SET_H_vec;
Ergebnisse_Flugleistung_1.SET_CL = SET_CL;
Ergebnisse_Flugleistung_1.TAS_SET_H_CL_vec = TAS_SE_H_CL_vec;
Ergebnisse_Flugleistung_1.SET_DEC = SET_DEC;
Ergebnisse_Flugleistung_1.TAS_SET_H_DEC_vec = TAS_SET_H_DEC_vec;
Ergebnisse_Flugleistung_1.SEP = SEP;
Ergebnisse_Flugleistung_1.TAS_SEP_H_vec = TAS_SEP_H_vec;
Ergebnisse_Flugleistung_1.SEP_CL = SEP_CL;
Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec = TAS_SEP_H_CL_vec;
Ergebnisse_Flugleistung_1.SEP_DEC = SEP_DEC;
Ergebnisse_Flugleistung_1.TAS_SEP_H_DEC_vec = TAS_SEP_H_DEC_vec;
Ergebnisse_Flugleistung_1.SR = SR;
Ergebnisse_Flugleistung_1.TAS_SR_H_vec = TAS_SR_H_vec;
Ergebnisse_Flugleistung_1.SR_CL = SR_CL;
Ergebnisse_Flugleistung_1.TAS_SR_H_CL_vec = TAS_SR_H_CL_vec;
Ergebnisse_Flugleistung_1.SR_DEC = SR_DEC;
Ergebnisse_Flugleistung_1.TAS_SR_H_DEC_vec = TAS_SR_H_DEC_vec;
Ergebnisse_Flugleistung_1.SE = SE;
Ergebnisse_Flugleistung_1.TAS_SE_H_vec = TAS_SE_H_vec;
Ergebnisse_Flugleistung_1.SE_CL = SE_CL; 
Ergebnisse_Flugleistung_1.TAS_SE_H_CL_vec = TAS_SE_H_CL_vec;
Ergebnisse_Flugleistung_1.SE_DEC = SE_DEC;
Ergebnisse_Flugleistung_1.TAS_SE_H_DEC_vec = TAS_SE_H_DEC_vec;
Ergebnisse_Flugleistung_1.Hochpunkte = Hochpunkte;
Ergebnisse_Flugleistung_1.v_TAS_HFD = v_TAS_HFD; 
Ergebnisse_Flugleistung_1.v_TAS_HFD_CL = v_TAS_HFD_CL;
Ergebnisse_Flugleistung_1.v_TAS_HFD_DEC = v_TAS_HFD_DEC;
Ergebnisse_Flugleistung_1.S_G_inter_HFD = S_G_inter_HFD;
Ergebnisse_Flugleistung_1.S_G_inter_HFD_CL = S_G_inter_HFD_CL;
Ergebnisse_Flugleistung_1.S_G_inter_HFD_DEC = S_G_inter_HFD_DEC;
Ergebnisse_Flugleistung_1.SEP_CL_ALT = SEP_CL_ALT;
Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec = TAS_SEP_H_CL_ALT_vec;

% fuer Flugbereichsdiagramm
        % CR
Ergebnisse_Flugleistung_1.v_s_1g = v_s_1g;
Ergebnisse_Flugleistung_1.v_s_min = v_s_min_DEV;
Ergebnisse_Flugleistung_1.m_BOP = m_BOP;
Ergebnisse_Flugleistung_1.M_BO = M_BO;
Ergebnisse_Flugleistung_1.v_BO = v_BO;
Ergebnisse_Flugleistung_1.v_MO = v_MO;
Ergebnisse_Flugleistung_1.v_min_HFD = v_min_HFD;
Ergebnisse_Flugleistung_1.v_max_HFD = v_max_HFD;
Ergebnisse_Flugleistung_1.H_Dienstgipfel_CR = H_Dienstgipfel_CR;
Ergebnisse_Flugleistung_1.H_Kabienendruck_CR = H_Kabienendruck_CR;
        % DEC
Ergebnisse_Flugleistung_1.v_s_1g_DEC = v_s_1g_DEC;
Ergebnisse_Flugleistung_1.v_s_min_DEC = v_s_min_DEC;
Ergebnisse_Flugleistung_1.M_BO_DEC = M_BO_DEC;
Ergebnisse_Flugleistung_1.v_BO_DEC = v_BO_DEC;
Ergebnisse_Flugleistung_1.v_MO_DEC = v_MO_DEC;
Ergebnisse_Flugleistung_1.v_min_HFD_DEC = v_min_HFD_DEC;
Ergebnisse_Flugleistung_1.v_max_HFD_DEC = v_max_HFD_DEC;
Ergebnisse_Flugleistung_1.H_Dienstgipfel_DEC = H_Dienstgipfel_CR;
Ergebnisse_Flugleistung_1.H_Kabienendruck_DEC = H_Kabienendruck_DEC;



save Ergebnisse_FLugleistung_1.mat Ergebnisse_Flugleistung_1

end





%% Funktion zur Berechnung der Atmospherischen Daten
function [rho_rho0_H, T_H, a_H, p_p0_H, rho] = Atmos_H(hoehe_m)
    load Ergebnisse_ISA_DATA.mat
    for zv = 1:size(hoehe_m)
    rho_rho0_H(zv,1) = ISA.rho(hoehe_m(zv,1)) / ISA.rho_0;
    T_H(zv,1) = ISA.T(hoehe_m(zv,1));
    a_H(zv,1)= ISA.a(hoehe_m(zv,1));
    p_p0_H(zv,1) = ISA.p(hoehe_m(zv,1)) / ISA.p0;
    rho(zv,1) = ISA.rho(hoehe_m(zv,1));
    end
end



%% Funktion zur berechnung von Gleichung 12 PS7 s.3 
function [S_S0_KF_j] = S_S0_KF_j(D, rho_rho0_H, Ma_j, p_p0, bypass)
    
    S_S0_KF_j = D .* rho_rho0_H .* exp(-0.35 .* Ma_j .* p_p0 .* sqrt(bypass));
end


%% 

%% PLOTS

% 
% % figure(1) % Horiontalflugdiagramm
% % hold on
% % grid on
% % title('Horizontalflugdiagramm G_{ICA}', 'FontSize',25);
% % ylabel('S/G', 'FontSize',20);
% % xlabel('v_{EAS}', 'FontSize',20);
% % xlim([0 300]);
% % ylim([0 1]);
% % 
% % p1_1(n_datensatz) = plot(v_EAS_j, S_G_vorh_j(n_datensatz, :), 'Color',colors_1(n_datensatz,:), 'LineStyle','--');
% % p1_2(n_datensatz) = plot(v_EAS_j, S_G_erf_j(n_datensatz, :),  'Color',colors_2(n_datensatz,:), 'LineStyle','-');
% % Legend{2*n_datensatz-1} =  '(S/G)_{erf} Hoehe ' + sprintf("%d",hoehe_m(n_datensatz)) + 'm' ;
% % Legend{2*n_datensatz} =  '(S/G)_{vorh} Hoehe ' + sprintf("%d",hoehe_m(n_datensatz)) + 'm';
% % legend(Legend(1:n_datensatz*2),Location='eastoutside', FontSize=15);
% % 
% % % plot(v_EAS_j(n_datensatz, :), eps_inkomp_j(n_datensatz, :), '--b');
% % hold off
% 
% %% Plot Luca
% %  P = InterX([v_EAS_j(n_datensatz, :);S_G_erf_j(n_datensatz, :)],[v_EAS_j(n_datensatz, :);S_G_vorh_j(n_datensatz, :)]);
% % 
% %     if(size(P,2)  == 1 )
% %         v_TAS_max(n_datensatz, :) = P(1,1) .* (sqrt(rho_rho0_H));
% %         v_TAS_min(n_datensatz, :) = NaN;
% %     elseif(size(P,2)  == 0 )
% %         v_TAS_max(n_datensatz, :) = NaN;
% %         v_TAS_min(n_datensatz, :) = NaN;
% %     else
% %         v_TAS_max(n_datensatz, :) = P(1,2) .* (sqrt(rho_rho0_H));
% %         v_TAS_min(n_datensatz, :) = P(1,1) .* (sqrt(rho_rho0_H));
% %     end
% % 
% % 
% % 
% %      figure(1)
% %     grid on 
% %     hold on
% %     xlim([0 350])
% %     ylim([0 0.6])
% %     
% % %     if((Step_Plot ~= 0 && mod(hoehe(n_datensatz),Step_Plot) == 0))
% %         k = k+1;
% %         plot(v_EAS_j(n_datensatz, :),S_G_erf_j(n_datensatz, :));%,"Color",colorgradient(n_datensatz,length(hoehe_m)))
% %         plot(v_EAS_j(n_datensatz, :),S_G_vorh_j(n_datensatz, :));%,"Color",colorgradient(n_datensatz,length(hoehe_m),LineStyle='--'))
% %         Legend{2*k-1} =  '(S/G)_{erf} FL' + sprintf("%d",hoehe(n_datensatz)/100) ;
% %         Legend{2*k} =  '(S/G)_{vorh} FL' + sprintf("%d",hoehe(n_datensatz)/100);
% %         
% %         
% % %     end
% % 
% % 
% % if(Step_Plot ~=0)
% %  legend(Legend(1:k*2),Location='eastoutside');
% % end
% % 
% % hold off
% 
% % figure(2) % SET
% % hold on 
% % grid on
% % title('Spezifischer Schubüberschuss G_{ICA}', 'FontSize',25)
% % ylabel('SET','FontSize',20);
% % xlabel('v_{TAS}','FontSize',20);
% % xlim([0 350]);
% % ylim([0 0.7]);
% % 
% % p2_1(n_datensatz) = plot(v_TAS_j(n_datensatz, :), SET(n_datensatz, :), 'Color',colors_2(n_datensatz,:), 'LineStyle','-');
% % p2_2(n_datensatz) = plot(v_TAS_j(n_datensatz, SET_x(n_datensatz,1)), SET_y(n_datensatz,1), '*k');
% % Legend2{n_datensatz*2-1} = 'SET in [rad] Hoehe ' + sprintf("%d",hoehe_m(n_datensatz)) + 'm';
% % Legend2{n_datensatz*2} = '';
% % legend(Legend2(1:n_datensatz*2), Location='eastoutside', FontSize=16)
% % 
% % hold off
% 
% % 
% % figure(3) % SEP
% % hold on 
% % grid on
% % title('Spezifischer Leistungsüberschuss G_{ICA}','FontSize',25)
% % ylabel('SEP in m/s','FontSize',20);
% % xlabel('v_{TAS} in m/s','FontSize',20);
% % xlim([0 350]);
% % ylim([0 100]);
% % 
% % p3_1(n_datensatz) = plot(v_TAS_j(n_datensatz, :), SEP(n_datensatz, :), 'Color',colors_2(n_datensatz,:), 'LineStyle','-');
% % p3_2(n_datensatz) = plot(v_TAS_j(n_datensatz, SEP_x(n_datensatz,1)), SEP_y(n_datensatz,1), '*b');
% % Legend3{n_datensatz*2-1} = 'SEP in m/s Hoehe ' + sprintf("%d",hoehe_m(n_datensatz)) + 'm' ;
% % Legend3{n_datensatz*2} = '';
% % legend(Legend3(1:n_datensatz*2), Location='eastoutside', FontSize=16)
% % 
% % hold off
% 
% 
% % figure(4) % SR
% % hold on 
% % grid on
% % title('Spezifische Reichweite G_{ICA}','FontSize',25)
% % ylabel('SR in m/kg','FontSize',20);
% % xlabel('v_{TAS} in m/s','FontSize',20);
% % xlim([0 300]);
% % ylim([0 20]);
% % 
% % p4_1(n_datensatz) = plot(v_TAS_j(n_datensatz, :), SR(n_datensatz, :), 'Color',colors_2(n_datensatz,:), 'LineStyle','-');
% % p4_2(n_datensatz) = plot(v_TAS_j(n_datensatz, SR_x(n_datensatz,1)), SR_y(n_datensatz,1), '*k');
% % Legend4{n_datensatz*2-1} = 'SEP in m/kg Hoehe ' + sprintf("%d",hoehe_m(n_datensatz)) + 'm' ;
% % Legend4{n_datensatz*2} = '';
% % legend(Legend4(1:n_datensatz*2), Location='eastoutside', FontSize=16)
% % hold off
% % 
% % 
% % 
% % 
% % % Spezifische Flugdauer
% % figure(6) 
% % 
% % hold on 
% % grid on 
% % title('Spezifische Flugdauer G_{ICA}','FontSize',25)
% % ylabel('SE in s/kg','FontSize',20);
% % xlabel('v_{TAS} in m/s','FontSize',20);
% % xlim([0 350]);
% % 
% % 
% % p5_1(n_datensatz) = plot(v_TAS_j(n_datensatz, :) , SE(n_datensatz,:), 'Color',colors_2(n_datensatz,:), 'LineStyle','-' );
% % p5_2(n_datensatz) = plot(v_TAS_j(n_datensatz, SE_x(n_datensatz,1)), SE_y(n_datensatz,1), '*k');
% % Legend5{n_datensatz*2-1} = ' in s/kg Hoehe ' + sprintf("%d",hoehe_m(n_datensatz)) + 'm' ;
% % Legend5{n_datensatz*2} = '';
% % legend(Legend5(1:n_datensatz*2), Location='eastoutside', FontSize=16)
% % hold off
% 
% % figure(5)
% hold on 
% grid on 
% 
% title('Optimalgeschwindigkeiten bei G_{LDG}','FontSize',25)
% 
% ylabel('H in m','FontSize',20);
% xlabel('v_{TAS} in m/s','FontSize',20);
% xlim([0 350]);
% 
% plot(TAS_SR_H_DEC_vec(1:10,1), TAS_SR_H_DEC_vec(1:10,3),'k');
% plot(TAS_SEP_H_DEC_vec(1:10,1), TAS_SEP_H_DEC_vec(1:10,3),'b');
% plot(TAS_SET_H_DEC_vec(1:10,1), TAS_SET_H_DEC_vec(1:10,3),'r');
% legend('SR', 'SEP', 'SET', 'FontSize',15)

