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
addpath('Unterfunktionen Widerstand');

%% Annahmen

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
Momentane_Masse_DEC = (Ergebnisse_Massen_FE2.M_TO - Ergebnisse_Massen_FE2.M_F)*specs.g;
Momentane_Masse_CL = TO_Masse * specs.g;
GTO_G_ICA = G_TO / Momentane_Masse_ICA;
GTO_G_CL = 1;
GTO_G_DEC = G_TO / Momentane_Masse_DEC;
S0 = k_CR * G_TO * (1/schub_CR.Eta) / (schub_CR.S_S0_CR * schub_CR.S_S0_E);
S0_GTO = S0/G_TO;

% Annahmen fuer Flugbereichsdiagramm

c_A_max = 3.58; %%%%%%%%%%%%%%% Achtung random wert, bitte von mac geben lassen
G = Momentane_Masse_ICA; %%%%%%%%%%%%%%%%%%% ACHTUNG


hoehe_CR = round(unitsratio('m','ft')*(specs.flight_level*10^2));

%% Horizontalflugdiagramm

% Flughoehe = specs.flight_level * 10^2 ;                         % in ft
% hoehe =round(unitsratio('m','ft')*Flughoehe);                     % in m

schritte = 10;
hoehe_plus = 100;
hoehe = (round(linspace(5, ((specs.flight_level + hoehe_plus)*10^2), schritte))); % 1000: 3000: specs.flight_level*10^2;
% hoehe = 1000: 1000: specs.flight_level*10^2;
hoehe_m = round(unitsratio('m','ft')*hoehe);
hoehe_m(1,length(hoehe_m)+1) = hoehe_CR;
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

[rho_rho0_H, T_H, a_H, p_p0_H, rho] = Atmos_H(hoehe_m.');


% PS7 S2 Formel 3
v_EAS_j = sqrt(( (2)./ (ISA.rho_0 .* c_A_j) ) .* (Momentane_Masse_ICA ./Ergebnisse_Fluegel.F));
v_EAS_j_CL = sqrt(( (2) ./ (ISA.rho_0 .* c_A_j) ) .* (Momentane_Masse_CL ./Ergebnisse_Fluegel.F));
v_EAS_j_DEC = sqrt(( (2) ./ (ISA.rho_0 .* c_A_j) ) .* (Momentane_Masse_DEC./Ergebnisse_Fluegel.F));




% PS7 S2 Formel 4 Flugmachzahl Ma_j in TAS
v_TAS_j = (v_EAS_j)./(sqrt(rho_rho0_H));
v_TAS_j_CL = (v_EAS_j_CL)./(sqrt(rho_rho0_H));
v_TAS_j_DEC = (v_EAS_j_DEC)./(sqrt(rho_rho0_H));

Ma_j = v_TAS_j ./ a_H;
Ma_j_CL = v_TAS_j_CL ./ a_H;
Ma_j_DEC = v_TAS_j_DEC ./ a_H;


% Vorhandenes Schub-Gewichts-Verhaeltnis

% PS7 S3 Formel 12
S_S0_KF_j_CR = S_S0_KF_j(specs.Drosselgrad(1,3), rho_rho0_H, Ma_j, p_p0_H, specs.bypass);
S_S0_KF_j_CL = S_S0_KF_j(specs.Drosselgrad(1,2), rho_rho0_H, Ma_j_CL, p_p0_H, specs.bypass);
S_S0_KF_j_DEC = S_S0_KF_j(specs.Drosselgrad(1,3), rho_rho0_H, Ma_j_DEC, p_p0_H, specs.bypass);



% Einlaufverluste PS7 S3 Formel 13
dp_p0 = 0.02;
S_S0_E = 1 - (1.3 + 0.25 * specs.bypass) * dp_p0; 

 
% PS7 S4 Formel 14
S_S0_j = S_S0_KF_j_CR + S_S0_E;
S_S0_j_CL = S_S0_KF_j_CL + S_S0_E;
S_S0_j_DEC = S_S0_KF_j_DEC + S_S0_E;


% PS7 S4 Formel 15

S_G_vorh_j = S_S0_j .* S0_GTO .* GTO_G_ICA;
S_G_vorh_j_CL = S_S0_j_CL .* S0_GTO .* GTO_G_CL;
S_G_vorh_j_DEC = S_S0_j_DEC .* S0_GTO .* GTO_G_DEC;

% SFC
[sfc_daNh_CR, ~, sfc_1PERs_Horizontalflug] = SFC_vec(hoehe_m, Ma_j, specs.bypass);
[sfc_daNh_CL, ~, sfc_1PERs_Horizontalflug_CL] = SFC_vec(hoehe_m, Ma_j_CL, specs.bypass);
[sfc_daNh_DEC, ~, sfc_1PERs_Horizontalflug_DEC] = SFC_vec(hoehe_m, Ma_j_DEC, specs.bypass);

for n_datensatz = 1:length(hoehe_m)


% Transsonischer Widerstand mit Funktion aus PS4
for zv = 1:j
delta_c_WM(n_datensatz,zv) = Transsonischer_W(Ma_j(n_datensatz,zv), c_A_F_j(1,zv));
delta_c_WM_CL(n_datensatz,zv) = Transsonischer_W(Ma_j_CL(n_datensatz,zv), c_A_F_j(1,zv));
delta_c_WM_DEC(n_datensatz,zv) = Transsonischer_W(Ma_j_DEC(n_datensatz,zv), c_A_F_j(1,zv));
end

% Gesamtwiderstansbeiwert PS7 S3 Formel 10
c_W_j(n_datensatz, :) = c_W_inkomp_j + delta_c_WM(n_datensatz, :);
c_W_j_CL(n_datensatz, :) = c_W_inkomp_j + delta_c_WM_CL(n_datensatz, :);
c_W_j_DEC(n_datensatz, :) = c_W_inkomp_j + delta_c_WM_DEC(n_datensatz, :);


% PS7 S3 Formel 11 Gleitzahl kompressibel
eps_kompr_j(n_datensatz, :) = c_W_j(n_datensatz, :) ./ c_A_j; % eps_kompr_j = (S/G)_erf_j
eps_kompr_j_CL(n_datensatz, :) = c_W_j_CL(n_datensatz, :) ./ c_A_j;
eps_kompr_j_DEC(n_datensatz, :) = c_W_j_DEC(n_datensatz, :) ./ c_A_j;

S_G_erf_j(n_datensatz, :) = c_W_j(n_datensatz, :) ./ c_A_j;
S_G_erf_j_CL(n_datensatz, :) = c_W_j_CL(n_datensatz, :) ./ c_A_j;
S_G_erf_j_DEC(n_datensatz, :) = c_W_j_DEC(n_datensatz, :) ./ c_A_j;




%% Optimale Leistungszustände
%% Spez Schubueberschuss SET

SET(n_datensatz, :) = S_G_vorh_j(n_datensatz, :) - eps_kompr_j(n_datensatz, :);
[SET_y(n_datensatz,:), SET_x(n_datensatz,:)] = max(SET(n_datensatz,:));
TAS_SET_H_vec(n_datensatz,:) = [v_TAS_j(n_datensatz, SET_x(n_datensatz,1)), SET_y(n_datensatz,1), hoehe_m(1,n_datensatz)];


SET_CL(n_datensatz, :) = S_G_vorh_j_CL(n_datensatz, :) - eps_kompr_j_CL(n_datensatz, :);
[SET_CL_y(n_datensatz,:), SET_CL_x(n_datensatz,:)] = max(SET_CL(n_datensatz,:));
TAS_SET_H_CL_vec(n_datensatz,:) = [v_TAS_j_CL(n_datensatz, SET_CL_x(n_datensatz,1)), SET_CL_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SET_DEC(n_datensatz, :) = S_G_vorh_j_DEC(n_datensatz, :) - eps_kompr_j_DEC(n_datensatz, :);
[SET_DEC_y(n_datensatz,:), SET_DEC_x(n_datensatz,:)] = max(SET_DEC(n_datensatz,:));
TAS_SET_H_DEC_vec(n_datensatz,:) = [v_TAS_j_DEC(n_datensatz, SET_DEC_x(n_datensatz,1)), SET_DEC_y(n_datensatz,1), hoehe_m(1,n_datensatz)];




%% Spezifischer Leistungueberschuss SEP

SEP(n_datensatz, :) = v_TAS_j(n_datensatz, :) .* (S_G_vorh_j(n_datensatz, :) - eps_kompr_j(n_datensatz, :));
[SEP_y(n_datensatz,:), SEP_x(n_datensatz,:)] = max(SEP(n_datensatz,:));
TAS_SEP_H_vec(n_datensatz,:) = [v_TAS_j(n_datensatz, SEP_x(n_datensatz,1)), SEP_y(n_datensatz,1), hoehe_m(1,n_datensatz)];


SEP_CL(n_datensatz, :) = v_TAS_j_CL(n_datensatz, :) .* (S_G_vorh_j_CL(n_datensatz, :) - eps_kompr_j_CL(n_datensatz, :));


[SEP_CL_y(n_datensatz,:), SEP_CL_x(n_datensatz,:)] = max(SEP_CL(n_datensatz,:));
TAS_SEP_H_CL_vec(n_datensatz,:) = [v_TAS_j_CL(n_datensatz, SEP_CL_x(n_datensatz,1)), SEP_CL_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SEP_DEC(n_datensatz, :) = v_TAS_j_DEC(n_datensatz, :) .* (S_G_vorh_j_DEC(n_datensatz, :) - eps_kompr_j_DEC(n_datensatz, :));
[SEP_DEC_y(n_datensatz,:), SEP_DEC_x(n_datensatz,:)] = max(SEP_DEC(n_datensatz,:));
TAS_SEP_H_DEC_vec(n_datensatz,:) = [v_TAS_j_DEC(n_datensatz, SEP_DEC_x(n_datensatz,1)), SEP_DEC_y(n_datensatz,1), hoehe_m(1,n_datensatz)];


%% Spezifische Reichweite SR

SR(n_datensatz,:) = (v_TAS_j(n_datensatz,:))./(sfc_1PERs_Horizontalflug(n_datensatz,:) .* eps_kompr_j(n_datensatz,:) .* Momentane_Masse_ICA);
[SR_y(n_datensatz,:), SR_x(n_datensatz,:)] = max(SR(n_datensatz,:));
TAS_SR_H_vec(n_datensatz,:) = [v_TAS_j(n_datensatz, SR_x(n_datensatz,1)), SR_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SR_CL(n_datensatz,:) = (v_TAS_j_CL(n_datensatz,:))./(sfc_1PERs_Horizontalflug_CL(n_datensatz,:) .* eps_kompr_j_CL(n_datensatz,:) .* Momentane_Masse_CL);
[SR_CL_y(n_datensatz,:), SR_CL_x(n_datensatz,:)] = max(SR_CL(n_datensatz,:));
TAS_SR_H_CL_vec(n_datensatz,:) = [v_TAS_j_CL(n_datensatz, SR_CL_x(n_datensatz,1)), SR_CL_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SR_DEC(n_datensatz,:) = (v_TAS_j_DEC(n_datensatz,:))./(sfc_1PERs_Horizontalflug_DEC(n_datensatz,:) .* eps_kompr_j_DEC(n_datensatz,:) .* Momentane_Masse_DEC);
[SR_DEC_y(n_datensatz,:), SR_DEC_x(n_datensatz,:)] = max(SR_DEC(n_datensatz,:));
TAS_SR_H_DEC_vec(n_datensatz,:) = [v_TAS_j_DEC(n_datensatz, SR_DEC_x(n_datensatz,1)), SR_DEC_y(n_datensatz,1), hoehe_m(1,n_datensatz)];


%% Spezifische Flugdauer SE

SE(n_datensatz,:) = 1 ./ (sfc_1PERs_Horizontalflug(n_datensatz,:) .* eps_kompr_j(n_datensatz,:) .* Momentane_Masse_ICA .* specs.g); 
[SE_y(n_datensatz,:), SE_x(n_datensatz,:)] = max(SE(n_datensatz,:));
TAS_SE_H_vec(n_datensatz,:) = [v_TAS_j(n_datensatz, SE_x(n_datensatz,1)), SE_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SE_CL(n_datensatz,:) = 1 ./ (sfc_1PERs_Horizontalflug_CL(n_datensatz,:) .* eps_kompr_j_CL(n_datensatz,:) .* Momentane_Masse_CL .* specs.g); 
[SE_CL_y(n_datensatz,:), SE_CL_x(n_datensatz,:)] = max(SE_CL(n_datensatz,:));
TAS_SE_H_CL_vec(n_datensatz,:) = [v_TAS_j_CL(n_datensatz, SE_CL_x(n_datensatz,1)), SE_CL_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

SE_DEC(n_datensatz,:) = 1 ./ (sfc_1PERs_Horizontalflug_DEC(n_datensatz,:) .* eps_kompr_j_DEC(n_datensatz,:) .* Momentane_Masse_DEC .* specs.g); 
[SE_DEC_y(n_datensatz,:), SE_DEC_x(n_datensatz,:)] = max(SE_DEC(n_datensatz,:));
TAS_SE_H_DEC_vec(n_datensatz,:) = [v_TAS_j_DEC(n_datensatz, SE_DEC_x(n_datensatz,1)), SE_DEC_y(n_datensatz,1), hoehe_m(1,n_datensatz)];

%% Fuer Flugbereichsdiagramm
% berechnung v_min v_max aus Horizontalflugdiagramm

% intersection_ICA = InterX([v_EAS_j(1,:); S_G_erf_j(n_datensatz,:)],[v_EAS_j(1,:); S_G_vorh_j(n_datensatz,:)]);
% v_TAS_HFD(n_datensatz,:) = intersection_ICA(1,:) ./ sqrt(rho_rho0_H(n_datensatz,:));
% S_G_inter_HFD(n_datensatz,:) = intersection_ICA(2,:);

% intersection_CL = InterX([v_EAS_j_CL(1,:); S_G_erf_j_CL(n_datensatz,:)],[v_EAS_j_CL(1,:); S_G_vorh_j_CL(n_datensatz,:)]);
% v_TAS_HFD_CL(n_datensatz,:) = intersection_CL(1,:) ./ sqrt(rho_rho0_H(n_datensatz,:));
% S_G_inter_HFD_CL(n_datensatz,:) = intersection_CL(2,:);
% 
% intersection_DEC = InterX([v_EAS_j(1,:); S_G_erf_j(n_datensatz,:)],[v_EAS_j(1,:); S_G_vorh_j(n_datensatz,:)]);
% v_TAS_HFD_DEC(n_datensatz,:) = intersection_DEC(1,:) ./ sqrt(rho_rho0_H(n_datensatz,:));
% S_G_inter_HFD_DEC(n_datensatz,:) = intersection_DEC(2,:);


%% PLOTS
figure(1) % Horiontalflugdiagramm
hold on
grid on

ylabel('S/G');
xlabel('v_{EAS}');
xlim([0 300]);
ylim([0 1]);

plot(v_EAS_j, S_G_vorh_j(n_datensatz, :), 'r');
plot(v_EAS_j, S_G_erf_j(n_datensatz, :), 'b');
% plot(v_EAS_j(n_datensatz, :), eps_inkomp_j(n_datensatz, :), '--b');
hold off

%% Plot Luca
%  P = InterX([v_EAS_j(n_datensatz, :);S_G_erf_j(n_datensatz, :)],[v_EAS_j(n_datensatz, :);S_G_vorh_j(n_datensatz, :)]);
% 
%     if(size(P,2)  == 1 )
%         v_TAS_max(n_datensatz, :) = P(1,1) .* (sqrt(rho_rho0_H));
%         v_TAS_min(n_datensatz, :) = NaN;
%     elseif(size(P,2)  == 0 )
%         v_TAS_max(n_datensatz, :) = NaN;
%         v_TAS_min(n_datensatz, :) = NaN;
%     else
%         v_TAS_max(n_datensatz, :) = P(1,2) .* (sqrt(rho_rho0_H));
%         v_TAS_min(n_datensatz, :) = P(1,1) .* (sqrt(rho_rho0_H));
%     end
% 
% 
% 
%      figure(1)
%     grid on 
%     hold on
%     xlim([0 350])
%     ylim([0 0.6])
%     
% %     if((Step_Plot ~= 0 && mod(hoehe(n_datensatz),Step_Plot) == 0))
%         k = k+1;
%         plot(v_EAS_j(n_datensatz, :),S_G_erf_j(n_datensatz, :));%,"Color",colorgradient(n_datensatz,length(hoehe_m)))
%         plot(v_EAS_j(n_datensatz, :),S_G_vorh_j(n_datensatz, :));%,"Color",colorgradient(n_datensatz,length(hoehe_m),LineStyle='--'))
%         Legend{2*k-1} =  '(S/G)_{erf} FL' + sprintf("%d",hoehe(n_datensatz)/100) ;
%         Legend{2*k} =  '(S/G)_{vorh} FL' + sprintf("%d",hoehe(n_datensatz)/100);
%         
%         
% %     end
% 
% 
% if(Step_Plot ~=0)
%  legend(Legend(1:k*2),Location='eastoutside');
% end
% 
% hold off

% figure(2) % SET
% hold on 
% grid on
% ylabel('SET');
% xlabel('v_{TAS}');
% xlim([0 350]);
% ylim([0 1]);
% 
% plot(v_TAS_j(n_datensatz, :), SET(n_datensatz, :), 'r');
% plot(v_TAS_j(n_datensatz, SET_x(n_datensatz,1)), SET_y(n_datensatz,1), '*b');
% hold off

% 
% figure(3) % SEP
% hold on 
% grid on
% ylabel('SEP');
% xlabel('v_{TAS}');
% xlim([0 350]);
% ylim([0 100]);
% 
% plot(v_TAS_j(n_datensatz, :), SEP(n_datensatz, :), 'r');
% plot(v_TAS_j(n_datensatz, SEP_x(n_datensatz,1)), SEP_y(n_datensatz,1), '*b');
% hold off


% figure(4) % SEP
% hold on 
% grid on
% ylabel('SR');
% xlabel('v_{TAS}');
% xlim([0 350]);
% ylim([0 300]);
% 
% plot(v_TAS_j(n_datensatz, :), SR(n_datensatz, :), 'r');
% plot(v_TAS_j(n_datensatz, SR_x(n_datensatz,1)), SR_y(n_datensatz,1), '*b');
% hold off




% Spezifische Flugdauer
% figure(6) 
% 
% hold on 
% grid on 
% 
% ylabel('SE in s/kg');
% xlabel('v_{TAS} in m/s');
% xlim([0 350]);
% 
% 
% plot(v_TAS_j(n_datensatz, :) , SE(n_datensatz,:), 'k');



end
% figure(5)
% hold on 
% grid on 
% 
% ylabel('H in m');
% xlabel('v_{TAS} in m/s');
% xlim([0 350]);
% 
% plot(TAS_SR_H_vec(:,1), TAS_SR_H_vec(:,3),'k');
% plot(TAS_SEP_H_vec(:,1), TAS_SEP_H_vec(:,3),'b');
% plot(TAS_SET_H_vec(:,1), TAS_SET_H_vec(:,3),'r');
% 

%% Flugbereichsdiagramm


% Physikalische Grenzen

% Formel 27 S10
v_s_1g = sqrt((2./(rho .* c_A_max)) .* (G ./ Ergebnisse_Fluegel.F));
v_s_min = 0.94 .* v_s_1g;

m_BOP = 0.795;
M_BO = m_BOP/sqrt(cos(Ergebnisse_Fluegel.phi_25_max));
v_BO = M_BO .* a_H;
v_MO = specs.Ma_MO .* a_H;


% test v min
% intersection = InterX([v_EAS_j(1,:); S_G_erf_j(1,:)],[v_EAS_j(1,:); S_G_vorh_j(1,:)])
% intersection_2 = InterX([v_EAS_j; S_G_erf_j],[v_EAS_j; S_G_vorh_j])
% v_min_HFD = v_TAS_HFD(:,1);
% v_max_HFD = v_TAS_HFD(:,2);



figure(7)
hold on
grid on

pl(1) = plot(v_s_1g, hoehe_m, '--b');
pl(2) = plot(v_s_min, hoehe_m, 'k');
pl(3) = plot(v_BO, hoehe_m, 'r');

pl(4) = plot(TAS_SR_H_vec(:,1), TAS_SR_H_vec(:,3),'.-k');
pl(5) = plot(TAS_SEP_H_vec(:,1), TAS_SEP_H_vec(:,3),'.-b');
pl(6) = plot(TAS_SET_H_vec(:,1), TAS_SET_H_vec(:,3),'.-r');
pl(7) = plot(TAS_SE_H_vec(:,1), TAS_SE_H_vec(:,3),'.-m');

pl(8) = plot(v_MO, hoehe_m, 'g');

% 
% pl(9) = plot(v_min_HFD, hoehe_m.', '.-g');
% pl(10) = plot(v_max_HFD, hoehe_m, '.-k');

pl(11) = plot(specs.Ma_CR * ISA.a(hoehe_CR), hoehe_CR, '*r'); % Design point

xlabel('v_{TAS} in m/s');
ylabel('H in m');

% 
% 





Ergebnisse_Flugleistung_1.hoehe_m = hoehe_m;
Ergebnisse_Flugleistung_1.c_A_j = c_A_j;
Ergebnisse_Flugleistung_1.c_A_F_j = c_A_F_j;
Ergebnisse_Flugleistung_1.c_W_inkomp_j = c_W_inkomp_j;
Ergebnisse_Flugleistung_1.j = j;
Ergebnisse_Flugleistung_1.v_EAS_j = v_EAS_j;
Ergebnisse_Flugleistung_1.v_EAS_j_CL = v_EAS_j_CL;
Ergebnisse_Flugleistung_1.v_EAS_j_DEC = v_EAS_j_DEC;
Ergebnisse_Flugleistung_1.Ma_j = Ma_j;
Ergebnisse_Flugleistung_1.Ma_j_CL = Ma_j_CL;
Ergebnisse_Flugleistung_1.Ma_J_DEC = Ma_j_DEC;
Ergebnisse_Flugleistung_1.S_S0_KF_j_CR = S_S0_KF_j_CR;
Ergebnisse_Flugleistung_1.S_S0_KF_j_CL = S_S0_KF_j_CL;
Ergebnisse_Flugleistung_1.S_S0_KF_j_DEC = S_S0_KF_j_DEC;
Ergebnisse_Flugleistung_1.S_S0_E = S_S0_E;
Ergebnisse_Flugleistung_1.S_S0_j_CR = S_S0_j;
Ergebnisse_Flugleistung_1.S_S0_j_CL = S_S0_j_CL;
Ergebnisse_Flugleistung_1.S_S0_j_DEC = S_S0_j_DEC;
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
Ergebnisse_Flugleistung_1.TAS_SET_H_DEC_vec = TAS_SEP_H_DEC_vec;
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


save Ergebnisse_FLugleistung_1.mat Ergebnisse_Flugleistung_1



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



