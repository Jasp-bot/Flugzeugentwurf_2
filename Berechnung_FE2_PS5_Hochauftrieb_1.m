function Berechnung_FE2_PS5_Hochauftrieb_1

clc
clear all
close all

load Projekt_specs.mat
load Ergebnisse_Auftrieb_Momente.mat
load Ergebnisse_Fluegel_Tank_NP.mat
load Ergebnisse_Auftrieb_Momente.mat
load Ergebnisse_Leitwerke.mat
load Ergebnisse_Start_Landeanforderungen.mat
load Zwischenergebnisse_PS5_Fluegelflaechen.mat
load Ergebnisse_Widerstand_FE2.mat
load Ergebnisse_ISA_DATA.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Feste Variablen
%phi_50 = atan(tan(Ergebnisse_Fluegel.phi_25_max) - (4/Ergebnisse_Fluegel.streckung_phi25_max) * (0.5-0.25) * ((1 - Ergebnisse_Fluegel.lambda)/(1 + Ergebnisse_Fluegel.lambda))) ;

%CA_H = Ergebnisse_Widerstand_FE2.c_A_H(1);
%
CAalpha_F = (pi * Ergebnisse_Fluegel.streckung_phi25_max) / (1+sqrt(1 + ((Ergebnisse_Fluegel.streckung_phi25_max/2)^2) * (tan(Ergebnisse_Fluegel.phi_50)^2 + (1-specs.Ma_CR^2))));

%CAalpha_F = Ergebnisse_Auftriebsverteilung.VWA.c_AF_anstieg;    % Altes CA alpha aus FE1 über Diederich?
% Welcher Wert

% [~,a,~,~,~] = atmosisa(convlength(50,'ft','m'));
M = (landeanvorderung.v_50) / ISA.a(round(convlength(50,'ft','m')));
%M = 0.3;
% CA F CR UND CA_H bestimmen über Machzahl

CA_ind = round(Ergebnisse_Widerstand_FE2.stuetzstellen / 2  * specs.Ma_CR);

%CA_CR_ges = Ergebnisse_Widerstand_FE2.c_A_ges(CA_ind);
CA_CR_ges = Ergebnisse_Widerstand_FE2.c_A_ges_off_D(CA_ind);
%CA_CR_ges_offD = Ergebnisse_Widerstand_FE2.c_A_ges_off_D(CA_ind);
%%%

%CA_CR_H = Ergebnisse_Widerstand_FE2.c_A_H(CA_ind);
CA_CR_H = Ergebnisse_Widerstand_FE2.c_A_H_off_D(CA_ind);
%CA_CR_H_offD = Ergebnisse_Widerstand_FE2.c_A_H_off_D(CA_ind);


% Probe ob Berechnung richtig :) PASST!
CA_F_CR_2 = Ergebnisse_Widerstand_FE2.c_A_F_off_D(CA_ind);


CA_F_CR = CA_CR_ges - CA_CR_H * 0.85 * ((2 * HLW.F + 31) / Ergebnisse_Fluegel.F);
%CA_F_CR = 0.5;

F_H = (2 * HLW.F) + 5 + 26; % HLW Fläche % Richtig @Japser?  PRÜFEN CAD!

F = Ergebnisse_Fluegel.F; % Flügelfläche

CA = 0.0; % Auftrieb im Cruise des gesamten AC für alpha 0 berechnung


%% 1. Winkel Bezugsflügel zur Nullauftriebsrichtung
% CA_F_CR = Flügelauftriebsbeiwert Reiseflug aus Trimmbetrachtung
% CAalpha_F = Flügelauftriebsgradient im Reiseflug -> Diederich/weissinger

alpha_MAC_F_CR_0 = CA_F_CR / CAalpha_F;
alpha_MAC_F_CR_0_deg = rad2deg(alpha_MAC_F_CR_0);

%% 2. Nullanstellwinkel Profil aus Katalog
% Aus Profilkatalog Alpha0 = -2.8° bei 0.82 Mach
alpha0profil = -2.8;

alpha_MAC_0_F = deg2rad(alpha0profil);
alpha_MAC_0_F_deg = rad2deg(alpha_MAC_0_F);

%% 3.Anstellwinkel Bezugsflügeltiefe

alpha_MAC_F_CR = alpha_MAC_F_CR_0 + alpha_MAC_0_F;
alpha_MAC_F_CR_deg = rad2deg(alpha_MAC_F_CR);


%% 4. Verwindungskorrektur zur Symmetriebene! ALSO MIT RUMPfANTEIL?!
% mittlere Verwindung % Aus diederich das Integral
% Integral
%X = 0:.001:1;
%Delta_epsilon_sym = -trapz(VWA.epsilon_eta,X);
%Delta_epsilon_sym_deg  = rad2deg(Delta_epsilon_sym);

% Neue Methode
int = 0;
for x = 1: length(Ergebnisse_Auftriebsverteilung.eta)
    int = int + Ergebnisse_Auftriebsverteilung.gamma_a_eta(x) * VWA.epsilon_eta(x) * 0.001;
end

Delta_epsilon_sym = int * -1;
Delta_epsilon_sym_deg  = rad2deg(Delta_epsilon_sym);




%% 5. Einbauwinkel
psi_sym_inst = alpha_MAC_F_CR + Delta_epsilon_sym;
psi_sym_inst_deg = rad2deg(psi_sym_inst);


%% 6.Einbauwinkel des Flügelwurzelprofils

%eta_root = length(Ergebnisse_Fluegel.Fluegeltiefen_eta_Ru);
%delta_epsilon_root = dEps_dEta * eta_root
Delta_epsilon_root = VWA.epsilon_eta_Ru(1,end);
deltaEpsRoot_deg = rad2deg(Delta_epsilon_root);

psi_root = psi_sym_inst + Delta_epsilon_root;
psiRootDeg= rad2deg(psi_root);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% B. Nullanstellwinkel Des Flugzeugs


%CA_F_cruise = CA - CA_H * 0.85 * (F_H/F);      % CA Flügel alleine im Cruise

% Braucht CA Flügel bei CAgesamt = 0
% Also Formel von oben mit CA= 0???????

CA_F_0 = CA - CA_CR_H * 0.85 * (F_H/F);

alpha_MAC_0 = alpha_MAC_0_F + (CA_F_0/CAalpha_F);

alpha_MAC_0_deg = rad2deg(alpha_MAC_0);


%Bezogen auf Flügelwurzelprofil

alpha_0_root = alpha_MAC_0 + Delta_epsilon_root + Delta_epsilon_sym;

alpha_0_root_deg = rad2deg(alpha_0_root);

% Alpha0  A/C bezüglich Rumpf

alpha_0 = alpha_0_root - psi_root;
alpha_0_deg = rad2deg(alpha_0);
%alpha_0_deg = -2.3;


%% C. Aufeglöste Polare

%CAalpha_F = Ergebnisse_Auftriebsverteilung.VWA.c_AF_anstieg;    %->Das
%wäre für Cruise

%Benötigt Werte:

% Alpha 0,F -> Oben bekannt aus Profilkatalog
% CAAlpha bei MA 0.20 -> Berechnet
CA_alpha_lowspeed =  (pi * Ergebnisse_Fluegel.streckung_phi25_max)/(1 + sqrt(1 + ((Ergebnisse_Fluegel.streckung_phi25_max/2)^2) * (tan(Ergebnisse_Fluegel.phi_50)^2 + (1 - (M^2))))); %->muss kleiner sein als bei Highspeed

% Alpha CA F Max -> Ablesen
delta_alpha_CA_F_max = deg2rad(3.6);
delta_alpha_CA_F_max_deg = 3.6;
% CA F Max
CA22DMax = 1.4;

for i = 1: length(Ergebnisse_Fluegel.Fluegeltiefen_eta)
    CA_F_max_temp(i) = (CA22DMax * (Ergebnisse_Fluegel.Fluegeltiefen_eta(i)/Ergebnisse_Fluegel.l_m) - Ergebnisse_Auftriebsverteilung.gamma_b_eta(i)) / Ergebnisse_Auftriebsverteilung.gamma_a_eta(i);
end
[mac_hat_2h_lang_diesen_fehler_gesucht,~] = min(CA_F_max_temp);

CA_F_max = mac_hat_2h_lang_diesen_fehler_gesucht;


%%%% Alles aufsummieren für Alpha_CA_F_Max

alpha_CA_F_MAX = (CA_F_max/CA_alpha_lowspeed) + delta_alpha_CA_F_max + alpha_MAC_0;

alpha_CA_F_MAX_deg = rad2deg(alpha_CA_F_MAX);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotten

alphas = -6:0.1:alpha_CA_F_MAX_deg-delta_alpha_CA_F_max_deg; % normal Plotten bis alphamax - delta alpha
for i=1:length(alphas)
    CA_s(i) = CA_alpha_lowspeed * (deg2rad(alphas(i)-alpha_MAC_0_deg));%*0.8309;
end


%% Finde CA bei 6 ° für Hochauftrieb 2 @ Chatgpt
target = alpha_MAC_0_F_deg + 6;

% Find the index of the value closest to the target
[~, index] = min(abs(alphas - target));

% Get the closest value
closestValue = alphas(index);


CA_BEI_6_grad = CA_s(index);


% ALLE Wichtige Variabeln Speichern

HA1.alpha_Ca_F_max = alpha_CA_F_MAX_deg;

HA1.CA_F_max = CA_F_max;

HA1.CA_alpha_F = CA_alpha_lowspeed;

HA1.M_LS = M;

HA1.alpha0 = alpha_0_deg;

HA1.alpha_MAC_0 = alpha_MAC_0_deg ; % Nicht flügel

HA1.alpha_MAC_0_F = alpha_MAC_0_F_deg;

HA1.CA_alpha_F_CR = CAalpha_F;

HA1.alpha_MAC_F_CR_0 = alpha_MAC_F_CR_0_deg;

HA1.alpha_MAC_F_0 = alpha_MAC_0_F_deg;

HA1.alpha_MAC_F_CR = alpha_MAC_F_CR_deg;

HA1.delta_eps_sym = Delta_epsilon_sym_deg;

HA1.psi_sym_inst = psi_sym_inst_deg;

HA1.psi_root_inst = psiRootDeg;

HA1.delta_eps_root = deltaEpsRoot_deg;

HA1.CA_bei6_deg = CA_BEI_6_grad;

HA1.delta_alpha_CA_F_max = delta_alpha_CA_F_max_deg;

HA1.alphas = alphas;

HA1.CAs = CA_s;

HA1.legend = "alles in Degree -> Zuerst umwandeln";

save Ergebnisse_Hochauftrieb_1.mat HA1

%Altes Speichern
%save Ergebnisse_Hochauftrieb_1.mat psiRootDeg psi_sym_inst_deg alpha_CA_F_MAX_deg CA_F_max CA_alpha_lowspeed delta_alpha_CA_F_max alpha_MAC_0_F_deg CA_F alpha_MAC_0_F CA_alpha_lowspeed CA_s alphas alpha_0 alpha_MAC_0_deg delta_alpha_CA_F_max_deg CA_BEI_6_grad;

% plot(alphas,CA_s,'blue','LineWidth',1.5)
%
%
% %P = polyfit(alphas,CA_s,1)
%
% hold on
% title("Aufgelöste Flügelpolare ohne Hochauftriebshilfen","FontSize",15)
% ylabel("Auftriebsbeiwert des Flügels C_{A} in [-]","FontWeight","bold")
% xlabel("Anstellwinkel \alpha_{F} in °","FontWeight","bold")
% P1 = [0 0];
% P2 = [-1 2];
% plot(P1,P2,'black')
% P3 = [-10 25];
% P4 = [0 0];
% plot(P3,P4,'black')
%
% % Kritische Points
% %Alpha MAX
% plot(alpha_CA_F_MAX_deg, 0, 'xred')
%
% %CA MAX
% plot([0 20],[CA_F_max CA_F_max], 'red')
% plot(0, CA_F_max,"xred")
%
% %Alpha 0
% plot(alpha_MAC_0_deg, 0,'xred')
%
% plot([alpha_CA_F_MAX_deg  alpha_CA_F_MAX_deg-delta_alpha_CA_F_max_deg],[CA_F_max CA_F_max],'xgreen')
%
% grid on

%
% % Parameter der quadratischen Funktion
% a = 0.0058; % Koeffizient von x^2
% h = alpha_CA_F_MAX_deg; % x-Koordinate des Maximums
% k = CA_F_max; % y-Koordinate des Maximums
%
% % Bereich der x-Achse
% x = linspace(h-7, h+7, 10); % Hier können Sie den Bereich anpassen
%
% % Quadratische Funktion berechnen
% y = -a * (x - h).^2 + k;
%
% % Plot erstellen
% plot(x, y,"blue--",'LineWidth',1.5)
% plot(h, k, 'bx','LineWidth',1.5)