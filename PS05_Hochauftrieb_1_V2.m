clc
clear
close

load Projekt_specs.mat
load Ergebnisse_Auftrieb_Momente.mat
load Ergebnisse_Fluegel_Tank_NP.mat
load Ergebnisse_Auftrieb_Momente.mat
load Ergebnisse_Leitwerke.mat
load Ergebnisse_Start_Landeanforderungen.mat
load Zwischenergebnisse_PS5_Fluegelflaechen.mat



%%%%%%%%%%%%%%%%%%%%%% Platzhalter-> Hier richtige Werte einfügen   von    PS4
CA_F_CR = 0.5;
CAalpha_F = Ergebnisse_Auftriebsverteilung.VWA.c_AF_anstieg;

%%%%%%%%%%%%%%%%%%%%%%

%% 1. Winkel Bezugsflügel zur Nullauftriebsrichtung
% CA_F_CR = Flügelauftriebsbeiwert Reiseflug aus Trimmbetrachtung
% CAalpha_F = Flügelauftriebsgradient im Reiseflug -> Diederich/weissinger

alpha_MAC_F_CR_0 = CA_F_CR / CAalpha_F;


%% 2. Nullanstellwinkel Profil aus Katalog 
% Aus Profilkatalog Alpha0 = -3° bei 0.83 Mach
alpha0profil = -2.6;

alpha_MAC_0_F = deg2rad(alpha0profil);
alpha_MAC_0_F_deg = rad2deg(alpha_MAC_0_F);

%% 3.Anstellwinkel Bezugsflügeltiefe

alpha_MAC_F_CR = alpha_MAC_F_CR_0 + alpha_MAC_0_F;
alpha_MAC_F_CR_deg = rad2deg(alpha_MAC_F_CR);


%% 4. Verwindungskorrektur zur Symmetriebene! ALSO MIT RUMPANTEIL?!
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
Delta_epsilon_root = VWA.epsilon_eta_Ru(1,93);
deltaEpsRoot_deg = rad2deg(Delta_epsilon_root);

psi_root = psi_sym_inst + Delta_epsilon_root;
psiRootDeg= rad2deg(psi_root);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Nullanstellwinkel Des Flugzeugs

F_H = HLW.F_aussen; % HLW Fläche % Richtig @Japser?

F = Ergebnisse_Fluegel.F; % Flügelfläche

CA = 1.5; % Auftrieb im Cruise 

CA_H = 0.3; % Aus Widerstand HLW CA


CA_F_cruise = CA - CA_H * 0.85 * (F_H/F);      % CA Flügel alleine im Cruise

% Braucht CA Flügel bei CAgesamt = 0
% Also Formel von oben mit CA= 0???????

CA_F = - CA_H * 0.85 * (F_H/F);

alpha_MAC_0 = alpha_MAC_0_F + (CA_F/CAalpha_F);

alpha_MAC_0_deg = rad2deg(alpha_MAC_0);


%Bezogen auf Flügelwurzelprofil

alpha_0_root = alpha_MAC_0 + Delta_epsilon_root + Delta_epsilon_sym;

alpha_0_root_deg = rad2deg(alpha_0_root);

% Alpha0  A/C bezüglich Rumpf

alpha_0 = alpha_0_root - psi_root;
alpha_0_deg = rad2deg(alpha_0);


%% Aufeglöste Polare ---->->--->->->-

streckung = Ergebnisse_Fluegel.streckung_phi25_max; % -> Richtig?

phi_50 = atan(tan(Ergebnisse_Fluegel.phi_25_max)-(4/streckung)* (0.5-0.25) * (1-Ergebnisse_Fluegel.lambda)/(1-Ergebnisse_Fluegel.lambda));
phi_50_deg = rad2deg(phi_50);

[~,a,~,~,~] = atmosisa(50);
M = (landeanvorderung.v_50 * 1.3)/a;

%CAalpha_F = Ergebnisse_Auftriebsverteilung.VWA.c_AF_anstieg;    %->Das
%wäre für Cruise

%Benötigt Werte:

% Alpha 0,F -> Oben bekannt aus Profilkatalog
alpha_0;
% CAAlpha bei MA 0.269 -> Berechnet
CA_alpha_lowspeed =  (pi * streckung)/(1 + sqrt(1 + ((streckung/2)^2) * (tan(phi_50)^2 + (1 - M^2)))); %->muss kleiner sein als bei Highspeed
% Alpha CA F Max -> Ablesen
delta_alpha_CA_F_max = deg2rad(3.6);
%CA F Max

CA22DMax = 1.4;

for i = 1: length(Ergebnisse_Fluegel.Fluegeltiefen_eta)
    CA_F_max_temp(i) = (CA22DMax * (Ergebnisse_Fluegel.Fluegeltiefen_eta(i)/Ergebnisse_Fluegel.l_m) - Ergebnisse_Auftriebsverteilung.gamma_b_eta(i)) / Ergebnisse_Auftriebsverteilung.gamma_a_eta(i);
end
[~,u] = min(CA_F_max_temp);

CA_F_max = Ergebnisse_Auftriebsverteilung.c_a_eta(u);


%%%% Alles aufsummieren für Alpha_CA_F_Max


alpha_CA_F_MAX = (CA_F_max/CA_alpha_lowspeed) + delta_alpha_CA_F_max + alpha_MAC_0_F;

alpha_CA_F_MAX_deg = rad2deg(alpha_CA_F_MAX);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotten

alphas = -6:0.01:alpha_CA_F_MAX_deg; % normal Plotten bis alphamax - delta alpha
CA_s = CA_alpha_lowspeed.*(deg2rad(alphas-alpha_0));
plot(alphas,CA_s)

% Noch den CA Abfall plotten !
hold on
title("Aufgelöste Flügelpolare ohne Hochauftriebshilfen")
ylabel("Auftriebsbeiwert CA")
xlabel("Anstellwinkel Alpha")
P1 = [0 0];
P2 = [-1 2];
plot(P1,P2,'black')
P3 = [-6 10];
P4 = [0 0];
plot(P3,P4,'black')

% Kritische Points
%Alpha MAX
plot(alpha_CA_F_MAX_deg, 0, 'xred')

%CA MAX
plot(0,CA_F_max, 'xred')

%Alpha 0
plot(alpha_0, 0,'xred')

plot([alpha_CA_F_MAX_deg , alpha_CA_F_MAX_deg-delta_alpha_CA_F_max],[CA_F_max CA_F_max],'xgreen')

grid on

save Ergebnisse_Hochauftrieb_1.mat psiRootDeg psi_sym_inst_deg alpha_CA_F_MAX_deg CA_F_max CA_alpha_lowspeed delta_alpha_CA_F_max alpha_MAC_0_F_deg CA_F alpha_MAC_0_F;



































