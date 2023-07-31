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

%%%%%%%%%%%%%%%%%%%%%% Platzhalter-> Hier richtige Werte einfügen von PS4
CA_F_CR = 0.5;
CAalpha_F = 6.36;
%%%%%%%%%%%%%%%%%%%%%%

% Variablen  %%%%%%%


%%%%%%%%%%%%%%%%%%%%

%% 1. Winkel Bezugsflügel zur Nullauftriebsrichtung
% CA_F_CR = Flügelauftriebsbeiwert Reiseflug aus Trimmbetrachtung
% CAalpha_F = Flügelauftriebsgradient im Reiseflug -> Diederich/weissinger

alpha_MAC_F_CR_0 = CA_F_CR / CAalpha_F;


%% 2. Nullanstellwinkel Profil aus Katalog 
% Aus Profilkatalog Alpha0 = -3°
alpha0profil = -3;

alpha_MAC_0_F = deg2rad(alpha0profil);


%% 3.Anstellwinkel Bezugsflügeltiefe

alpha_MAC_F_CR = alpha_MAC_F_CR_0 + alpha_MAC_0_F;
%test = rad2deg(alpha_MAC_F_CR)


%% 4. Verwindungskorrektur zur Symmetriebene! ALSO MIT RUMPANTEIL?!
% mittlere Verwindung % Aus diederich das Integral
% Integral 
Delta_epsilon_sym = -trapz(VWA.epsilon_eta); %?
%test = rad2deg(Delta_epsilon_sym)

%% 5. Einbauwinkel 
psi_sym_inst = alpha_MAC_F_CR + Delta_epsilon_sym;
%test = rad2deg(psi_sym_inst)


%% 6.Einbauwinkel des Flügelwurzelprofils 

eta_root = length(Ergebnisse_Fluegel.Fluegeltiefen_eta_Ru);
%delta_epsilon_root = dEps_dEta * eta_root 
Delta_epsilon_root = VWA.epsilon_eta_Ru(1,93);
deltaEpsRoot = rad2deg(Delta_epsilon_root)

psi_root = psi_sym_inst + Delta_epsilon_root;
psiRoot= rad2deg(psi_root)



%% Nullanstellwinkel Flugzeug
% Für zwei Szenarien -> CG max front / CG max rear

Fh = HLW.F_aussen;
F = Ergebnisse_Auftriebsverteilung.VWA.c_AF_anstieg;

CA_F = CA - CA_h * 0.85 * (Fh/F);


alpha_0_MAC = alpha_MAC_0_F + (CA_F/CAalpha_F);


alpha_0_root = alpha_0_MAC + Delta_epsilon_sym + Delta_epsilon_root;


alpha_0 = alpha_0_root - psi_root



%% Teil 2.
%Aufgelöste Polare über Weissinger Formel

CAalpha_F




