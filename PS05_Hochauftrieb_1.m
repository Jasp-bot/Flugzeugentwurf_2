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

% Variablen  %%%%%%%


%%%%%%%%%%%%%%%%%%%%

%% 1. Winkel Bezugsflügel zur Nullauftriebsrichtung
% CA_F_CR = Flügelauftriebsbeiwert Reiseflug aus Trimmbetrachtung
% CAalpha_F = Flügelauftriebsgradient im Reiseflug -> Diederich/weissinger

alpha_MAC_F_CR_0 = CA_F_CR / CAalpha_F;


%% 2. Nullanstellwinkel Profil aus Katalog 
% Aus Profilkatalog Alpha0 = -3° bei 0.83 Mach
alpha0profil = -3;

alpha_MAC_0_F = deg2rad(alpha0profil);


%% 3.Anstellwinkel Bezugsflügeltiefe

alpha_MAC_F_CR = alpha_MAC_F_CR_0 + alpha_MAC_0_F;
%test = rad2deg(alpha_MAC_F_CR)


%% 4. Verwindungskorrektur zur Symmetriebene! ALSO MIT RUMPANTEIL?!
% mittlere Verwindung % Aus diederich das Integral
% Integral 
X = 0:.001:1;
Delta_epsilon_sym = -trapz(VWA.epsilon_eta,X); %% Fehler????
%test = rad2deg(Delta_epsilon_sym)

%% 5. Einbauwinkel 
psi_sym_inst = alpha_MAC_F_CR + Delta_epsilon_sym;
%test = rad2deg(psi_sym_inst)


%% 6.Einbauwinkel des Flügelwurzelprofils 

eta_root = length(Ergebnisse_Fluegel.Fluegeltiefen_eta_Ru);
%delta_epsilon_root = dEps_dEta * eta_root 
Delta_epsilon_root = VWA.epsilon_eta_Ru(1,93);
%deltaEpsRoot = rad2deg(Delta_epsilon_root)

psi_root = psi_sym_inst + Delta_epsilon_root;
psiRootDeg= rad2deg(psi_root);



%% Nullanstellwinkel Flugzeug
% Für zwei Szenarien -> CG max front / CG max rear

Fh = HLW.F_aussen ; % -> Ist das der richtige Wert?
F = Ergebnisse_Fluegel.F; % -> Richtig ? Ganze Flügelfläche ohne Rumpf soll das sein1

CA_front = 0.1;   % Hier noch richtige Werte -> Muss bei CAGesamt = 0!!!!!!!
CA_aft = -0.1;   %
CA_h = 0.5;     %

%% Für CG Front
CA_F_front = CA_front - CA_h * 0.85 * (Fh/F);

alpha_0_MAC_front = alpha_MAC_0_F + (CA_F_front/CAalpha_F);


alpha_0_root_front = alpha_0_MAC_front + Delta_epsilon_sym + Delta_epsilon_root;


alpha_0_front = alpha_0_root_front - psi_root;
alpha_0_aft = rad2deg(alpha_0_front)

%% Für CG AFT
CA_F_aft = CA_aft - (CA_h * 0.85 * (Fh/F));

alpha_0_MAC_aft = alpha_MAC_0_F + (CA_F_aft/CAalpha_F);

alpha_0_root_aft = alpha_0_MAC_aft + Delta_epsilon_sym + Delta_epsilon_root;


alpha_0_aft = alpha_0_root_aft - psi_root;
alpha_o_aft_deg = rad2deg(alpha_0_aft)

%% Teil 2.
%Aufgelöste Polare über Weissinger Formel -> Schritte wie in Übungsfolien

%% 1. CA Alpha

streckung = Ergebnisse_Fluegel.streckung_phi25_max; % -> Richtig?
phi_50 = atan(tan(Ergebnisse_Fluegel.phi_25_max)-(4/streckung)* (0.5-0.25) * (1-Ergebnisse_Fluegel.lambda)/(1-Ergebnisse_Fluegel.lambda));
phi_50_deg = rad2deg(phi_50);
[~,a,~,~,~] = atmosisa(50);
M = (landeanvorderung.v_50 * 1.3)/a;

%Fällt weg weil schon berechnet?
CA_alpha =  (pi * streckung)/(1 + sqrt(1 + ((streckung/2)^2) * (tan(phi_50)^2 + (1 - M^2))));
%CA_alpha = Ergebnisse_Auftriebsverteilung.VWA.c_AF_anstieg;    -> Richtig?

%% 2. delta CAf_max

% Wird aus Plot abgelesen für phi_VK = 32.2
rad2deg(Ergebnisse_Fluegel.phi_VK_max)
delta_alpha_CA_F_max = deg2rad(3.6); 

%% 3. CA_F_max
% Maximaler Auftriebsbeiwert am Flügel
% = min-Wert der folgenden Formel
Ca_2D = 1.4;

for i = 1: length(Ergebnisse_Fluegel.Fluegeltiefen_eta)
    CA_F_max(i) = (Ergebnisse_Auftriebsverteilung.c_a_eta(i) * (Ergebnisse_Fluegel.Fluegeltiefen_eta(i)) - Ergebnisse_Auftriebsverteilung.gamma_b_eta(i)) / Ergebnisse_Auftriebsverteilung.gamma_a_eta(i);
end
[~,u] = min(CA_F_max);

CA_F_max = Ergebnisse_Auftriebsverteilung.c_a_eta(u)

%% 4. alphaCaf_max

alpha_Ca_f_max = (CA_F_max/ CA_alpha) * alpha_0_aft - delta_alpha_CA_F_max
deg2rad(alpha_Ca_f_max)


%% PLotting Aufgelöste Polare

alphas = -6:0.1:(alpha_Ca_f_max-delta_alpha_CA_F_max); % normal Plotten bis alphamax - delta alpha

plot(alphas,(CA_alpha.*alphas))
title("Aufgelöste Flügelpolare ohne Hochauftriebshilfen")
ylabel("Auftriebsbeiwert CA")
xlabel("Anstellwinkel Alpha")
grid on