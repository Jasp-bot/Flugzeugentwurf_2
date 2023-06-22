load Projekt_specs.mat

specs.Ma_CR

% CA_F_CR = Flügelauftriebsbeiwert Reiseflug aus Trimmbetrachtung
% CAalpha_F = Flügelauftriebsgradient im Reiseflug -> Diederich/weissinger

alpha_MAC_F_CR_0 = CA_F_CR / CAalpha_F;

% Aus Profilkatalog Alpha0 = -3°
alpha0profil = -3;

alpha_MAC_0_F = deg2rad(alpha0profil);

%Anstellwinkel Bezugsflügeltiefe

alpha_MAC_F_CR = alpha_MAC_F_CR_0 + alpha_MAC_0_F;

% mittlere Verwindung % Aus diederich das Integral

Delta_epsilon_sym = 

% Einbauwinkel 
psi_sym_inst = alpha_MAC_F_CR + Delta_epsilon_sym;

%Einbauwinkel des Flügelwurzelprofils 



