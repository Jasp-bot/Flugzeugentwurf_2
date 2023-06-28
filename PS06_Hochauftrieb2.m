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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Berehcnung F_K
% Tiefdecker
b_k_a = Ergebnisse_Fluegel.b * 0.7;
b_k_i = specs.D_rumpf + 0.03 * Ergebnisse_Fluegel.b; %  3% Abstand zum Rump über Halbspannweite

% Fläche Berechnen
F_k = 180;  % -trapz(VWA.epsilon_eta,X);
F = Ergebnisse_Fluegel.F;


% Alpha delta aus Abbildung 4 -> Unser Profil zwischen 16 und 18%
% -> Ablesen aus Bild: bei Klappentiefe 0,3(Vorgebeben+Video)#
%Berechnung mit 30° Ausschlag Fowler
delta_CA_MAX_K =  1.4;

%Formel 1

delta_CA_F_MAX_SF = delta_CA_MAX_K * (F_k/F) * cos(Ergebnisse_Fluegel.phi_25_max)^2;


% Ablesen aus Bild 5 -> Für 0,3 ck/c bei delta 20° Takeoff und
% Dickenverhältnis 17%
delta_CA_FK = 1.21;


%Formel 2
delta_CA_F_SF = delta_CA_FK * (F_k/F) * cos(Ergebnisse_Fluegel.phi_25_max)^2;


% Formel 3
%Korrektur für Verlängerung des Flügels

CA_alpha_F_SF = Ergebnisse_Auftriebsverteilung.VWA.c_AF_anstieg * (( (c_/c) - 1 ) * (F_k/F) + 1);

% Formel CAF_MAX VFFK aus Übungsfolie
% Aus Grafiken entnehmen
delta_CA_F_MAX_VF = 1;


CA_F_MAX_VFFK = CA_F_MAX + delta_CA_F_MAX_SF + delta_CA_F_MAX_VF;


%Letzter Einfluss - Übung
delta_alpha_F_MAX_VFFK = delta_alpha_CA_MAX; % ----> AUS EINBAUWINKEL/AUFGELÖSTE POLARE!!

%% ALLE RECHNUNGEN FÜR TAKEOFF UND LANDING KLAPPENSETTINGS UND CA´s


%% Momentenänderung - Trimmung noch möglich?





