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

F_k = 

% Formel 1
alpha_f_max_VFFK = 0;


