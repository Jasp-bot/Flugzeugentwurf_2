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
load Ergebnisse_Hochauftrieb_1.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Berechnung F_K
% Tiefdecker
b_k_a = Ergebnisse_Fluegel.b * 0.7;
b_k_i = specs.D_rumpf + 0.03 * Ergebnisse_Fluegel.b; %  3% Abstand zum Rump über Halbspannweite

% Fläche Berechnen
F_k = 180;%<-Annahmne erstezen mit ->  % -trapz(VWA.epsilon_eta,X);
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
c_ = 0.3;
c = 1;
    %Eigentlich VK
CA_alpha_F_SF = Ergebnisse_Auftriebsverteilung.VWA.c_AF_anstieg * (( (c_/c) - 1 ) * (F_k/F) + 1);


%% Vorflügel
% Formel CAF_MAX VFFK aus Übungsfolie
% Aus Grafiken entnehmen
fsv = 1;
fdv = 1;
ffv = 1;
fphiv = 1;
delta_CA_MAX_0 = 1; %???????????

%Graphen 3
delta_CA_F_MAX_VF = delta_CA_MAX_0 * fsv * fdv * ffv * fphiv;

% Große Formel 2
CA_F_MAX_VFFK = CA_F_max + delta_CA_F_MAX_SF + delta_CA_F_MAX_VF;

%Einflüsse aus PS05
CA_F_max;
CA_F;
alpha_MAC_0_F;
delta_alpha_CA_F_max;



%Große Formel 1 -> Maximales alpha mit owler und Slats
alpha_F_MAX_VFFK = (CA_F_MAX_VFFK/CA_alpha_F_SF) + (CA_F*(alpha_MAC_0_F + deg2rad(6)) + delta_CA_F_SF)/(CA_alpha_F_SF) + deg2rad(6) + alpha_MAC_0_F + delta_alpha_CA_F_max;



%% ALLE RECHNUNGEN FÜR TAKEOFF UND LANDING KLAPPENSETTINGS





%% Momentenänderung - Trimmung noch möglich?


%% Widerstandszuwachs durch Klappen

%C_w_F = C_W_oK + delta_C_W_K

% delta_C_W_K = delta_C_W_P + delta_C_W_ind + delta_C_W_int +
% delta_C_W_Fahrwerk + delta_C_W_VF


% Profilwiderstandszuwachs

F_K_F_rumpf = 300;
F_K_F_ohnerumpf = 200;

F_K_F = F_K_F_rumpf - F_K_F_ohnerumpf;
                % Wird abglesen Landing und Takeoff einzeln rechnen!!!
delta_C_W_P = 0.065 * F_K_F * cos(Ergebnisse_Fluegel.phi_25_max);


% induzierter Widerstand -> w,v ablesen woher delta_...???
w = 0.3;
v = 0.5;
delta_C_A_K_phi0 = 1;

CA_F = CA * (1 - (delta_x_SP / l_mue)/(rh/l_mue)) - (c_m_0 + delta_C_M_K)/(rh/l_mue);

delta_C_W_ind = CA_F * delta_C_A_K_phi0 * v + delta_C_A_K_phi0^2 * w;


% inteferenz Widerstand
delta_C_W_int = (1/3) * delta_C_W_P;


% Fahrwerkswiderstand

delta_C_W_Fahrwerk = ((1.5 * F_vorder + 0.75 * F_hinter)/ F) * (1- 0.04 * ((CA_F + delta_CA_F_0 * (1.5 * (F/F_K)-1))/(l_HFW/l_mue)))^2;


% Vorflügelwiderstand

delta_C_W_VF = 


