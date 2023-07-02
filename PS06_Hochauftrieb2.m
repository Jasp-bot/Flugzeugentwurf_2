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
load Ergebnisse_CG.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Berechnung F_K
% Tiefdecker
b_k_a = Ergebnisse_Fluegel.b * 0.7;
b_k_i = specs.D_rumpf + 0.03 * Ergebnisse_Fluegel.b; %  3% Abstand zum Rump über Halbspannweite

% Fläche Berechnen
X = linspace(0,0.7,700);
Fluegel = Ergebnisse_Fluegel.Fluegeltiefen_eta(1,1:700);
F_k = 2 * (Ergebnisse_Fluegel.b/2)*trapz(X,Fluegel); % MIT RUMPF!

% Rumpf + extra stücke 3%
Fluegel2 = Ergebnisse_Fluegel.Fluegeltiefen_eta_oR(1,1:30);
X2 = linspace(0,0.03,30);
F_rumpf_rechteck = Ergebnisse_Fluegel.Fluegeltiefen_eta_Ru(1) * 6.21;
F_rumpf_Ecken = 2 * ((Ergebnisse_Fluegel.b)/2)*trapz(X2,Fluegel2);
F_rumpf_dreieck = Flaeche_im_Rumpf_oberes_dreieck;

F_rumpf = F_rumpf_dreieck + F_rumpf_Ecken + F_rumpf_rechteck;

% Somit F_k - F Rumpf = Klappenfläche

F_klappen = F_k - F_rumpf;


%% LANDING

% Kleine Formel CA_F_MAX -> AUS PS05
CA_F_max; %-> Übernehmen

% kleine Formel deltaCAFMAX,SF,phi
delta_Ca_max_SF_phi = 1.59;  % Ablesen aus Grafik bei Klappenausschlag ~45°
%                                                                       in Rad oder Grad?
delta_CA_F_max_SF_phi = delta_Ca_max_SF_phi * (F_klappen/Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max)^2;
    
%Annahme eta_opt für slat = 55°
% Vorflügel           % 0.93 = Faktor für Vorflügel/Slat        % Annahme
% 74% des Flügels mit Slats belegt      Annahme Slats tiefe = 1m -> 1/l_m=6,57 = 0.152 
delta_CA_F_max_VF = 0.93 * 0.64 * 0.69 * 1.0 * (cos(Ergebnisse_Fluegel.phi_25_max)^2);
                                                        %Rad oder Degree ? 


% Damit kann kleine Formel 1 berechnet werden!

CA_F_max_VFFK = CA_F_max + delta_CA_F_max_SF_phi + delta_CA_F_max_VF;

% Neuer Auftriebsanstieg
c = Ergebnisse_Fluegel.l_m; % Mittlere Flügeltiefe benutzen
c_ = 2.5 + Ergebnisse_Fluegel.l_m; % Quelle GPT und ACAMP
test2=c_/c;

CA_alpha_F_FK_phi = CA_alpha_lowspeed * (((c_/c)-1) * (F_klappen/Ergebnisse_Fluegel.F)+1);


%delta_CA_F_SF_phi
c_k = c - (Ergebnisse_Fluegel.l_m * 0.65);
test = c_k / c; % -> Für Landing ist dann der Faktor = 1.84
delta_CA_F_SF_phi = (F_k / Ergebnisse_Fluegel.F) * (cos(Ergebnisse_Fluegel.phi_25_max)^2) *  1.84;

% Große Formel VFFK
                                                    %In Übungfolien ist
                                                    %hier ein +
alpha_F_max_VFFK = (CA_F_max_VFFK / CA_alpha_F_FK_phi) - (((CA_F * (alpha_MAC_0_F + deg2rad(6)) + delta_CA_F_SF_phi)) / CA_alpha_F_FK_phi) + deg2rad(6) + delta_alpha_CA_F_max;

alpha_F_max_VFFK_deg = rad2deg(alpha_F_max_VFFK);

plot(alphas,CA_s,'blue')
hold on

title("Aufgelöste Flügelpolare mit Hochauftriebshilfen","FontSize",15)
ylabel("Auftriebsbeiwert des Flügels C_{A} in [-]","FontWeight","bold")
xlabel("Anstellwinkel \alpha_{F} in °","FontWeight","bold")


alphas = -14:0.01:alpha_F_max_VFFK_deg-delta_alpha_CA_F_max_deg; % normal Plotten bis alphamax - delta alpha
CA_s = CA_alpha_lowspeed.*(deg2rad(alphas-alpha_MAC_0_deg)) + delta_CA_F_max_SF_phi;
plot(alphas,CA_s,'green')


%% TAKEOFF
%%%%%%%%%%%%%%%%%%%% Gleiche Formeln andere Faktoren 


% Kleine Formel CA_F_MAX -> AUS PS05
CA_F_max; %-> Übernehmen

% kleine Formel deltaCAFMAX,SF,phi
delta_Ca_max_SF_phi_TO = 1.1;  % Ablesen aus Grafik bei Klappenausschlag ~45°
%                                                                       in Rad oder Grad?
delta_CA_F_max_SF_phi_TO = delta_Ca_max_SF_phi_TO * (F_klappen/Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max)^2;
    
%Annahme eta_opt für slat = 55°
% Vorflügel           % 0.93 = Faktor für Vorflügel/Slat        % Annahme
% 74% des Flügels mit Slats belegt      Annahme Slats tiefe = 1m -> 1/l_m=6,57 = 0.152 
delta_CA_F_max_VF_TO = 0.93 * 0.64 * 0.69 * 1.0 * (cos(Ergebnisse_Fluegel.phi_25_max)^2);
                                                        %Rad oder Degree ? 


% Damit kann kleine Formel 1 berechnet werden!

CA_F_max_VFFK_TO = CA_F_max + delta_CA_F_max_SF_phi_TO + delta_CA_F_max_VF_TO;

% Neuer Auftriebsanstieg
c_TO = Ergebnisse_Fluegel.l_m; % Mittlere Flügeltiefe benutzen
c__TO = 2+Ergebnisse_Fluegel.l_m; % Quelle GPT und ACAMP
test2=c_/c;

CA_alpha_F_FK_phi_TO = CA_alpha_lowspeed * (((c__TO/c_TO)-1) * (F_klappen/Ergebnisse_Fluegel.F)+1);


%delta_CA_F_SF_phi
c_k = c - (Ergebnisse_Fluegel.l_m * 0.65);
test = c_k / c; % -> Für Takeoff ist dann der Faktor = 1.3
delta_CA_F_SF_phi_TO = (F_k / Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max)^2 *  1.3;

% Große Formel VFFK
                                                        %IN ÜBUNG IST HIER
                                                        %+
alpha_F_max_VFFK_TO = (CA_F_max_VFFK_TO/CA_alpha_F_FK_phi_TO) - (((CA_F * (alpha_MAC_0_F + deg2rad(6)) + delta_CA_F_SF_phi_TO))/CA_alpha_F_FK_phi_TO) + deg2rad(6) + delta_alpha_CA_F_max;

alpha_F_max_VFFK_deg_TO = rad2deg(alpha_F_max_VFFK_TO);


%%%%%%%%%%%%%%%%%%%%%

alphas = -12:0.01:alpha_F_max_VFFK_deg_TO-delta_alpha_CA_F_max_deg; % normal Plotten bis alphamax - delta alpha
CA_s = CA_alpha_lowspeed.*(deg2rad(alphas-alpha_MAC_0_deg)) + delta_CA_F_max_SF_phi_TO;
plot(alphas,CA_s,'red')



grid on

% Kritische Points
%Alpha MAX
plot(alpha_CA_F_MAX_deg, 0, 'xred')

%CA MAX

%plot(0, CA_F_max,"xred")

%Alpha 0
%plot(alpha_MAC_0_F_deg, 0,'xblue')

plot([alpha_CA_F_MAX_deg  alpha_CA_F_MAX_deg-delta_alpha_CA_F_max_deg],[CA_F_max CA_F_max],'xgreen')

P1 = [0 0];
P2 = [-1 3];
plot(P1,P2,'black')
P3 = [-14 30];
P4 = [0 0];
plot(P3,P4,'black')

legend("0° Klappenausschlag","Landing mit 45° Klappenausschlag","Takeoff mit 20° Klappenausschlag",'','','','')



%% Momentenänderung - Trimmung noch möglich?
CM0 = -1;
deltaCM_HKK = 1;
r_h = 10;
l_mue = Ergebnisse_Fluegel.l_mue;
deltaXSP = 3;

% Formel 4
CA_MAX = ( CA_F_max + ( (CM0 + deltaCM_HKK)/( r_h/l_mue ) ) ) / (1 - ( deltaXSP/l_mue )/( r_h/l_mue ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Aus Formel.m
lambda       = 0.26;   %Zuspitzung des Flügels [-]
bf_s         = 0.7;    %prozentuale Spannweite der Hinterkantenklappen
%nach Skript Teil D Seite 114
dcMk_dcmK  = ( 3.739701 * lambda^4 - 10.762986 * lambda^3 + 12.791164 * lambda^2 - 8.860305 * lambda + 4.093244) * bf_s + ...
                            (-5.744622 * lambda^4 + 17.487598 * lambda^3 - 21.515763 * lambda^2 + 13.990786 * lambda - 4.219788) * bf_s^2 +...
                                ( 1.960553 * lambda^4 -  7.330970 * lambda^3 + 10.252347 * lambda^2 -  6.347223 * lambda + 1.466302) * bf_s^3;

%nach Skript Teil D Seite 114
dcMk_dcAK =  (-0.254693 * lambda^3 + 0.343337 * lambda^2 - 0.172780 * lambda + 0.333119) * bf_s +...
                            ( 0.986155 * lambda^3 - 1.620515 * lambda^2 + 1.049655 * lambda - 0.661691) * bf_s^2 + ...
                                (-0.735627 * lambda^3 + 1.282365 * lambda^2 - 0.877425 * lambda + 0.328347) * bf_s^3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Formel 7 + 8 

CA_REF_TO = CA_F_max_VFFK_TO/(1.2^2);
CA_REF_LDG = CA_F_max_VFFK / (1.3^2);

%Formel 10
ck_1 = c - c * 0.65;
theta = acos(2 * (ck_1/c) -1);

%Formel 9

dCM_dCA_deltaCA = 0.5 * (1 - (ck_1/c)) * (sin(theta)) / (pi - (theta - sin(theta)));  


%Formel 6
                                                            %delta Ca_FK
                                                            %bei Landing
                                                            %und Takeoff ->
                                                            %Abgelesen aus
                                                            %Grafik!!!!!!
delta_Cm_HKK_LDG = dCM_dCA_deltaCA * (c_/c) - ((CA_REF_LDG + 1.84 * (1-(F_k/Ergebnisse_Fluegel.F))) / (8)) * (c_/c) * ((c_/c)-1);
delta_Cm_HKK_TO = dCM_dCA_deltaCA * (c_/c) - ((CA_REF_TO + 1.3 * (1-(F_k/Ergebnisse_Fluegel.F))) / (8)) * (c_/c) * ((c_/c)-1);

%-> Für Landing ist dann der Faktor = 1.84
%-> Für takeoff ist dann der Faktor = 1.3


% Jetzt kann Formel 5 berechnet werden für Takeoff und Landing
                                                                                                                                                                                %SRichtig?
delta_CM_HKK_LDG = dcMk_dcmK * delta_Cm_HKK_LDG*(CA_REF_LDG) + 0.7 * ((Ergebnisse_Fluegel.streckung_phi25_max)/(1+2/Ergebnisse_Fluegel.streckung_phi25_max)) * dcMk_dcAK * delta_CA_F_SF_phi + tan(Ergebnisse_Fluegel.phi_25_max);

delta_CM_HKK_TO = dcMk_dcmK * delta_Cm_HKK_TO*(CA_REF_TO) + 0.7 * ((Ergebnisse_Fluegel.streckung_phi25_max)/(1+2/Ergebnisse_Fluegel.streckung_phi25_max)) * dcMk_dcAK * delta_CA_F_SF_phi_TO + tan(Ergebnisse_Fluegel.phi_25_max);


% WAS JETZT mit diesen Werten ?!?!?





%% Widerstandszuwachs durch Klappen

C_w_F = C_W_oK + delta_C_W_K

delta_C_W_K = delta_C_W_P + delta_C_W_ind + delta_C_W_int;
delta_C_W_Fahrwerk + delta_C_W_VF


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

delta_C_W_VF = 0;













%% OLD/ ARCHIV
% 
% 
% % Alpha delta aus Abbildung 4 -> Unser Profil zwischen 16 und 18%
% % % -> Ablesen aus Bild: bei Klappentiefe 0,3(Vorgebeben+Video)#
% % %Berechnung mit 30° Ausschlag Fowler
% % delta_CA_MAX_K =  1.4;
% % 
% % %Formel 1
% % delta_CA_F_MAX_SF = delta_CA_MAX_K * (F_k/F) * cos(Ergebnisse_Fluegel.phi_25_max)^2;
% % 
% % 
% % % Ablesen aus Bild 5 -> Für 0,3 ck/c bei delta 20° Takeoff und
% % % Dickenverhältnis 17%
% % delta_CA_FK = 1.21;
% % 
% % 
% % %Formel 2
% % delta_CA_F_SF = delta_CA_FK * (F_k/F) * cos(Ergebnisse_Fluegel.phi_25_max)^2;
% % 
% % 
% % % Formel 3
% % %Korrektur für Verlängerung des Flügels
% % c_ = 0.3;
% % c = 1;
% %     %Eigentlich VK
% % CA_alpha_F_SF = Ergebnisse_Auftriebsverteilung.VWA.c_AF_anstieg * (( (c_/c) - 1 ) * (F_k/F) + 1);
% % 
% % 
% % %% Vorflügel
% % % Formel CAF_MAX VFFK aus Übungsfolie
% % % Aus Grafiken entnehmen
% % fsv = 1;
% % fdv = 1;
% % ffv = 1;
% % fphiv = 1;
% % delta_CA_MAX_0 = 1; %???????????
% % 
% % %Graphen 3
% % delta_CA_F_MAX_VF = delta_CA_MAX_0 * fsv * fdv * ffv * fphiv;
% % 
% % % Große Formel 2
% % CA_F_MAX_VFFK = CA_F_max + delta_CA_F_MAX_SF + delta_CA_F_MAX_VF;
% % 
% % %Einflüsse aus PS05
% % CA_F_max;
% % CA_F;
% % alpha_MAC_0_F;
% % delta_alpha_CA_F_max;
% % 
% % 
% % 
% % %Große Formel 1 -> Maximales alpha mit owler und Slats
% % alpha_F_MAX_VFFK = (CA_F_MAX_VFFK/CA_alpha_F_SF) + (CA_F*(alpha_MAC_0_F + deg2rad(6)) + delta_CA_F_SF)/(CA_alpha_F_SF) + deg2rad(6) + alpha_MAC_0_F + delta_alpha_CA_F_max;
% % 
% % 
% % 
% % %% ALLE RECHNUNGEN FÜR TAKEOFF UND LANDING KLAPPENSETTINGS
% % 
% % 
% % 
