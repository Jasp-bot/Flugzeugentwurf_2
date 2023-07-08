clc
clear
close

%% Aus FE 2

load Projekt_specs.mat
load Ergebnisse_Auftrieb_Momente.mat
load Ergebnisse_Fluegel_Tank_NP.mat
load Ergebnisse_Auftrieb_Momente.mat
load Ergebnisse_Leitwerke.mat
load Ergebnisse_Start_Landeanforderungen.mat
load Zwischenergebnisse_PS5_Fluegelflaechen.mat
load Ergebnisse_Hochauftrieb_1.mat
load Ergebnisse_CG.mat

%% Aus FE 1

load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat
load Ergebnisse_Basis_stat_m.mat
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat
load Ergebnisse_Start_Landeanforderungen.mat
load Ergebnisse_Fluegel_Tank_NP.mat;

load Ergebnisse_Endwerte_Iteration_V1.mat
Endwerte_Iteration = Berechnungen_PS10_Widerstand;

%% Laden der Dateien

load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Basis_stat_m.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Start_Landeanforderungen.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Widerstand.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Berechnung F_K
% Tiefdecker
b_k_a = Ergebnisse_Fluegel.b * 0.7;
b_k_i = specs.D_rumpf + 0.03 * Ergebnisse_Fluegel.b; %  3% Abstand zum Rump über Halbspannweite

% Fläche Berechnen
X = linspace(0,0.7,600);
Fluegel = Ergebnisse_Fluegel.Fluegeltiefen_eta(1,1:600);
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
delta_Ca_max_SF_phi = 1.57;  % Ablesen aus Grafik bei Klappenausschlag ~45°
%                                                                       in Rad oder Grad?
delta_CA_F_max_SF_phi = delta_Ca_max_SF_phi * (F_klappen/Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max^2);
    
%Annahme eta_opt für slat = 55°
% Vorflügel           % 0.93 = Faktor für Vorflügel/Slat        % Annahme
% 74% des Flügels mit Slats belegt      Annahme Slats tiefe = 1m -> 1/l_m=6,57 = 0.152 
delta_CA_F_max_VF = 0.93 * 0.94 * 0.99 * 0.8 * (cos(rad2deg(Ergebnisse_Fluegel.phi_25_max)-5))^2;
                                                        %Rad oder Degree ? 


% Damit kann kleine Formel 1 berechnet werden!

CA_F_max_VFFK = CA_F_max + delta_CA_F_max_SF_phi + delta_CA_F_max_VF;

% Neuer Auftriebsanstieg
c = Ergebnisse_Fluegel.l_m; % Mittlere Flügeltiefe benutzen
c_ = 2.5 + Ergebnisse_Fluegel.l_m; % Quelle GPT und ACAMP
%test2=c_/c;

CA_alpha_F_FK_phi = CA_alpha_lowspeed * (((c_/c)-1) * (F_klappen/Ergebnisse_Fluegel.F)+1);


%delta_CA_F_SF_phi
c_k = c - (Ergebnisse_Fluegel.l_m * 0.65);
test = c_k / c; % -> Für Landing ist dann der Faktor = 1.84
delta_C_a_FK = 1.84;
delta_CA_F_SF_phi = (F_klappen / Ergebnisse_Fluegel.F) * (cos(Ergebnisse_Fluegel.phi_25_max^2)) *  delta_C_a_FK;

% Große Formel VFFK
                                                    %In Übungfolien ist
                                                    %hier ein +                                                                             %+
alpha_F_max_VFFK = (CA_F_max_VFFK / CA_alpha_F_FK_phi) - (((CA_F * (alpha_MAC_0_F + deg2rad(6)) + delta_CA_F_SF_phi)) / CA_alpha_F_FK_phi) + deg2rad(6) + delta_alpha_CA_F_max + alpha_MAC_0_F;

alpha_F_max_VFFK_deg = rad2deg(alpha_F_max_VFFK);

plot(alphas,CA_s,'blue')
hold on

title("Aufgelöste Flügelpolare mit Hochauftriebshilfen","FontSize",15)
ylabel("Auftriebsbeiwert des Flügels C_{A} in [-]","FontWeight","bold")
xlabel("Anstellwinkel \alpha_{F} in °","FontWeight","bold")


alphas = -15:0.01:alpha_F_max_VFFK_deg-delta_alpha_CA_F_max_deg; % normal Plotten bis alphamax - delta alpha
CA_s = CA_alpha_lowspeed.*(deg2rad(alphas-alpha_MAC_0_deg)) + delta_CA_F_max_SF_phi;
plot(alphas,CA_s,'green')


%% TAKEOFF
%%%%%%%%%%%%%%%%%%%% Gleiche Formeln andere Faktoren 


% Kleine Formel CA_F_MAX -> AUS PS05
CA_F_max; %-> Übernehmen

% kleine Formel deltaCAFMAX,SF,phi
delta_Ca_max_SF_phi_TO = 1.1;  % Ablesen aus Grafik bei Klappenausschlag ~45°
%                                                                       in Rad oder Grad?
delta_CA_F_max_SF_phi_TO = delta_Ca_max_SF_phi_TO * (F_klappen/Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max^2);
    
%Annahme eta_opt für slat = 55°
% Vorflügel           % 0.93 = Faktor für Vorflügel/Slat        % Annahme
% 74% des Flügels mit Slats belegt      Annahme Slats tiefe = 1m -> 1/l_m=6,57 = 0.152 
delta_CA_F_max_VF_TO = 0.93 * 0.94 * 0.99 * 0.8 * (cos(rad2deg(Ergebnisse_Fluegel.phi_25_max)-5))^2; %0.93 * 0.64 * 0.69 * 1.0 * (cos(Ergebnisse_Fluegel.phi_25_max)^2);
                                                        %Rad oder Degree ? 


% Damit kann kleine Formel 1 berechnet werden!

CA_F_max_VFFK_TO = CA_F_max + delta_CA_F_max_SF_phi_TO + delta_CA_F_max_VF_TO;

% Neuer Auftriebsanstieg
c_TO = 0.5; %Ergebnisse_Fluegel.l_m; % Mittlere Flügeltiefe benutzen
c__TO = .686; %1.5+Ergebnisse_Fluegel.l_m; % Quelle GPT und ACAMP
test2= c__TO/c_TO;

CA_alpha_F_FK_phi_TO = CA_alpha_lowspeed * (((c__TO/c_TO)-1) * (F_klappen/Ergebnisse_Fluegel.F) + 1);


%delta_CA_F_SF_phi
c_k = c - (Ergebnisse_Fluegel.l_m * 0.65);
test = c_k / c; % -> Für Takeoff ist dann der Faktor = 1.3
delta_C_a_FK_TO = 1.3;
delta_CA_F_SF_phi_TO = (F_klappen / Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max^2) *  delta_C_a_FK_TO;

% Große Formel VFFK
                                                        %IN ÜBUNG IST HIER
                                                        %+                                                                                           %+ 
alpha_F_max_VFFK_TO = (CA_F_max_VFFK_TO/CA_alpha_F_FK_phi_TO) - (((CA_F * (alpha_MAC_0_F + deg2rad(6)) + delta_CA_F_SF_phi_TO))/CA_alpha_F_FK_phi_TO) + deg2rad(6) + delta_alpha_CA_F_max + alpha_MAC_0_F;

alpha_F_max_VFFK_deg_TO = rad2deg(alpha_F_max_VFFK_TO);


%%%%%%%%%%%%%%%%%%%%%

alphas = -12:0.01:alpha_F_max_VFFK_deg_TO-delta_alpha_CA_F_max_deg; % normal Plotten bis alphamax - delta alpha
CA_s = CA_alpha_lowspeed.*(deg2rad(alphas-alpha_MAC_0_deg)) + delta_CA_F_max_SF_phi_TO;
plot(alphas,CA_s,'red')
grid on
plot([alpha_F_max_VFFK_deg_TO-delta_alpha_CA_F_max_deg],[0],'redx')
plot([alpha_F_max_VFFK_deg-delta_alpha_CA_F_max_deg],[0],'greenx')

plot([0],[CA_F_max_VFFK_TO],'redx')
plot([0],[CA_F_max_VFFK],'greenx')


plot(alpha_CA_F_MAX_deg, 0, 'xblue')
%CA MAX
plot(0, CA_F_max,"xblue")

%% Abfall der Polaren Plotten

% Parameter der quadratischen Funktion
a = 0.0057; % Koeffizient von x^2
h = alpha_F_max_VFFK_deg; % x-Koordinate des Maximums
k = CA_F_max_VFFK; % y-Koordinate des Maximums

% Bereich der x-Achse
x = linspace(h-6, h+7, 10); % Hier können Sie den Bereich anpassen

% Quadratische Funktion berechnen
y = -a * (x - h).^2 + k;

% Plot erstellen
plot(x, y,"green--")
plot(h, k, 'gx')


% Parameter der quadratischen Funktion
a = 0.0052; % Koeffizient von x^2
h = alpha_F_max_VFFK_deg_TO; % x-Koordinate des Maximums
k = CA_F_max_VFFK_TO; % y-Koordinate des Maximums

% Bereich der x-Achse
x = linspace(h-7, h+7, 10); % Hier können Sie den Bereich anpassen

% Quadratische Funktion berechnen
y = -a * (x - h).^2 + k;

% Plot erstellen
plot(x, y,"red--")
plot(h, k, 'rx')



% Parameter der quadratischen Funktion
a = 0.0058; % Koeffizient von x^2
h = alpha_CA_F_MAX_deg; % x-Koordinate des Maximums
k = CA_F_max; % y-Koordinate des Maximums

% Bereich der x-Achse
x = linspace(h-7, h+7, 10); % Hier können Sie den Bereich anpassen

% Quadratische Funktion berechnen
y = -a * (x - h).^2 + k;

% Plot erstellen
plot(x, y,"blue--")
plot(h, k, 'bx')

% Kritische Points
%Alpha MAX
%plot(alpha_CA_F_MAX_deg, 0, 'xred')

%CA MAX

%plot(0, CA_F_max,"xred")

%Alpha 0
%plot(alpha_MAC_0_F_deg, 0,'xblue')

%plot([alpha_CA_F_MAX_deg  alpha_CA_F_MAX_deg-delta_alpha_CA_F_max_deg],[CA_F_max CA_F_max],'xgreen')

P1 = [0 0];
P2 = [-1 3];
plot(P1,P2,'black')
P3 = [-14 30];
P4 = [0 0];
plot(P3,P4,'black')

ylim([-1, 3])

legend("Clean Konfiguration mit 0° Klappenausschlag","Landing mit 45° Klappenausschlag","Takeoff mit 20° Klappenausschlag",'','','','','Location', 'southeast')



%% Momentenänderung - Trimmung noch möglich?
CM0 = -1;
deltaCM_HKK = 1;
r_h = 10;
l_mue = Ergebnisse_Fluegel.l_mue;
deltaXSP = 3;

% Formel 4      -> Muss hier für TO und LDG berechnet werden und bei 7u.8
% eingesetzt werden
CA_MAX = ( CA_F_max + ( (CM0 + deltaCM_HKK)/( r_h/l_mue ) ) ) / (1 - ( deltaXSP/l_mue )/( r_h/l_mue ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Aus Formel.m
lambda       = 1;   %Zuspitzung des Flügels [-]      %% JASPER HILFE WELCHE WERTE RICHTIG?
bf_s         = 0.7;    %prozentuale Spannweite der Hinterkantenklappen      %% Iterabel?
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

CA_REF_TO = CA_F_max_VFFK_TO/(1.2^2);       % Etwas zu klein
CA_REF_LDG = CA_F_max_VFFK / (1.3^2);

%Formel 10
ck_1 = c - c * 0.68;            % Nach 65% Holm etwas Toleranz
theta = acos(2 * (ck_1/c) -1);

%Formel 9

dCM_dCA_deltaCA = 0.5 * (1 - (ck_1/c)) * (sin(theta)) / (pi - (theta - sin(theta)));  


%Formel 6
                                                            %delta Ca_FK
                                                            %bei Landing
                                                            %und Takeoff ->
                                                            %Abgelesen aus
                                                            %Grafik!!!!!!
delta_Cm_HKK_LDG = dCM_dCA_deltaCA * (c_/c) - ((CA_REF_LDG + 1.84 * (1-(F_klappen/Ergebnisse_Fluegel.F))) / (8)) * (c_/c) * ((c_/c)-1);
delta_Cm_HKK_TO = dCM_dCA_deltaCA * (c_/c) - ((CA_REF_TO + 1.3 * (1-(F_klappen/Ergebnisse_Fluegel.F))) / (8)) * (c_/c) * ((c_/c)-1);

%-> Für Landing ist dann der Faktor = 1.84
%-> Für takeoff ist dann der Faktor = 1.3


% Jetzt kann Formel 5 berechnet werden für Takeoff und Landing
                                                                                                                                                                                %SRichtig?
delta_CM_HKK_LDG = dcMk_dcmK * delta_Cm_HKK_LDG*(CA_REF_LDG) + 0.7 * ((Ergebnisse_Fluegel.streckung_phi25_max)/(1+2/Ergebnisse_Fluegel.streckung_phi25_max)) * dcMk_dcAK * delta_CA_F_SF_phi + tan(Ergebnisse_Fluegel.phi_25_max);

delta_CM_HKK_TO = dcMk_dcmK * delta_Cm_HKK_TO*(CA_REF_TO) + 0.7 * ((Ergebnisse_Fluegel.streckung_phi25_max)/(1+2/Ergebnisse_Fluegel.streckung_phi25_max)) * dcMk_dcAK * delta_CA_F_SF_phi_TO + tan(Ergebnisse_Fluegel.phi_25_max);


% WAS JETZT mit diesen Werten ?!?!?





%% Widerstandszuwachs durch Klappen


% Profilwiderstand % Unterschiedlich bei TO /LDG
delta_C_W_P_phi = 0.081;
delta_C_W_P_phi_TO = 0.015;

delta_C_W_P = delta_C_W_P_phi * (F_klappen/Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max);
delta_C_W_P_TO = delta_C_W_P_phi_TO * (F_klappen/Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max);

%Induzierter Widerstand % Unterschiedlich für TO/LDG

delta_Ca_max_SF_phi_TO;
delta_Ca_max_SF_phi;

delta_Ca_max_SF_phi;


% Klappenspannweite durch Spannweite 
temp = 0.6;
Ergebnisse_Fluegel.lambda;
v = 0.001;                    % Richtig abgelesen?!?!
w = 0.0073;

delta_CM_K = 0.5; % Wie kommt man auf den Wert ?
CA = 0.5;
CA_F = CA * (1- ((deltaXSP/Ergebnisse_Fluegel.l_mue)/(r_h/Ergebnisse_Fluegel.l_mue))) - ((CM0+delta_CM_K)/(r_h/Ergebnisse_Fluegel.l_mue))/(r_h/Ergebnisse_Fluegel.l_mue);


delta_C_W_Ind = CA_F * delta_C_a_FK * v + delta_Ca_max_SF_phi^2 *w;

delta_C_W_Ind_TO = CA_F * delta_C_a_FK * v + delta_Ca_max_SF_phi_TO^2 *w;


% Inteferenzwiderstand % Unterschiedlich bei TO/LDG


delta_C_W_Inf = (1/3) * delta_C_W_P;
delta_C_W_Inf_TO = (1/3) * delta_C_W_P_TO;

% Vorflügelwiderstand % Gleich für TO/LDG weil gleiche Settings gewählt

C_W_P_Min_RE = 0.3; %% JASPER WERT PLS!

%Fläche Vorflügel ist von Holm vorne bis Vorderkante    % Über mittlere
%Flügeltiefe ok?
tiefe_Slats = Ergebnisse_Fluegel.l_m * 0.3;
laenge_Slats = 24; % m
laenge_Fluegel = Ergebnisse_Fluegel.b/2 - 3.15; %Ohne Rumpf?

% A*b ansatz ? Genau genug
F_VF = tiefe_Slats * laenge_Slats;


delta_CW_VF = C_W_P_Min_RE * (F_VF/Ergebnisse_Fluegel.F) * (laenge_Slats/laenge_Fluegel) *cos(Ergebnisse_Fluegel.phi_25_max);


% Fahrwerkswiderstand % Bleibt gleich bei TO/LDG

% Braucht finale Werte aus Fahrwerk 
durchmesser = 2;
breite = 0.5;
F_vorder = durchmesser*breite*4;   % Gleiche Werte weil vorne und hinten gleich groß sind
F_hinter = durchmesser*breite*4;
l_HFW = 30; % Ríchtiger Wert aus CG?
delta_CA_F_0 = 0.3; % Richtiger Wert aus wo?

delta_C_W_Fahrwerk = ((1.5 * F_vorder + 0.75 * F_hinter)/Ergebnisse_Fluegel.F) * (1 - 0.04 * ((CA_F + delta_CA_F_0 *(1.5 * (Ergebnisse_Fluegel.F / F_klappen)-1))/(l_HFW/Ergebnisse_Fluegel.l_m)));


%% Gesamtwiderstand durch klappen

delta_CW_klappen_LDG = delta_CW_VF + delta_C_W_Ind + delta_C_W_Inf + delta_C_W_P;
delta_CW_klappen_TO = delta_CW_VF + delta_C_W_Ind_TO + delta_C_W_Inf_TO + delta_C_W_P_TO;

delta_CW_klappen_LDG_fahrwerk = delta_CW_VF + delta_C_W_Ind + delta_C_W_Inf + delta_C_W_P + delta_C_W_Fahrwerk;
delta_CW_klappen_TO_fahrwerk = delta_CW_VF + delta_C_W_Ind_TO + delta_C_W_Inf_TO + delta_C_W_P_TO + delta_C_W_Fahrwerk;

%% Aufaddieren zu dem vorhandenen Widerstand

Widerstand.C_w_clean_all = Widerstand.C_w_clean_all + delta_CW_VF + delta_C_W_Ind + delta_C_W_Inf + delta_C_W_P;



%% Prüfen ob CA max passt für alle Flugphasen

% check ob StartCA erreicht
startschub.c_A_max_thrust_match < CA_F_max_VFFK_TO

% Check ob lande CA erreicht
landeanvorderung.c_A_max_LDG < CA_F_max_VFFK










%% Plotten der Reziproken Gleitzahlen mit und ohne Fahrwerk für TO,LDG und Clean
% Gleitzahl -> x ist Cw y ist Ca

% Alle CA´s bestimmen für verschiedene Alpha
% Alle CW_s bestimmen für verschiedene Alpha

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flughoehe = specs.flight_level * 10^2 ;                         % in ft
% hoehe = round(unitsratio('m','ft')*Flughoehe);     % in m
% 
% p_CR = ISA.p(hoehe); % p_Basis Druckvector aus ISA+15
% rho_CR = ISA.rho(hoehe);  % rho_Basis Dichtevector aus ISA+15
% 
% G_To = Ergebnisse_stat_Flaechenbelastung.G_To;
% 
% %% Berechnung Reynolds
% 
% Ind_W.u_cr = specs.Ma_CR * ISA.a(hoehe) ; % fuer reynolds Umstroemungsgeachw
% 
% 
% v_kin = ISA.kin_visk(hoehe); % kin vis
% Ind_W.Re_crit = (Ind_W.u_cr * specs.l_rumpf)/v_kin;
% % cf_turb = 0.455 / (log(Re_crit)^2.58); % const
% 
% %% Flügel
%                                     %d/l Aus Wingdata
% Ind_W.CwF_Fluegel_I = 0.0054 * 1 * (1 + 3 * 0.13 * (cos(NP.phi_25_I))^2) * DT.F_I;
% Ind_W.CwF_Fluegel_A = 0.0054 * 1 * (1 + 3 * 0.13 * (cos(NP.phi_25_A))^2) * DT.F_A;
% %CwF_Fluegel_R = 0.0054 * 1 * (1 + 3 * 0.13 * (cos(phi_25_R))^2) * F_R; %
% %Hat keinen Einfluss auf den Cw wert
% 
% Ind_W.CwF_Fluegel = (Ind_W.CwF_Fluegel_I + Ind_W.CwF_Fluegel_A)*2;% + CwF_Fluegel_R;
% 
% %% Rumpf
% Ind_W.b_rumpf = specs.D_rumpf;
% Ind_W.h_rumpf = specs.D_rumpf;
% Ind_W.r_rumpf_faktor = 1; 
% 
% Ind_W.CwF_Rumpf = 0.0031 * Ind_W.r_rumpf_faktor * specs.l_rumpf * (Ind_W.b_rumpf + Ind_W.h_rumpf); % Formel 6 
% 
% %% Leitwerk
% 
% % Ind_W.CwF_LWT = (Ind_W.CwF_Rumpf + Ind_W.CwF_Fluegel) * 0.24;  % 1.24
% Ind_W.CwF_LWT = 1.24; % Annahme aud Aufgabenstellung
% %% TW
%     % Formfaktor r_n, berücksichtigt Widerstand von Pylonen sowie Interferenz
%         %  1,50  alle Triebwerke in Gondeln am Flügel 
%         %  1,65  dito, jedoch ein Triebwerk im Rumpfheck
%         %  1,25  Triebwerke in Gondeln am Rumpfheck
%         %  1,00  interne Triebwerke im Rumpfeinlauf
%         %  0,30  Triebwerke in der Flügelwurzel 
%     % Formfaktor für Schubumkehr
%         %  1,00  mit Schubumkehranlage
% 
% % Muss nochmal überprüft werden Ind_W.S_0
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ind_W.S_0 = startschub.S0_GTo_To(1,startschub.Startstrecke) * G_To; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Ind_W.PSI = 30; % (24s - 40s) [s] %%%%%%% kein plan !!!!!!!!!!! nachfragen
% 
% Ind_W.CwF_TW = 1.72 * 1.5 * 1 * ((5 + specs.bypass)/(1+specs.bypass)) * ((Ind_W.S_0)/(ISA.p0 * Ind_W.PSI));
% 
% %% Fahrwerk
% Ind_W.r_uc = 1; % Formfaktor Fahrwerk 
%                 %  1,00  voll einziehbares Fahrwerk ohne Verkleidung
%                 %  1,03  Hauptfahrwerk in PTL-Triebwerksgondel
%                 %        eingezogen
%                 %  1,08  Hauptfahrwerk in Rumpfgondel eingezogen
%                 %        (Transall)
%                 %  1,25  nicht einziehbares Fahrwerk
%                 %  1,35  nicht einziehbares Fahrwerk, unverkleidet
% 
% %% Cw0 
% 
% N_W.r_re = 47 * Ind_W.Re_crit^(-0.2); % const
% 
% N_W.Cw0F = N_W.r_re * Ind_W.r_uc * ( Ind_W.CwF_LWT * (Ind_W.CwF_Fluegel + Ind_W.CwF_Rumpf) + Ind_W.CwF_TW);
% N_W.Cw0 = N_W.Cw0F / ((DT.F_I + DT.F_A )*2);%+ F_R ); %% Frage die ganze Flügelfläche(F/2) oder nur die Umspülte Flügelfläche 
% % Cw0 = 0.016
% 
% %% CA Werte
% 
% 
% Widerstand.y_CR = 0:0.001:0.85; % gewählter wert aus profilPDF Seite 4     %% REISEFLUG -> CA Anpassen an maxmialen Wert CA_max
% oswald_clean = 0.9; % von 0.85 bis 0.9
% d_cW_compr_clean = 0;
% 
% oswald_comp_LR = 0.8;
% d_cW_compr_LR = 0.0005;
% 
% oswald_comp_HS = 0.75;
% d_cW_compr_HS = 0.002;
% 
% 
% for x = 1 : length(Widerstand.y_CR)       
%  % CLEAN
%     Widerstand.C_w_clean = N_W.Cw0 + ((Widerstand.y_CR(x)^2)/...
%         (pi*Ergebnisse_Fluegel.streckung_phi25_max * oswald_clean)) + d_cW_compr_clean; 
%     Widerstand.C_w_clean_all(x) = Widerstand.C_w_clean;
%  % COMPRSIBILITY LR
%     Widerstand.C_w_LR = N_W.Cw0 + ((Widerstand.y_CR(x)^2)/...
%         (pi*Ergebnisse_Fluegel.streckung_phi25_max*oswald_comp_LR)) + d_cW_compr_LR;
%     Widerstand.C_w_LR_all(x) = Widerstand.C_w_LR;
%  % COMPRESIBILITY HIGH SPEED
%     Widerstand.C_w_HS = N_W.Cw0 + ((Widerstand.y_CR(x)^2)/...
%         (pi*Ergebnisse_Fluegel.streckung_phi25_max*oswald_comp_HS)) + d_cW_compr_HS; 
%     Widerstand.C_w_HS_all(x) = Widerstand.C_w_HS;
% 
% end
% 
% %% LANGSAMFLUG / Lsndeanflug / Takeoff
% 
% %for x = 1 : length(y)       %CLEAN zum Vergleich
% 
% %    C_w_clean = Cw0 + ((y(x)^2)/(pi*10.9*0.9))+0;
% %    C_w_clean_all(x) = C_w_clean;
% %end
% 
% %plot(C_w_clean_all,y)
% Widerstand.y_to = 0:0.001:startschub.c_A_max_thrust_match;
% 
% oswald_comp_TO = 0.85; % 0.8 - 0.85
% d_cW_compr_TO = 0.02; % 0.01 bis 0.02
% 
% d_cW_compr_GD = 0.0115; % 0.0115 - 0.025
% 
% for x = 1 : length(Widerstand.y_to)       %TO CLEAN ohne Fahrwerk
% 
%     Widerstand.C_w_TO_clean = N_W.Cw0 + ((Widerstand.y_to(x)^2)/(pi * Ergebnisse_Fluegel.streckung_phi25_max * oswald_comp_TO)) + d_cW_compr_TO;
%     Widerstand.C_w_TO_clean_all(x) = Widerstand.C_w_TO_clean;
%     %TO mit Fahrwerk
%     Widerstand.C_w_TO = N_W.Cw0 + ((Widerstand.y_to(x)^2)/(pi * Ergebnisse_Fluegel.streckung_phi25_max * oswald_comp_TO)) + d_cW_compr_TO + d_cW_compr_GD;
%     Widerstand.C_w_TO_all(x) = Widerstand.C_w_TO;
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% LANGSAMFLUG LANDING
% 
% Widerstand.y = 0:0.001:landeanvorderung.c_A_max_LDG;
% oswald_comp_LDG = 0.8; % 0.75 bis 0.8
% d_cW_compr_LDG = 0.065; % 0.055 - 0.065
% 
% for x = 1 : length(Widerstand.y)       %LDG ohne Fahrwerk
% 
%     Widerstand.C_w_LDG_clean = N_W.Cw0 + ((Widerstand.y(x)^2)/(pi * Ergebnisse_Fluegel.streckung_phi25_max * oswald_comp_LDG)) + d_cW_compr_LDG;
%     Widerstand.C_w_LDG_clean_all(x) = Widerstand.C_w_LDG_clean;
%     %LDG mit Fahrwerk
%     Widerstand.C_w_LDG = N_W.Cw0 + ((Widerstand.y(x)^2)/(pi * Ergebnisse_Fluegel.streckung_phi25_max * oswald_comp_LDG)) + d_cW_compr_LDG + d_cW_compr_GD;
%     Widerstand.C_w_LDG_all(x) = Widerstand.C_w_LDG;
% end
% 
% %%%%%%%%%%%%%%%% REZIPROKE GLEITZAHLEN! 
% 
% 
% % E = CA / CW
% GZ.CA_CW_Clean = Widerstand.y_CR./Widerstand.C_w_clean_all;
% GZ.CA_CW_LR = Widerstand.y_CR./Widerstand.C_w_LR_all;
% GZ.CA_CW_HS = Widerstand.y_CR./Widerstand.C_w_HS_all;
% GZ.CA_CW_TO_Clean = Widerstand.y_to./Widerstand.C_w_TO_clean_all;
% GZ.CA_CW_TO = Widerstand.y_to./Widerstand.C_w_TO_all;
% GZ.CA_CW_LDG_Clean = Widerstand.y./Widerstand.C_w_LDG_clean_all;
% GZ.CA_CW_LDG = Widerstand.y./Widerstand.C_w_LDG_all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




figure(2)
hold on
grid on
title('Lilienthalpolare')
plot(Widerstand.C_w_clean_all,Widerstand.y_CR,'k')
plot(Widerstand.C_w_LR_all,Widerstand.y_CR,'b')
plot(Widerstand.C_w_HS_all,Widerstand.y_CR,'--k')
legend('Cruise Clean Configuration','Cruise Long Range','Cruise High Speed','Location','southeast')

xlabel('Widerstandsbeiwert C_w')
ylabel('Auftriebsbeiwert C_A')
hold off



figure(3)
hold on
grid on
title('Lilienthalpolare')
plot(Widerstand.C_w_clean_all,Widerstand.y_CR,'k')
plot(Widerstand.C_w_TO_clean_all,Widerstand.y_to)
plot(Widerstand.C_w_TO_all,Widerstand.y_to)

% Landung

plot(Widerstand.C_w_LDG_clean_all,Widerstand.y)
plot(Widerstand.C_w_LDG_all,Widerstand.y)

legend('Cruise Clean Configuration','TO Clean Configuration','TO Gear Down','LD Clean Configuration','LD Gear Down', 'Location','southeast')

xlabel('Widerstandsbeiwert C_w')
ylabel('Auftriebsbeiwert C_A')
ylim([0 2.7])
hold off

figure(4)
hold on
grid on
title('Reziproke Gleitzahl')
plot(Widerstand.y_CR,GZ.CA_CW_Clean, 'b');
plot(Widerstand.y_CR,GZ.CA_CW_LR,'b--');
plot(Widerstand.y_CR,GZ.CA_CW_HS,'b-.');
plot(Widerstand.y_to,GZ.CA_CW_TO_Clean,'r');
plot(Widerstand.y_to,GZ.CA_CW_TO,'r--');
plot(Widerstand.y,GZ.CA_CW_LDG_Clean,'m');
plot(Widerstand.y,GZ.CA_CW_LDG,'m--')
xlabel("Auftriebsbeiwert C_A")
ylabel("Reziproke Gleitzahl E")

xlim([0 2.7])


plot(Ergebnisse_stat_Flaechenbelastung.C_A_CR, schub_CR.Eta, 'oblack','MarkerSize', 8)
plot(startschub.c_A_max_thrust_match, startschub.Eta_To_inv(3,1),'xblack','MarkerSize', 8)
plot(landeanvorderung.c_A_max_LDG,(1/landeanvorderung.Eta_LDG),'*black','MarkerSize', 8)

plot(Widerstand.y_CR(1,(round(Ergebnisse_stat_Flaechenbelastung.C_A_CR * 10^3))),...
    GZ.CA_CW_LR(1,(round(Ergebnisse_stat_Flaechenbelastung.C_A_CR * 10^3))),'or','MarkerSize', 8) % 523
plot(Widerstand.y_to(1,(round(startschub.c_A_max_thrust_match * 10^3))),...
    GZ.CA_CW_TO(1,(round(startschub.c_A_max_thrust_match * 10^3))),'xr','MarkerSize', 8) % 1801
plot(Widerstand.y(1,(round(landeanvorderung.c_A_max_LDG * 10^3))),...
    GZ.CA_CW_LDG(1,(round(landeanvorderung.c_A_max_LDG * 10^3))),'*red','MarkerSize', 8) % 2401

legend("Cruise Clean","Cruise LR","Cruise HS","TO Clean","TO Gear Down","LDG clean","LDG Gear Down",...
    'Reale Gleitzahl CR','Reale Gleitzahl TO','Reale Gleitzahl LDG','Ideale Gleitzahl CR','Ideale Gleitzahl TO','Ideale Gleitzahl LDG')

hold off













%% Vergleichen mit Schubanforderung FE1































 
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
