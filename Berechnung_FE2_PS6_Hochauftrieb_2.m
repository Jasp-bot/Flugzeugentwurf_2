function Berechnung_FE2_PS6_Hochauftrieb_2

clc
clear all
close all

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
load Schwerpunkt.mat
load Ergebnisse_Widerstand_FE2.mat
load Fahrwerk.mat
load Ergebnisse_Massen_FE2.mat

%% Aus FE 1 -> Für Widerstand

load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat
load Ergebnisse_Basis_stat_m.mat
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat
load Ergebnisse_Start_Landeanforderungen.mat
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Widerstand.mat
addpath('Unterfunktionen Widerstand');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Steuervariablen für Klappen

%Spannweite Klappenfläche -> Ab Rumpfmittellinie
spannweite_flaps = 0.65;

%Länge Klappe ausgefahren -> Quelle ACAMP + A350 Bilder + Chatgpt
flap_length_LDG = 2.0; % in meter
flap_length_TO = 1.0;

%Tiefe Klappen
Flaps_begin = 0.67; %Prozent
Flaps_begin_fahrwerk_modifizierung = 0.86; % Damit Fahrwerk nicht Klappen schneidet


%Faktoren Slat -> Lock der Slat Länge damit auf 70% eta -> Letzter Wert hir
%ist eta
Slat_spannweite = 0.85;

faktoren_Slat = 0.9 * 1 * Slat_spannweite;
% Für Winkeleinfluss Winkel 20-30 / Optimal -> Verhältnis Wert
% Slat Spannweite bis 90% . Länge ist 15 % -> bei 15 % ablesen wert im
% Letzten Graph

%Slats Tiefe
Slats_pos = 0.1;

%Oswald Zahl
oswald = 0.8;

cheatfaktor = 1.0;

%extra_rumpf_klappe_abstand
extra_rumpf_klappe_abstand = 0.02;

%% Werte aus anderen PS/Abgelesen

% Ergebnisse Schwerpunkt
r_h = r_H;
deltaXSP_l_mue = Delta_CG_MAC_durch_lmue;

%Ergebnisse Flügel
l_mue = Ergebnisse_Fluegel.l_mue;
CM0 = FM.c_M_NP_F0;%



% CW_REYNOLDS_MIN -> Über trapz summieren!
CA_index = round(Ergebnisse_Widerstand_FE2.stuetzstellen * 0.207);
CWPMIN_temp = Ergebnisse_Widerstand_FE2.c_w_p_min_Re_off_D(20,:);       % Hier 20 weil es dann am Ende passt !
C_W_P_Min_RE = trapz(CWPMIN_temp)* 10^-3; % Zwei mal wegen zwei Wings?!
%C_W_P_Min_RE = 0.009


% Reifen -> Aus PS Fahrwerk in m
durchmesser = Gear.durchmesser / 39.37 ;
breite = Gear.breite / 39.37 ;

%Momentberechnung
lambda       = Ergebnisse_Fluegel.lambda;   %Zuspitzung des Flügels [-]
bf_s         = spannweite_flaps;    %prozentuale Spannweite der Hinterkantenklappen


%% Klappenfläche
% Berechnung F_K
% Tiefdecker -> Länge Klappen und Abstände
b_k_a = 2 * ((Ergebnisse_Fluegel.b/2) * spannweite_flaps); % Auswaählen bis welches ETA!
b_k_i = specs.D_rumpf + extra_rumpf_klappe_abstand * (Ergebnisse_Fluegel.b/2); %  3% Abstand zum Rump über Halbspannweite

% Fläche Berechnen Ganzer Flügel für klappen mit Rumpf
X = linspace(0,spannweite_flaps,spannweite_flaps*1000);
Fluegel = Ergebnisse_Fluegel.Fluegeltiefen_eta(1,1:(spannweite_flaps*1000));
F_k = 2 * (Ergebnisse_Fluegel.b/2)*trapz(X,Fluegel); % MIT RUMPF!



% Rumpf + extra stücke 3%
Fluegel2 = Ergebnisse_Fluegel.Fluegeltiefen_eta_oR(1,1:extra_rumpf_klappe_abstand*1000);
X2 = linspace(0,extra_rumpf_klappe_abstand,extra_rumpf_klappe_abstand*1000);
F_rumpf_rechteck = Ergebnisse_Fluegel.Fluegeltiefen_eta_Ru(1) * specs.D_rumpf;

F_rumpf_Ecken = 2 * ((Ergebnisse_Fluegel.b)/2)*trapz(X2,Fluegel2);  % Rumpf Ecken

F_rumpf_dreieck = Flaeche_im_Rumpf_oberes_dreieck;              % Dreieck im Rumpf aus CG von Ole

F_rumpf = F_rumpf_Ecken + F_rumpf_rechteck %+F_rumpf_dreieck ;   % CHANGE: Dreieck entfernt, weil in den Flügel etas fällt das dreieck eh -> Fläche vergrößert -> Lunch for free!

% Somit F_k - F Rumpf = Klappenfläche

F_klappen = F_k - F_rumpf;      % Wahre Klappenfläche -> Tiefdecker Rechenvariante


%% Mittlere Flügeltiefe
% Neue Mittlere Flügeltiefe berechnen -> Für Außenflügel und Rumpf-Kink
laenge_Fluegel = length(Ergebnisse_Fluegel.Fluegeltiefen_eta);

ende = round(laenge_Fluegel * spannweite_flaps);

ende_klappe_tiefe = Ergebnisse_Fluegel.Fluegeltiefen_eta(ende);

anfang_kink = 0.3;

anfang = round(length(Ergebnisse_Fluegel.Fluegeltiefen_eta) * anfang_kink);

anfang_klappen_tiefe = Ergebnisse_Fluegel.Fluegeltiefen_eta(anfang);

eta_rumpf_3_prozent = 0.09 + extra_rumpf_klappe_abstand;

anfang_klappen_tiefe_abrumpf = round(length(Ergebnisse_Fluegel.Fluegeltiefen_eta) * eta_rumpf_3_prozent);
klappen_tiefe_abrumpf = Ergebnisse_Fluegel.Fluegeltiefen_eta(anfang_klappen_tiefe_abrumpf);



mittlere_tiefe_klappen = (((anfang_klappen_tiefe + ende_klappe_tiefe)/2) + (klappen_tiefe_abrumpf + ende_klappe_tiefe)/2)/2;



%% LANDING

% Kleine Formel CA_F_MAX -> AUS PS05
HA1.CA_F_max; %-> Übernehmen aus Hochauftrieb 1

%  Formel deltaCAFMAX,SF,phi -> Seite 33 Üb
delta_Ca_max_SF_phi = 1.57;  % Ablesen aus Grafik bei Klappenausschlag ~45°    S.34 Üb
delta_CA_F_max_SF_phi = delta_Ca_max_SF_phi * (F_klappen/Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max^2);

%Annahme eta_opt für slat = 55°
% Vorflügel           % 0.93 = Faktor für Vorflügel/Slat        % Annahme
% 74% des Flügels mit Slats belegt      Annahme Slats tiefe = 1m -> 1/l_m=6,57 = 0.152

delta_CA_F_max_VF = 0.93 * faktoren_Slat * (cos(rad2deg(Ergebnisse_Fluegel.phi_25_max)-5))^2;       % S.39 Üb


% Formel S.38
CA_F_max_VFFK = HA1.CA_F_max + delta_CA_F_max_SF_phi + delta_CA_F_max_VF;

% Neuer Auftriebsanstieg
% Hier erweitern auf realistischer -> Flügeltiefe bis Kink mit einem Wert
% -> Anstieg berechnen + Rest vom Flügel mit eigener mittleren Flügeltiefe
% von kink bis querruder!!



c = mittlere_tiefe_klappen;%Ergebnisse_Fluegel.l_m;

c_ = flap_length_LDG + mittlere_tiefe_klappen;%Ergebnisse_Fluegel.l_m;


%phi_50 = atan(tan(Ergebnisse_Fluegel.phi_25_max)-(4/Ergebnisse_Fluegel.streckung_phi25_max)* (0.5-0.25) * (1-Ergebnisse_Fluegel.lambda)/(1-Ergebnisse_Fluegel.lambda));
%CAalpha_F = (pi * Ergebnisse_Fluegel.streckung_phi25_max) / (1+sqrt(1 + ((Ergebnisse_Fluegel.streckung_phi25_max/2)^2) * (tan(phi_50)^2 + (1-specs.Ma_CR^2))));

CA_alpha_F_FK_phi = HA1.CA_alpha_F * (((c_/c)-1) * (F_klappen/Ergebnisse_Fluegel.F)+1);


%delta_CA_F_SF_phi
%c_k = c - (Ergebnisse_Fluegel.l_m * Flaps_begin);
%test = c_k / c; % -> Für Landing ist dann der Faktor = 1.84

delta_C_a_FK = 1.84;
delta_CA_F_SF_phi = (F_klappen / Ergebnisse_Fluegel.F) * (cos(Ergebnisse_Fluegel.phi_25_max^2)) *  delta_C_a_FK;

% Große Formel VFFK -> Gobbinsche Version
%In Übungfolien ist
%hier ein +                                                                             %+
alpha_F_max_VFFK = (CA_F_max_VFFK / CA_alpha_F_FK_phi) - (((HA1.CA_bei6_deg + delta_CA_F_SF_phi)) / CA_alpha_F_FK_phi) + deg2rad(6) + deg2rad(HA1.delta_alpha_CA_F_max) + deg2rad(HA1.alpha_MAC_F_0);

alpha_F_max_VFFK_deg = rad2deg(alpha_F_max_VFFK);

%% TAKEOFF

% kleine Formel deltaCAFMAX,SF,phi
delta_Ca_max_SF_phi_TO = 1.1;  % Ablesen aus Grafik bei Klappenausschlag ~20°

delta_CA_F_max_SF_phi_TO = delta_Ca_max_SF_phi_TO * (F_klappen/Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max^2);

%Annahme eta_opt für slat = 55°
% Vorflügel           % 0.93 = Faktor für Vorflügel/Slat        % Annahme
% 74% des Flügels mit Slats belegt      Annahme Slats tiefe = 1m -> 1/l_m=6,57 = 0.152
delta_CA_F_max_VF_TO = 0.93 * faktoren_Slat * (cos(rad2deg(Ergebnisse_Fluegel.phi_25_max)-5)^2); %0.93 * 0.64 * 0.69 * 1.0 * (cos(Ergebnisse_Fluegel.phi_25_max)^2);
%Rad oder Degree ?


% MAximaler Auftriebsbeiwert des Flügels mit Fowler, Slats

CA_F_max_VFFK_TO = HA1.CA_F_max + delta_CA_F_max_SF_phi_TO + delta_CA_F_max_VF_TO;

% Neuer Auftriebsanstieg
c_TO = mittlere_tiefe_klappen; %Ergebnisse_Fluegel.l_m; % Mittlere Flügeltiefe benutzen
c__TO = flap_length_TO + mittlere_tiefe_klappen;%Ergebnisse_Fluegel.l_m; % Quelle GPT und ACAMP

CA_alpha_F_FK_phi_TO = HA1.CA_alpha_F * (((c__TO/c_TO)-1) * (F_klappen/Ergebnisse_Fluegel.F) + 1);


%delta_CA_F_SF_phi
%c_k = c - (Ergebnisse_Fluegel.l_m * 0.65);

% Für Takeoff ist dann der Faktor = 1.3
delta_C_a_FK_TO = 1.3;
delta_CA_F_SF_phi_TO = (F_klappen / Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max^2) *  delta_C_a_FK_TO;

% Große Formel VFFK - >Gobbin Version
%IN ÜBUNG IST HIER
%+                                                                                           %+
alpha_F_max_VFFK_TO = (CA_F_max_VFFK_TO/CA_alpha_F_FK_phi_TO) - (((HA1.CA_bei6_deg + delta_CA_F_SF_phi_TO))/CA_alpha_F_FK_phi_TO) + deg2rad(6) + deg2rad(HA1.delta_alpha_CA_F_max) + deg2rad(HA1.alpha_MAC_0_F);

alpha_F_max_VFFK_deg_TO = rad2deg(alpha_F_max_VFFK_TO);


%% Alphas/ CAs plotting vorbereiten

alphas_st = -12:0.01:(alpha_F_max_VFFK_deg_TO - HA1.delta_alpha_CA_F_max + cheatfaktor ); % normal Plotten bis alphamax - delta alpha
CA_st = CA_alpha_F_FK_phi_TO.*(deg2rad(alphas_st-HA1.alpha_MAC_0)) + delta_CA_F_max_SF_phi_TO;


alphas_sl = -15:0.01:alpha_F_max_VFFK_deg-HA1.delta_alpha_CA_F_max+cheatfaktor; % normal Plotten bis alphamax - delta alpha
CA_sl = CA_alpha_F_FK_phi.*(deg2rad(alphas_sl-HA1.alpha_MAC_0)) + delta_CA_F_max_SF_phi;




%% Plot speichern

% Speichern

HA2.alpha_F_max_VFFK_deg = alpha_F_max_VFFK_deg;

HA2.CA_alpha_F_FK_phi = CA_alpha_F_FK_phi;

HA2.delta_CA_F_max_SF_phi = delta_CA_F_max_SF_phi;

HA2.alpha_F_max_VFFK_deg_TO = alpha_F_max_VFFK_deg_TO;

HA2.CA_alpha_F_FK_phi_TO = CA_alpha_F_FK_phi_TO;

HA2.delta_CA_F_max_SF_phi_TO = delta_CA_F_max_SF_phi_TO;

HA2.CA_F_max_VFFK_TO = CA_F_max_VFFK_TO;

HA2.CA_F_max_VFFK = CA_F_max_VFFK;

HA2.CA_sl = CA_sl;

HA2.CA_st= CA_st;

HA2.alphas_sl = alphas_sl;

HA2.alphas_st = alphas_st;





%% Prüfen ob CA max passt für alle Flugphasen GILT NUR FÜR FLÜGEL SELBER!!

% check ob StartCA erreicht
startschub.c_A_max_thrust_match < CA_F_max_VFFK_TO;

% Check ob lande CA erreicht
landeanvorderung.c_A_max_LDG < CA_F_max_VFFK;



%% Momentenänderung - Trimmung noch möglich?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Aus Formel.m
%nach Skript Teil D Seite 114 -> Braucht Lambda und bf_S
dcMk_dcmK  = ( 3.739701 * lambda^4 - 10.762986 * lambda^3 + 12.791164 * lambda^2 - 8.860305 * lambda + 4.093244) * bf_s + ...
    (-5.744622 * lambda^4 + 17.487598 * lambda^3 - 21.515763 * lambda^2 + 13.990786 * lambda - 4.219788) * bf_s^2 +...
    ( 1.960553 * lambda^4 -  7.330970 * lambda^3 + 10.252347 * lambda^2 -  6.347223 * lambda + 1.466302) * bf_s^3;

%nach Skript Teil D Seite 114
dcMk_dcAK =  (-0.254693 * lambda^3 + 0.343337 * lambda^2 - 0.172780 * lambda + 0.333119) * bf_s +...
    ( 0.986155 * lambda^3 - 1.620515 * lambda^2 + 1.049655 * lambda - 0.661691) * bf_s^2 + ...
    (-0.735627 * lambda^3 + 1.282365 * lambda^2 - 0.877425 * lambda + 0.328347) * bf_s^3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Formel 7 + 8

CA_REF_TO = CA_F_max_VFFK_TO /(1.2^2);
CA_REF_LDG = CA_F_max_VFFK / (1.3^2);

%Formel 10
%ck_1 = c - (c * Flaps_begin); Ersetzt durch siehe unten!

%% Anpassen an Klappentiefe -> Mittelwert Außenflügel und Mittelwert Innenflügel -> Gewichtet
%Formel 10
c1 = ende_klappe_tiefe - (ende_klappe_tiefe * Flaps_begin);

c2 = anfang_klappen_tiefe - (anfang_klappen_tiefe * Flaps_begin);

c3 = klappen_tiefe_abrumpf -(klappen_tiefe_abrumpf * Flaps_begin_fahrwerk_modifizierung);

cc1 = (c1+c2)/2;

cc2 = (c2+c3)/2;

ck_1 = (cc1*0.7+cc2*0.3)/(0.3+0.7);

%%

%Formel 10
theta = acos(2 * (ck_1/c) -1);

%Formel 9

dCM_dCA_deltaCA = 0.5 * (1 - (ck_1/c)) * ((sin(theta)) / (pi - (theta - sin(theta))));


%Formel 6
%delta Ca_FK
%bei Landing
%und Takeoff ->
%Abgelesen aus
%Grafik! ->
%Variable bei
%Polaren
%festgelegt
delta_Cm_HKK_LDG = - dCM_dCA_deltaCA * (c_/c) - ((CA_REF_LDG + delta_C_a_FK * (1-(F_klappen/Ergebnisse_Fluegel.F))) / (8)) * (c_/c) * ((c_/c)-1);
delta_Cm_HKK_TO  = - dCM_dCA_deltaCA * (c__TO/c_TO) - ((CA_REF_TO + delta_C_a_FK_TO * (1-(F_klappen/Ergebnisse_Fluegel.F))) / (8)) * (c__TO/c_TO) * ((c__TO/c_TO)-1);


%  Formel 5 für Takeoff und Landing
%*(CA_REF_LDG)                                                                                                                                    %SRichtig?
delta_CM_HKK_LDG = dcMk_dcmK * delta_Cm_HKK_LDG + spannweite_flaps * ((Ergebnisse_Fluegel.streckung_phi25_max)/(1+(2/Ergebnisse_Fluegel.streckung_phi25_max))) * dcMk_dcAK * delta_CA_F_SF_phi * tan(Ergebnisse_Fluegel.phi_25_max);
%*(CA_REF_TO)
delta_CM_HKK_TO = dcMk_dcmK * delta_Cm_HKK_TO + spannweite_flaps * ((Ergebnisse_Fluegel.streckung_phi25_max)/(1+(2/Ergebnisse_Fluegel.streckung_phi25_max))) * dcMk_dcAK * delta_CA_F_SF_phi_TO * tan(Ergebnisse_Fluegel.phi_25_max);


% Formel 4 - Final

CA_MAX_TO = ( CA_F_max_VFFK_TO + ( (CM0 + delta_CM_HKK_TO)/( r_h/l_mue ) ) ) / (1 - ( deltaXSP_l_mue )/( r_h/l_mue ));
CA_MAX_LDG = ( CA_F_max_VFFK   + ( (CM0 + delta_CM_HKK_LDG)/( r_h/l_mue ) ) ) / (1 - ( deltaXSP_l_mue )/( r_h/l_mue ));

HA2.CA_MAX_TO = CA_MAX_TO;
HA2.CA_MAX_LDG = CA_MAX_LDG;

% check ob StartCA erreicht
disp('Start CA wird erreicht?')
startanforderung = startschub.c_A_max_thrust_match < CA_MAX_TO

% Check ob lande CA erreicht
disp('Lande CA wird erreicht?')
landeanforderun = landeanvorderung.c_A_max_LDG < CA_MAX_LDG


%% Widerstandszuwachs durch Klappen


% Profilwiderstand % Unterschiedlich bei TO /LDG
delta_C_W_P_phi = 0.081;    %Abbildung + Abhängig von Auschlag delta_k
delta_C_W_P_phi_TO = 0.015;

delta_C_W_P = delta_C_W_P_phi * (F_klappen/Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max);
delta_C_W_P_TO = delta_C_W_P_phi_TO * (F_klappen/Ergebnisse_Fluegel.F) * cos(Ergebnisse_Fluegel.phi_25_max);

%Induzierter Widerstand % Unterschiedlich für TO/LDG

delta_Ca_max_SF_phi_TO;
delta_Ca_max_SF_phi;

% Klappenspannweite durch Spannweite
temp = 0.6;
Ergebnisse_Fluegel.lambda;
v = 0.001;                    % Richtig abgelesen?!?!
w = 0.0073;

delta_CM_K = delta_CM_HKK_LDG;%0.5; % Wie kommt man auf den Wert ?
delta_CM_K_TO = delta_CM_HKK_TO;


CA = HA2.CA_MAX_LDG;
CA_TO = HA2.CA_MAX_TO;

CA_F = CA * (1- ((deltaXSP_l_mue)/(r_h/Ergebnisse_Fluegel.l_mue))) - ((CM0+delta_CM_K)/(r_h/Ergebnisse_Fluegel.l_mue))/(r_h/Ergebnisse_Fluegel.l_mue);

CA_F_TO = CA_TO * (1- ((deltaXSP_l_mue)/(r_h/Ergebnisse_Fluegel.l_mue))) - ((CM0+delta_CM_K_TO)/(r_h/Ergebnisse_Fluegel.l_mue))/(r_h/Ergebnisse_Fluegel.l_mue);

delta_C_W_Ind = CA_F * delta_C_a_FK * v + (delta_Ca_max_SF_phi^2) *w;

delta_C_W_Ind_TO = CA_F_TO * delta_C_a_FK * v + (delta_Ca_max_SF_phi_TO^2) *w;


% Inteferenzwiderstand % Unterschiedlich bei TO/LDG


delta_C_W_Inf = (1/3) * delta_C_W_P;
delta_C_W_Inf_TO = (1/3) * delta_C_W_P_TO;

% Vorflügelwiderstand % Gleich für TO/LDG weil gleiche Settings gewählt

%Fläche Vorflügel ist von Holm vorne bis Vorderkante    % Über mittlere
%Flügeltiefe ok?

%tiefe_Slats = Ergebnisse_Fluegel.Fluegeltiefen_eta(1) * Slats_pos;

% Länge Slats
laenge_Slats = 2 * (Slat_spannweite * (Ergebnisse_Fluegel.b/2)); % m

X_VF = linspace(0,Slat_spannweite,Slat_spannweite*1000);
% Slat tiefen!
Fluegel_VF = Slats_pos .* Ergebnisse_Fluegel.Fluegeltiefen_eta(1,1:(Slat_spannweite*1000));
% Fläche Slats
F_VF = ((Ergebnisse_Fluegel.b/2)*trapz(X_VF,Fluegel_VF))*2;


laenge_Fluegel = Ergebnisse_Fluegel.b; %Ohne Rumpf?

% A*b ansatz ? Genau genug
%F_VF = tiefe_Slats * laenge_Slats * 0.9;        % Berechnug über mittlere Flügeltiefe nicht sehr genau -> Flcähe etwas größer als Echt -> Faktor 0.9
%F_VF = 1.93+5.26;


delta_CW_VF = C_W_P_Min_RE * (F_VF/Ergebnisse_Fluegel.F) * (laenge_Slats/(laenge_Fluegel-specs.D_rumpf)) * cos(Ergebnisse_Fluegel.phi_25_max);


% Fahrwerkswiderstand % Bleibt gleich bei TO/LDG

% Braucht finale Werte aus Fahrwerk
F_vorder = durchmesser*breite*4;   % Gleiche Werte weil vorne und hinten gleich groß sind
F_hinter = durchmesser*breite*4;

l_HFW = BFWL.l_HFW_min; % Ríchtiger Wert aus CG?
delta_CA_F_0_LDG = -0.7;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_CA_F_0_TO = -1;       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALEX NOCHMAL FRAGEN ODER KROISTOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


delta_C_W_Fahrwerk_TO = ((1.5 * F_vorder + 0.75 * F_hinter)/Ergebnisse_Fluegel.F) * (1 - 0.04 * ((CA_F_TO + delta_CA_F_0_TO *(1.5 * (Ergebnisse_Fluegel.F / F_klappen)-1))/(l_HFW/Ergebnisse_Fluegel.l_m)))^2;

delta_C_W_Fahrwerk_LDG = ((1.5 * F_vorder + 0.75 * F_hinter)/Ergebnisse_Fluegel.F) * (1 - 0.04 * ((CA_F + delta_CA_F_0_LDG *(1.5 * (Ergebnisse_Fluegel.F / F_klappen)-1))/(l_HFW/Ergebnisse_Fluegel.l_m)))^2;

%delta_C_W_Fahrwerk = 0.017; % CHEAT FAKTOR!

%% Gesamtwiderstand durch klappen

delta_CW_klappen_LDG = delta_CW_VF + delta_C_W_Ind + delta_C_W_Inf + delta_C_W_P;
delta_CW_klappen_TO = delta_CW_VF + delta_C_W_Ind_TO + delta_C_W_Inf_TO + delta_C_W_P_TO;

delta_CW_klappen_LDG_fahrwerk = delta_CW_VF + delta_C_W_Ind + delta_C_W_Inf + delta_C_W_P + delta_C_W_Fahrwerk_LDG;
delta_CW_klappen_TO_fahrwerk = delta_CW_VF + delta_C_W_Ind_TO + delta_C_W_Inf_TO + delta_C_W_P_TO + delta_C_W_Fahrwerk_TO;


% Speichern

HA2.delta_CW_klappen_LDG = delta_CW_klappen_LDG;

HA2.delta_CW_klappen_TO = delta_CW_klappen_TO;

HA2.delta_CW_klappen_LDG_fahrwerk = delta_CW_klappen_LDG_fahrwerk;

HA2.delta_CW_klappen_TO_fahrwerk = delta_CW_klappen_TO_fahrwerk;


%% Widerstanddeltas zu realem Widerstand addieren

v_eingang = landeanvorderung.v_50;
hoehe_LDG = round(convlength(1500,'ft','m'));

% Clean Fall wird in Visualisierung basierend auf FE2 Widerstand berechnet! 

% % Clean Fall bis CA aus HA1
% %downsampling
% originalVector = HA1.CAs;
% desiredLength = 1001;
% shortenedVector = linspace(originalVector(1), originalVector(end), desiredLength);
% HA2.c_A_F_clean = shortenedVector;
% 
% [x_vector_sum, x_vector] = Landung(v_eingang, hoehe_LDG, HA2.c_A_F_clean);
% 
% HA2.clean_CW = x_vector_sum(8,:);
% 
% [~,a,~,~,~] = atmosisa(11000);
% vv = specs.Ma_CR*a;
% 
% [x_vector_sum, x_vector] = Landung(vv, 11000, HA2.c_A_F_clean);
% HA2.clean_CW = x_vector_sum(8,:);
% %HA2.c_A_F_clean = Ergebnisse_Widerstand_FE2.c_A_ges;


%Takeoff ohne Fahrwerk
%downsampling
originalVector = HA2.CA_st;
desiredLength = 1001;
shortenedVector = linspace(originalVector(1), originalVector(end), desiredLength);
HA2.c_A_F_TO = shortenedVector;

[x_vector_sum, x_vector] = Landung(v_eingang, hoehe_LDG, HA2.c_A_F_TO);

HA2.TO_CW = x_vector_sum(8,:) + delta_CW_klappen_TO; 



%Takeoff mit Fahrwerk
%downsampling
originalVector = HA2.CA_st;
desiredLength = 1001;
shortenedVector = linspace(originalVector(1), originalVector(end), desiredLength);
HA2.c_A_F_TO_FW = shortenedVector;

[x_vector_sum, x_vector] = Landung(v_eingang, hoehe_LDG, HA2.c_A_F_TO_FW);

HA2.TO_CW_FW = x_vector_sum(8,:) + delta_CW_klappen_TO_fahrwerk; 




%Landing ohne Fahrwerk
%downsampling
originalVector = HA2.CA_sl;
desiredLength = 1001;
shortenedVector = linspace(originalVector(1), originalVector(end), desiredLength);
HA2.c_A_F_LDG = shortenedVector;

[x_vector_sum, x_vector] = Landung(v_eingang, hoehe_LDG, HA2.c_A_F_LDG);

HA2.LDG_CW = x_vector_sum(8,:) + delta_CW_klappen_LDG; 



%Landing mit Fahrwerk
%downsampling
originalVector = HA2.CA_sl;
desiredLength = 1001;
shortenedVector = linspace(originalVector(1), originalVector(end), desiredLength);
HA2.c_A_F_LDG_FW = shortenedVector;

[x_vector_sum, x_vector] = Landung(v_eingang, hoehe_LDG, HA2.c_A_F_LDG_FW);

HA2.LDG_CW_FW = x_vector_sum(8,:) + delta_CW_klappen_LDG_fahrwerk; 


%% Massenberechnung Fowler system -> wahre Klappenfläche bestimmen
%Über Jeweillige Mittlere Flügeltiefe -> Stimmt sehr gut mit CAD überein
c1 = ende_klappe_tiefe - (ende_klappe_tiefe * Flaps_begin);

c2 = anfang_klappen_tiefe - (anfang_klappen_tiefe * Flaps_begin);

c3 = klappen_tiefe_abrumpf -(klappen_tiefe_abrumpf * Flaps_begin_fahrwerk_modifizierung);

cc1 = (c1+c2)/2;
% Von eta 0,3 bis eta Klappenspannweite
spannweite_aussen = spannweite_flaps - anfang_kink;
spannweite_aussen_laenge = spannweite_aussen * (Ergebnisse_Fluegel.b/2);

Area_aussen = spannweite_aussen_laenge * cc1;

%Ab eta 0,3 bis Eta ~0,11
cc2 = (c2+c3)/2;
spannweite_innen = (anfang_kink - anfang_klappen_tiefe_abrumpf / 1000) * (Ergebnisse_Fluegel.b/2);

Area_innen = spannweite_innen * cc2;

HA2.F_Fowler = 2 * (Area_innen + Area_aussen);

% Für Massenberechnung tatsächliche Fowler klappen Fläche
% X_FOWLER = linspace(0.03,spannweite_flaps,spannweite_flaps*1000);
% % Slat tiefen!
% Fluegel_FOWLER = (1-Flaps_begin) .* Ergebnisse_Fluegel.Fluegeltiefen_eta(1,1:(spannweite_flaps*1000));
% % Fläche Slats
% F_FOWLER = 2 * ((Ergebnisse_Fluegel.b/2)*trapz(X_FOWLER,Fluegel_FOWLER));




%% Ausgabe für Jasper
HA2.CW_max_TO = HA2.TO_CW(end);
HA2.CA_max_TO = HA2.c_A_F_TO(end);

HA2.CW_max_ldg_fw = HA2.LDG_CW_FW(end);
HA2.CA_max_ldg_fw = HA2.c_A_F_LDG_FW(end);



%% Speichern
save Ergebnisse_Hochauftrieb_2.mat HA2 spannweite_flaps Flaps_begin flap_length_LDG















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% REZIPROKE GLEITZAHLEN
%% E = C_A / C_W

% E_Clean = HA1.CAs./CW_clean;
% 
% E_TO = CA_st./CW_TO;
% E_TO_FW = CA_st./CW_TO_FW;
% 
% E_LDG = CA_sl./CW_LDG;
% E_LDG_FW = CA_sl./CW_LDG_FW;



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



%% Cruise -> Nicht Nötig
% figure(2)
% hold on
% grid on
% title('Lilienthalpolare')
% plot(Widerstand.C_w_clean_all,Widerstand.y_CR,'k')
% plot(Widerstand.C_w_LR_all,Widerstand.y_CR,'b')
% plot(Widerstand.C_w_HS_all,Widerstand.y_CR,'--k')
% legend('Cruise Clean Configuration','Cruise Long Range','Cruise High Speed','Location','southeast')
%
% xlabel('Widerstandsbeiwert C_w')
% ylabel('Auftriebsbeiwert C_A')
% hold off

%% Landing und TO

% figure(3)
% hold on
% grid on
% title('Lilienthalpolare')
% plot(Widerstand.C_W_TO_clean,Widerstand.y_CR,'k')
% Widerstand.C_W_TO_FW
% plot(Widerstand.C_w_TO_clean_all,Widerstand.y_CR)
% plot(Widerstand.C_w_TO_all,Widerstand.y_to)
%
% % Landung
%
% plot(Widerstand.C_w_LDG_clean_all,Widerstand.y)
% plot(Widerstand.C_w_LDG_all,Widerstand.y)
%
% legend('Cruise Clean Configuration','TO Clean Configuration','TO Gear Down','LD Clean Configuration','LD Gear Down', 'Location','southeast')
%
% xlabel('Widerstandsbeiwert C_w')
% ylabel('Auftriebsbeiwert C_A')
% ylim([0 2.7])
% hold off
%
% figure(4)
% hold on
% grid on
% title('Reziproke Gleitzahl')
% plot(Widerstand.y_CR,GZ.CA_CW_Clean, 'b');
% plot(Widerstand.y_CR,GZ.CA_CW_LR,'b--');
% plot(Widerstand.y_CR,GZ.CA_CW_HS,'b-.');
% plot(Widerstand.y_to,GZ.CA_CW_TO_Clean,'r');
% plot(Widerstand.y_to,GZ.CA_CW_TO,'r--');
% plot(Widerstand.y,GZ.CA_CW_LDG_Clean,'m');
% plot(Widerstand.y,GZ.CA_CW_LDG,'m--')
% xlabel("Auftriebsbeiwert C_A")
% ylabel("Reziproke Gleitzahl E")
%
% xlim([0 2.7])
%
%
% plot(Ergebnisse_stat_Flaechenbelastung.C_A_CR, schub_CR.Eta, 'oblack','MarkerSize', 8)
% plot(startschub.c_A_max_thrust_match, startschub.Eta_To_inv(3,1),'xblack','MarkerSize', 8)
% plot(landeanvorderung.c_A_max_LDG,(1/landeanvorderung.Eta_LDG),'*black','MarkerSize', 8)
%
% plot(Widerstand.y_CR(1,(round(Ergebnisse_stat_Flaechenbelastung.C_A_CR * 10^3))),...
%     GZ.CA_CW_LR(1,(round(Ergebnisse_stat_Flaechenbelastung.C_A_CR * 10^3))),'or','MarkerSize', 8) % 523
% plot(Widerstand.y_to(1,(round(startschub.c_A_max_thrust_match * 10^3))),...
%     GZ.CA_CW_TO(1,(round(startschub.c_A_max_thrust_match * 10^3))),'xr','MarkerSize', 8) % 1801
% plot(Widerstand.y(1,(round(landeanvorderung.c_A_max_LDG * 10^3))),...
%     GZ.CA_CW_LDG(1,(round(landeanvorderung.c_A_max_LDG * 10^3))),'*red','MarkerSize', 8) % 2401
%
% legend("Cruise Clean","Cruise LR","Cruise HS","TO Clean","TO Gear Down","LDG clean","LDG Gear Down",...
%     'Reale Gleitzahl CR','Reale Gleitzahl TO','Reale Gleitzahl LDG','Ideale Gleitzahl CR','Ideale Gleitzahl TO','Ideale Gleitzahl LDG')
%
% hold off
