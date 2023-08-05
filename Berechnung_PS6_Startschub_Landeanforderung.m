% berechnungen fuer Startschub und Landeanforderung

function Berechnung_PS6_Startschub_Landeanforderung(Startwerte_Iteration, Eingabewert_Iteration)

clc
%clear all
close all

%% Inputs

load Projekt_specs.mat;
load Ergebnisse_Basis_stat_m.mat;
load Ergebnisse_ISA_DATA.mat
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat
load Ergebnisse_Endwerte_Iteration_V1.mat
addpath('Unterfunktionen Widerstand');

%% Berechnungen

% Deklaration
Flughoehe = specs.flight_level * 10^2 ;                         % in ft
hoehe =round(unitsratio('m','ft')*Flughoehe);                     % in m

p_CR_B = ISA.p(hoehe,1);                                        % p_Basis Druckvector aus ISA+15
rho_CR_B = ISA.rho(hoehe,1);                                    % rho_Basis Dichtevector aus ISA+15

G_To = Ergebnis_basis_m.m_To * specs.g;

D = [1; 0.95; 0.9; 0.8];                                        % Drosselgrad TW

if Eingabewert_Iteration == 0
    load Ergebnisse_Widerstand.mat;
    % Dieser wert wird itteriert in PS10
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    schub_CR.Eta = Startwerte_Iteration.CA_CW_LR; %17.35;                               %[18; 19; 20; 21];       % mittlere erwartete Gelitzahl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Itteraion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    startschub.Eta_To = [1/8; 1/9; 1/Startwerte_Iteration.CA_CW_TO ];%10.558;];% 1/11; 1/12];
    startschub.Eta_To_thustmatch = 1/Startwerte_Iteration.CA_CW_TO; %10.558;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iteration 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Eta_LDG_vec = [1/6; 1/7.43; 1/8]; % gleitberhaeltnisse für LDG [3x1]
    Eta_LDG_vec = [1/6; 1/Startwerte_Iteration.CA_CW_LDG; 1/8];
    landeanvorderung.Eta_LDG = Eta_LDG_vec(2,1); % E = 7
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif Eingabewert_Iteration == 1
    clear Startwerte_Iteration;
    load Ergebnisse_Widerstand_FE2.mat;
    load Ergebnisse_Hochauftrieb_1.mat;
    load Ergebnisse_Hochauftrieb_2.mat;

    % CR
    schub_CR.Eta = 1/ Ergebnisse_Widerstand_FE2.cW_cA_off_D;
    % TO
    startschub.Eta_To_thustmatch = 1 / (HA2.CA_max_TO/ HA2.CW_max_TO);
    startschub.Eta_To = [1/8; 1/9; startschub.Eta_To_thustmatch];%10.558;];% 1/11; 1/12];
    % LDG
    landeanvorderung.Eta_LDG = 1/ (HA2.CA_max_ldg_fw/ HA2.CW_max_ldg_fw); % E = 7
    Eta_LDG_vec = [1/6; landeanvorderung.Eta_LDG; 1/8];
end

k_CR = 0.98;                                        % G_ICA / G_To


%% Rechnung fuer Cruise Zustand

% Berechnung aus Parameterstudie PS06 
% Formel 6 fuer S_S0_E 
% Formel 4 fuer S_S0_CR
% Formel 8 fuer S0_GTo_CR
% Formel 7 fuer S0_CR Erforderlicher Gesamtschub des zu intallierenden TW

% Schubverlust durch Einlass
schub_CR.S_S0_E = 1 - (1.3 + 0.25 * specs.bypass) * 0.02;


% Schubverhaeltnis Reiseflug
% schub_CR.S_S0_CR = D(4,1) .* (rho_CR_B/ISA.rho_0) .* exp(-0.35 * specs.Ma_CR * (p_CR_B/ISA.p0) * sqrt(specs.bypass));
schub_CR.S_S0_CR = S_S0_KF_j(D(4,1), (rho_CR_B/ISA.rho_0), specs.Ma_CR, (p_CR_B/ISA.p0), specs.bypass);
% Gesamtschubverhaeltnis
schub_CR.S0_GTo_CR =specs.Schubfaktor .* (k_CR .* 1./schub_CR.Eta)./(schub_CR.S_S0_CR .* schub_CR.S_S0_E);                                        % realistisch

% erforderlicher Gesamtschub
schub_CR.S_CR = G_To .* (k_CR * 1/schub_CR.Eta)./(schub_CR.S_S0_CR * schub_CR.S_S0_E) .* specs.Schubfaktor;


%% Rechnung für Startstrecke PS7
%% Startschub und Startstrecke s1

% Deklaration
a0 = ISA.a(1,1);                                        % Schallgeschwindigkeit auf Meereshoehe ISA+15
rho_To = ISA.rho(1,1);
p_To = ISA.p(1,1); % Pa 
startschub.c_A_max = [1.8; 1.9; 2]; % Vector für verschiedenne CA max 1.8; 1.9; 2; 
startschub.c_A_max_thrust_match = 1.8;
m_To_Flugzeug_thrust_match = Ergebnis_basis_m.m_To;
startschub.s1 = 1:1:2600;
reibung_boden = 0.02;

% Berechnungen
% Thurst matching veränderung von basis c_A_max(3,1) auf x
% v_LOF = 1.15* v_stall
% v_stall = v_s = sqrt((2 * m_To_Flugzeug_thrust_match * specs.g)./(rho_To .* c_A_max_thrust_match .* Ergebnisse_stat_Flaechenbelastung.F))
startschub.v_s = sqrt((2 * m_To_Flugzeug_thrust_match * specs.g)...
   ./(rho_To .* startschub.c_A_max_thrust_match .* Ergebnisse_stat_Flaechenbelastung.F));

startschub.v_LOF = 1.15 .* startschub.v_s;

startschub.Ma_To = (0.71 * startschub.v_LOF)/ a0; % Machzahl auf Meereshoehe bei ISA+15 für Startvorgang aus  PS Startstrecke und Landestrecke

% Startschubverhaeltnis fuer Takeoff
% startschub.S_S0_TO = D(1,1) .* (rho_To /ISA.rho_0) .* exp(-0.35 * startschub.Ma_To * (p_To/ISA.p0) * sqrt(specs.bypass)); 
startschub.S_S0_TO = S_S0_KF_j(D(1,1), (rho_To /ISA.rho_0), startschub.Ma_To, (p_To/ISA.p0), specs.bypass);

% Startschub/Gewichts verhaeltnis mit Startstrecke s1 (Formel 9)

term1 = 1.32 ./(specs.g .* rho_To .* (1 - 1./(2 .*specs.n_TW)) .*startschub.s1 .* startschub.c_A_max);
term2 = 1/((1-reibung_boden) * startschub.S_S0_TO * schub_CR.S_S0_E) * (m_To_Flugzeug_thrust_match * specs.g)/Ergebnisse_stat_Flaechenbelastung.F;
startschub.S0_GTo_To = term1 .* term2 .* specs.Schubfaktor; % Zusammensetzung der Formel um Startschubanforderung zu berechnen
% Ergebnis 6x2600 matrix zusammengesetzt aus s1 und c_A_max Werten

%% Steigstrecke s3

s3 = -1*(-2600:1:-1);

startschub.Eta_To_inv = 1./(startschub.Eta_To); % fuer legende
%c_A_max = [1.8; 1.9; 2;] % 2.1; 2.2; 2.25]; % Vector für verschiedenne CA max
h_2 = unitsratio('m','ft')*35;            % zu überfliegendes Hindernis (hoehe)

startschub.v_2 = 1.2 * startschub.v_s;
startschub.v_CL = (startschub.v_LOF + startschub.v_2) / 2; % mittlere Geschwindigkeit bei 35ft hoehe

startschub.Ma_CL = startschub.v_CL / a0; % mit der Annahme, dass in 35ft hoehe die gleichen bedingungen hersdchen wie bei h = 1m, sonst a(h_2)
% startschub.S_S0_CL = D(3,1) * (rho_To/ISA.rho_0) * exp(-0.35 * startschub.Ma_CL * (p_To/ISA.p0) * sqrt(specs.bypass)); % Verhaeltnis Schub zu Startschub Climb
startschub.S_S0_CL = S_S0_KF_j(D(3,1), (rho_To/ISA.rho_0), startschub.Ma_CL, (p_To/ISA.p0), specs.bypass);

% Ergebnis 5x2600 = s3 x Eta_To (Formel 14)
startschub.S0_GTo_CL = (sin(atan(h_2 ./ s3)) + startschub.Eta_To )./((1 - 1./specs.n_TW) .* startschub.S_S0_CL .* schub_CR.S_S0_E) .* specs.Schubfaktor;


%Startstrecke
startschub.intersections_1 = InterX([startschub.s1; startschub.S0_GTo_To(1,:)],[startschub.s1; startschub.S0_GTo_CL(3,:)]);
startschub.Startstrecke = round(startschub.intersections_1(1,1));

startschub.S0 = startschub.intersections_1(2,1) * G_To;
%% Landeanforderungen

% Deklaration von Werten

h_50 = unitsratio('m','ft')*50; % zu überfliegende Sicherheitshoehe bei Landung in m

b_m_vec = -(0.3:0.01:0.4)*specs.g; % Bremsverzögerung [1x11]
landeanvorderung.b_m = b_m_vec(1,11); % = -3.9227   Bremsverzögerung
c_A_max_LDG_vec = (2.2:0.1:2.8).'; % Spaltenvektor mit Werten für c_A_max_LDG [9x1]

landeanvorderung.c_A_max_LDG = c_A_max_LDG_vec(3,1); % c_A_max_LDG 2.3

rho_LDG = ISA.rho(1,1); % rho_B(round(h_50),1); wenn man genau sein will, Unterschied ist aber maginal

%% max Landeflaechenbelastung PS07 formel 24

% betimmung der realen Flaechenbelastung bei der Landung (annahme, dass der Reisetreibstoff aufgebraucht ist)

landeanvorderung.G_LDG_F_real = (G_To / Ergebnisse_stat_Flaechenbelastung.F) * (1-Ergebnis_basis_m.kappa_DP);
landeanvorderung.G_LDG_real = landeanvorderung.G_LDG_F_real * Ergebnisse_stat_Flaechenbelastung.F ; % Gewichtskraft

% maximale Landestrecke

landeanvorderung.G_LDG_F_max_s= (specs.s_LDG_safty - (h_50 ./ landeanvorderung.Eta_LDG))...
    .* (rho_LDG .* landeanvorderung.c_A_max_LDG)./((1/(4*specs.g*landeanvorderung.Eta_LDG) - (1.44 ./landeanvorderung.b_m)));
landeanvorderung.v_s_LDG = sqrt((2 * landeanvorderung.G_LDG_real )/(rho_To * Ergebnisse_stat_Flaechenbelastung.F * landeanvorderung.c_A_max_LDG));
landeanvorderung.v_50 = 1.3 * landeanvorderung.v_s_LDG; % !!! wir duerfen nicht gleich nach dem Start wieder LAnden, sondern muessen fuel dumpen

%% Maximale Landeflaechenbelastung aus v_max_LDG
% Gleichung 26/27 PS7

landeanvorderung.v_max_LDG = landeanvorderung.v_50; % 72:1:84; % in m/s

% Maximale Landeflaechenbelastung bei v_max
landeanvorderung.G_LDG_F_max_v = (rho_LDG./2) .* landeanvorderung.c_A_max_LDG.* (landeanvorderung.v_max_LDG).^2; % Ergebnis 9x31 Matrix, Zeilen unterschiedl. c_A Spalten v_max_LDG

% max Landegewicht
landeanvorderung.G_max_LDG = (rho_LDG./2) .* landeanvorderung.c_A_max_LDG .* (landeanvorderung.v_max_LDG).^2 .* Ergebnisse_stat_Flaechenbelastung.F;
landeanvorderung.m_max_LDG = landeanvorderung.G_max_LDG ./ specs.g;


% Überprüfung der moeglichen Landestrecke nach Formel 28 PS7
% einsetzen von GL26 in GL28 
% s_max_v_LDG =  G_LDG_F_max_v.* (1./(rho_LDG .* c_A_max_LDG)) .* ((1./4.*g.*Eta_LDG(2,1) - 1.44./b_m(1,6))) + (h_50./Eta_LDG(2,1)); 
% einsetzten unf kuertzen (26 in 28)

landeanvorderung.s_max_v_LDG = landeanvorderung.G_LDG_F_max_v .*...
    ((1./(rho_LDG .* landeanvorderung.c_A_max_LDG)) .* ((1./(4 .* specs.g .* landeanvorderung.Eta_LDG) - (1.44 ./ landeanvorderung.b_m))))...
    + (h_50 ./ landeanvorderung.Eta_LDG);








%% Speichern der Daten
save Ergebnisse_Start_Landeanforderungen.mat schub_CR startschub landeanvorderung

