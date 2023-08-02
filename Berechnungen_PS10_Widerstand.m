% Berechnungen für den Widerstand

function [Startwerte_Iteration] = Berechnungen_PS10_Widerstand(Eingabewert_Iteration)



%% Laden der Werte

load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat
load Ergebnisse_Basis_stat_m.mat
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat
load Ergebnisse_Start_Landeanforderungen.mat
load Ergebnisse_Fluegel_Tank_NP.mat;


%% Bestimmung CWi

% phi25_1 = phi_25_I;
% F_1 = F_I; % Fluegelflaeche innen


% phi25_2 = phi_25_A;
% F_2 = F_A; % Fluegelflaeche außen


% test Streckung um Gleitzahl zu ermitteln
% streckung_phi25_max = 19;

Flughoehe = specs.flight_level * 10^2 ;                         % in ft
hoehe = round(unitsratio('m','ft')*Flughoehe);     % in m

p_CR = ISA.p(hoehe); % p_Basis Druckvector aus ISA+15
rho_CR = ISA.rho(hoehe);  % rho_Basis Dichtevector aus ISA+15

G_To = Ergebnisse_stat_Flaechenbelastung.G_To;

%% Berechnung Reynolds

Ind_W.u_cr = specs.Ma_CR * ISA.a(hoehe) ; % fuer reynods Umstroemungsgeachw


v_kin = ISA.kin_visk(hoehe); % kin vis
Ind_W.Re_crit = (Ind_W.u_cr * specs.l_rumpf)/v_kin;
% cf_turb = 0.455 / (log(Re_crit)^2.58); % const





%% Flügel
                                    %d/l Aus Wingdata
Ind_W.CwF_Fluegel_I = 0.0054 * 1 * (1 + 3 * 0.13 * (cos(NP.phi_25_I))^2) * DT.F_I;
Ind_W.CwF_Fluegel_A = 0.0054 * 1 * (1 + 3 * 0.13 * (cos(NP.phi_25_A))^2) * DT.F_A;
%CwF_Fluegel_R = 0.0054 * 1 * (1 + 3 * 0.13 * (cos(phi_25_R))^2) * F_R; %
%Hat keinen Einfluss auf den Cw wert

Ind_W.CwF_Fluegel = (Ind_W.CwF_Fluegel_I + Ind_W.CwF_Fluegel_A)*2;% + CwF_Fluegel_R;

%% Rumpf
Ind_W.b_rumpf = specs.D_rumpf;
Ind_W.h_rumpf = specs.D_rumpf;
Ind_W.r_rumpf_faktor = 1; 

Ind_W.CwF_Rumpf = 0.0031 * Ind_W.r_rumpf_faktor * specs.l_rumpf * (Ind_W.b_rumpf + Ind_W.h_rumpf); % Formel 6 

%% Leitwerk

% Ind_W.CwF_LWT = (Ind_W.CwF_Rumpf + Ind_W.CwF_Fluegel) * 0.24;  % 1.24
Ind_W.CwF_LWT = 1.24; % Annahme aud Aufgabenstellung
%% TW
    % Formfaktor r_n, berücksichtigt Widerstand von Pylonen sowie Interferenz
        %  1,50  alle Triebwerke in Gondeln am Flügel 
        %  1,65  dito, jedoch ein Triebwerk im Rumpfheck
        %  1,25  Triebwerke in Gondeln am Rumpfheck
        %  1,00  interne Triebwerke im Rumpfeinlauf
        %  0,30  Triebwerke in der Flügelwurzel 
    % Formfaktor für Schubumkehr
        %  1,00  mit Schubumkehranlage

% Muss nochmal überprüft werden Ind_W.S_0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ind_W.S_0 = startschub.S0_GTo_To(1,startschub.Startstrecke) * G_To; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ind_W.PSI = 30; % (24s - 40s) [s] %%%%%%% kein plan !!!!!!!!!!! nachfragen

Ind_W.CwF_TW = 1.72 * 1.5 * 1 * ((5 + specs.bypass)/(1+specs.bypass)) * ((Ind_W.S_0)/(ISA.p0 * Ind_W.PSI));

%% Fahrwerk
Ind_W.r_uc = 1; % Formfaktor Fahrwerk 
                %  1,00  voll einziehbares Fahrwerk ohne Verkleidung
                %  1,03  Hauptfahrwerk in PTL-Triebwerksgondel
                %        eingezogen
                %  1,08  Hauptfahrwerk in Rumpfgondel eingezogen
                %        (Transall)
                %  1,25  nicht einziehbares Fahrwerk
                %  1,35  nicht einziehbares Fahrwerk, unverkleidet

%% Cw0 

N_W.r_re = 47 * Ind_W.Re_crit^(-0.2); % const

N_W.Cw0F = N_W.r_re * Ind_W.r_uc * ( Ind_W.CwF_LWT * (Ind_W.CwF_Fluegel + Ind_W.CwF_Rumpf) + Ind_W.CwF_TW);
N_W.Cw0 = N_W.Cw0F / ((DT.F_I + DT.F_A )*2);%+ F_R ); %% Frage die ganze Flügelfläche(F/2) oder nur die Umspülte Flügelfläche 
% Cw0 = 0.016

%% CA Werte


Widerstand.y_CR = 0:0.001:0.85; % gewählter wert aus profilPDF Seite 4     %% REISEFLUG -> CA Anpassen an maxmialen Wert CA_max
oswald_clean = 0.9; % von 0.85 bis 0.9
d_cW_compr_clean = 0;

oswald_comp_LR = 0.8;
d_cW_compr_LR = 0.0005;

oswald_comp_HS = 0.75;
d_cW_compr_HS = 0.002;


for x = 1 : length(Widerstand.y_CR)       
 % CLEAN
    Widerstand.C_w_clean = N_W.Cw0 + ((Widerstand.y_CR(x)^2)/...
        (pi*Ergebnisse_Fluegel.streckung_phi25_max * oswald_clean)) + d_cW_compr_clean; 
    Widerstand.C_w_clean_all(x) = Widerstand.C_w_clean;
 % COMPRSIBILITY LR
    Widerstand.C_w_LR = N_W.Cw0 + ((Widerstand.y_CR(x)^2)/...
        (pi*Ergebnisse_Fluegel.streckung_phi25_max*oswald_comp_LR)) + d_cW_compr_LR;
    Widerstand.C_w_LR_all(x) = Widerstand.C_w_LR;
 % COMPRESIBILITY HIGH SPEED
    Widerstand.C_w_HS = N_W.Cw0 + ((Widerstand.y_CR(x)^2)/...
        (pi*Ergebnisse_Fluegel.streckung_phi25_max*oswald_comp_HS)) + d_cW_compr_HS; 
    Widerstand.C_w_HS_all(x) = Widerstand.C_w_HS;
    
end

%% LANGSAMFLUG / Lsndeanflug / Takeoff

%for x = 1 : length(y)       %CLEAN zum Vergleich
    
%    C_w_clean = Cw0 + ((y(x)^2)/(pi*10.9*0.9))+0;
%    C_w_clean_all(x) = C_w_clean;
%end

%plot(C_w_clean_all,y)
Widerstand.y_to = 0:0.001:startschub.c_A_max_thrust_match;

oswald_comp_TO = 0.85; % 0.8 - 0.85
d_cW_compr_TO = 0.02; % 0.01 bis 0.02

d_cW_compr_GD = 0.0115; % 0.0115 - 0.025

for x = 1 : length(Widerstand.y_to)       %TO CLEAN oC_w_clean_allhne Fahrwerk
    
    Widerstand.C_w_TO_clean = N_W.Cw0 + ((Widerstand.y_to(x)^2)/(pi * Ergebnisse_Fluegel.streckung_phi25_max * oswald_comp_TO)) + d_cW_compr_TO;
    Widerstand.C_w_TO_clean_all(x) = Widerstand.C_w_TO_clean;
    %TO mit Fahrwerk
    Widerstand.C_w_TO = N_W.Cw0 + ((Widerstand.y_to(x)^2)/(pi * Ergebnisse_Fluegel.streckung_phi25_max * oswald_comp_TO)) + d_cW_compr_TO + d_cW_compr_GD;
    Widerstand.C_w_TO_all(x) = Widerstand.C_w_TO;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% LANGSAMFLUG LANDING

Widerstand.y = 0:0.001:landeanvorderung.c_A_max_LDG;
oswald_comp_LDG = 0.8; % 0.75 bis 0.8
d_cW_compr_LDG = 0.065; % 0.055 - 0.065

for x = 1 : length(Widerstand.y)       %LDG ohne Fahrwerk
    
    Widerstand.C_w_LDG_clean = N_W.Cw0 + ((Widerstand.y(x)^2)/(pi * Ergebnisse_Fluegel.streckung_phi25_max * oswald_comp_LDG)) + d_cW_compr_LDG;
    Widerstand.C_w_LDG_clean_all(x) = Widerstand.C_w_LDG_clean;
    %LDG mit Fahrwerk
    Widerstand.C_w_LDG = N_W.Cw0 + ((Widerstand.y(x)^2)/(pi * Ergebnisse_Fluegel.streckung_phi25_max * oswald_comp_LDG)) + d_cW_compr_LDG + d_cW_compr_GD;
    Widerstand.C_w_LDG_all(x) = Widerstand.C_w_LDG;
end

%%%%%%%%%%%%%%%% REZIPROKE GLEITZAHLEN! 


    % E = CA / CW
    GZ.CA_CW_Clean = Widerstand.y_CR./Widerstand.C_w_clean_all;
    GZ.CA_CW_LR = Widerstand.y_CR./Widerstand.C_w_LR_all;
    GZ.CA_CW_HS = Widerstand.y_CR./Widerstand.C_w_HS_all;
    GZ.CA_CW_TO_Clean = Widerstand.y_to./Widerstand.C_w_TO_clean_all;
    GZ.CA_CW_TO = Widerstand.y_to./Widerstand.C_w_TO_all;
    GZ.CA_CW_LDG_Clean = Widerstand.y./Widerstand.C_w_LDG_clean_all;
    GZ.CA_CW_LDG = Widerstand.y./Widerstand.C_w_LDG_all;


if Eingabewert_Iteration == 0
    %% Definition der startwerte für Iteration nach E
    
    Startwerte_Iteration.CA_CW_LR = GZ.CA_CW_LR(1, round(Ergebnisse_stat_Flaechenbelastung.C_A_CR * 10^3));
    Startwerte_Iteration.CA_CW_TO = GZ.CA_CW_TO(1, round(startschub.c_A_max_thrust_match * 10^3));
    Startwerte_Iteration.CA_CW_LDG = GZ.CA_CW_LDG(1,round(landeanvorderung.c_A_max_LDG * 10^3)) ;
    Startwerte_Iteration.CA_CW_Clean = GZ.CA_CW_Clean(1, round(Ergebnisse_stat_Flaechenbelastung.C_A_CR * 10^3));


elseif Eingabewert_Iteration == 1 
    load Ergebnisse_Widerstand_FE2.mat;
    load Ergebnisse_Hochauftrieb_2.mat;


    Startwerte_Iteration.CA_CW_LR = 1 / Ergebnisse_Widerstand_FE2.cW_cA_off_D;
    Startwerte_Iteration.CA_CW_TO = (HA2.CA_max_TO/ HA2.CW_max_TO);
    Startwerte_Iteration.CA_CW_LDG =(HA2.CA_max_ldg_fw/ HA2.CW_max_ldg_fw);
    Startwerte_Iteration.CA_CW_Clean = 1/ Ergebnisse_Widerstand_FE2.cW_cA_off_D;


end

save Ergebnisse_Widerstand.mat GZ Widerstand N_W Ind_W Startwerte_Iteration





