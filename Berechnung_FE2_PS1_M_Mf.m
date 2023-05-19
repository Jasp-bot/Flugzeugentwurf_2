%% PS 1 FE2 Fuel Fraction und Erweiterte MAssenabschaetzung nach Toerenbeck

% function Berechnung_FE2_PS1_M_Mf;
clc
clear all
close all



%% Daten aus Grafik Torenbeek S.454 digitalisiert
% Data_High_lift.daten_specific_flap_weight = readtable('Specific_flap_weight_torenbeek.csv'); % Fig. c-2
% Data_High_lift.daten_secific_weight_leading_edge_HL_div = readtable('Leading_edge_HL_Div_torenbeek.csv'); % Fig. c-3
% 
% % Erste y-Staplte = M_To; 
% % Zweite Spalte = specific weight in [kg/(m^2)]
% Data_High_lift.I_slot_splitflap = table2array(Data_High_lift.daten_specific_flap_weight(1:42, 1:2));
% Data_High_lift.II_2Slot = table2array(Data_High_lift.daten_specific_flap_weight(1:45, 3:4));
% Data_High_lift.III_Flower = table2array(Data_High_lift.daten_specific_flap_weight(1:36, 5:6));
% Data_High_lift.IV_2Slot_Flower = table2array(Data_High_lift.daten_specific_flap_weight(1:31, 7:8));
% Data_High_lift.V_3Slot_Flower = table2array(Data_High_lift.daten_specific_flap_weight(1:31, 9:10));
% 
% Data_High_lift.Hinged_Nose = table2array(Data_High_lift.daten_secific_weight_leading_edge_HL_div(1:18, 11:12));
% Data_High_lift.Slats = table2array(Data_High_lift.daten_secific_weight_leading_edge_HL_div(1:27, 13:14));
% 
% % hold on %% test um InterX zu validieren
% % plot(Data_High_lift.I_slot_splitflap(:,1), Data_High_lift.I_slot_splitflap(:,2));
% % plot([Ergebnis_basis_m.m_To, Ergebnis_basis_m.m_To], [0, 50])
% 
% % Berechnugn der spezifischen Massen der Hochauftriebshilfen
% 
% % Fig. c-2
% % Spec_Weight_I_slot_splitflap 
% Data_High_lift.SP_I_slot_splitflap = InterX([Data_High_lift.I_slot_splitflap(:,1).'; Data_High_lift.I_slot_splitflap(:,2).'], [[Ergebnis_basis_m.m_To, Ergebnis_basis_m.m_To]; [0, 100]]);
% Data_High_lift.Spec_Weight_I_slot_splitflap = Data_High_lift.SP_I_slot_splitflap(2,1);
% 
% % Fig. c-2
% % Spec_Weight_SP_II_2Slot 
% Data_High_lift.SP_II_2Slot = InterX([Data_High_lift.II_2Slot(:,1).'; Data_High_lift.II_2Slot(:,2).'], [[Ergebnis_basis_m.m_To, Ergebnis_basis_m.m_To]; [0, 100]]);
% Data_High_lift.Spec_Weight_SP_II_2Slot = Data_High_lift.SP_II_2Slot(2,1);
% 
% % Fig. c-2
% % Spec_Weight_III_Flower
% Data_High_lift.SP_III_Flower = InterX([Data_High_lift.III_Flower(:,1).'; Data_High_lift.III_Flower(:,2).'], [[Ergebnis_basis_m.m_To, Ergebnis_basis_m.m_To]; [0, 100]]);
% Data_High_lift.Spec_Weight_III_Flower = Data_High_lift.SP_III_Flower(2,1);
% 
% % Fig. c-2
% % Spec_Weight_IV_2Slot_Flower
% Data_High_lift.SP_IV_2Slot_Flower = InterX([Data_High_lift.IV_2Slot_Flower(:,1).'; Data_High_lift.IV_2Slot_Flower(:,2).'], [[Ergebnis_basis_m.m_To, Ergebnis_basis_m.m_To]; [0, 100]]);
% Data_High_lift.Spec_Weight_IV_2Slot_Flower = Data_High_lift.SP_IV_2Slot_Flower(2,1);
% 
% % Fig. c-2
% % Spec_Weight_V_3Slot_Flower
% Data_High_lift.SP_V_3Slot_Flower = InterX([Data_High_lift.V_3Slot_Flower(:,1).'; Data_High_lift.V_3Slot_Flower(:,2).'], [[Ergebnis_basis_m.m_To, Ergebnis_basis_m.m_To]; [0, 100]]);
% Data_High_lift.Spec_Weight_V_3Slot_Flower = Data_High_lift.SP_V_3Slot_Flower(2,1);
% 
% % Fig. c-3
% % Spec_Weight_Hinged_Nose
% %SP_Hinged_Nose = InterX([Data_High_lift.Hinged_Nose(:,1).'; Data_High_lift.Hinged_Nose(:,2).'], [[Ergebnis_basis_m.m_To, Ergebnis_basis_m.m_To]; [0, 100]]);
% Data_High_lift.Spec_Weight_Hinged_Nose = Data_High_lift.Hinged_Nose(18, 2); % Die Grafik hoert bei dem Punkt auf und sieht so aus, als würde der verlauf konstant sein
% 
% % Fig. c-3
% % Spec_Weight_Slats
% Data_High_lift.SP_Slats = InterX([Data_High_lift.Slats(:,1).'; Data_High_lift.Slats(:,2).'], [[Ergebnis_basis_m.m_To, Ergebnis_basis_m.m_To]; [0, 100]]);
% Data_High_lift.Spec_Weight_Slats = Data_High_lift.SP_Slats(2,1);
% 
% save Data_PS1_High_lift_divice_Torenbeek.mat Data_High_lift


%% Loads

load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Widerstand.mat;
load Ergebnisse_Start_Landeanforderungen.mat;
load Ergebnisse_Endwerte_Iteration_V1.mat;
load Ergebnisse_Basis_stat_m.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Leitwerke.mat;
load Data_PS1_High_lift_divice_Torenbeek.mat;   % Laden der oben erstellten daten zur schnelleren benutzung

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Berechnungen Fuel Fraction %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Flughoehe_CR = specs.flight_level * 10^2 ;     % in ft

hoehe_CR = round(distdim(Flughoehe_CR ,'ft','m'));     % in m
hoehe_CL = round(distdim(Flughoehe_CR*(2/3) ,'ft','m'));     % in m
hoehe_CL_ALT = round(distdim(specs.flight_level_ALT * 10^2*(2/3) ,'ft','m'));     % in m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ?????
hoehe_ALT = round(distdim(specs.flight_level_ALT* 10^2,'ft','m'));     % in m
hoehe_HLD = round(distdim(1500, 'ft', 'm'));

%% Climb
% Berechnungen für FuelFraction Methde nach Roskamp

% Formel 9 PS1
FF.S_S0_KF_CL = 0.9 * (ISA.rho(hoehe_CL)/ISA.rho_0) * ...
    exp(-0.35 * specs.Ma_CR*(2/3) * (ISA.p(hoehe_CL)/ISA.p0) * sqrt(specs.bypass));

% Fromel 11 PS1
FF.S_S0_CL = FF.S_S0_KF_CL * schub_CR.S_S0_E;

% Formel 12 PS1
FF.S_G_CL = FF.S_S0_CL * schub_CR.S0_GTo_CR *0.99;          %%%%%%%%%% schub_CR.S0_GTo_CR !!!!!!!!!!!!!!!!!! nicht sicher

FF.eps_CL = 1.1 * (1/(Endwerte_Iteration.CA_CW_Clean));
FF.v_CL = (2/3) * specs.Ma_CR * ISA.a(hoehe_CL);

% Formel 8 PS1
FF.t_CL = hoehe_CR / ((FF.S_G_CL - FF.eps_CL) * FF.v_CL);

% SFC bestimmen für CL
[FF.sfc_CL_daNh, FF.sfc_CL_1PERh, FF.sfc_CL_1PERs] = SFC(hoehe_CL, ((2/3) * specs.Ma_CR), specs.bypass);

% Formel 7 PS1 (umgestellte Formel 6 PS1 um mf3 zu erhalten) 
FF.mf3 = exp(-FF.t_CL * FF.sfc_CL_1PERs * FF.S_G_CL);

%% Cruise

FF.R_CR = specs.max_range_basis_km * 1000 - (FF.v_CL * FF.t_CL);
FF.v_CR = specs.Ma_CR * ISA.a(hoehe_CR);
[FF.sfc_CR_daNh, FF.sfc_CR_1PERh, FF.sfc_CR_1PERs] = SFC(hoehe_CR, (specs.Ma_CR), specs.bypass);

FF.mf4 = 1/exp((FF.R_CR * FF.sfc_CR_1PERs * (1/(Endwerte_Iteration.CA_CW_Clean))/(FF.v_CR)));

FF.mf5 = 1;
%% Diversion

% Climb

FF.S_S0_KF_CL_ALT = 0.9 * (ISA.rho(hoehe_CL_ALT)/ISA.rho_0) * ...
    exp(-0.35 * specs.Ma_CR*(2/3) * (ISA.p(hoehe_CL_ALT)/ISA.p0) * sqrt(specs.bypass));     %%%%%%%%%%% Nicht sicher ob richtige Machzahl

FF.S_S0_CL_ALT = FF.S_S0_KF_CL_ALT * schub_CR.S_S0_E;

FF.S_G_CL_ALT = FF.S_S0_CL_ALT * schub_CR.S0_GTo_CR *0.99;          %%%%%%%%%% schub_CR.S0_GTo_CR !!!!!!!!!!!!!!!!!! nicht sicher


FF.Ma_CL_ALT = ((2/3) * specs.v_HLD) / ISA.a(hoehe_CL_ALT);
FF.v_CL_ALT = ((2/3) * specs.v_HLD);

FF.t_CL_ALT = hoehe_ALT / ((FF.S_G_CL_ALT - ((FF.eps_CL))) * FF.v_CL_ALT);  %% nicht sicher ob mit v_CL rechnen oder neuer geschwindigkeit
 

[FF.sfc_CL_ALT_daNh, FF.sfc_CL_ALT_1PERh, FF.sfc_CL_ALT_1PERs] = SFC(hoehe_CL_ALT, FF.Ma_CL_ALT, specs.bypass);

FF.mf6 = exp(-FF.t_CL_ALT * FF.sfc_CL_ALT_1PERs * FF.S_G_CL_ALT);

% Diversion CR

FF.R_ALT = specs.R_ALT ;% - (FF.v_CL_ALT * FF.t_CL_ALT); Es werden mit 200 nm im CR DIV gerechnet
FF.v_ALT = specs.v_HLD;
FF.Ma_CR_ALT = FF.v_ALT / ISA.a(hoehe_ALT);

[FF.sfc_CR_ALT_daNh, FF.sfc_CR_ALT_1PERh, FF.sfc_CR_ALT_1PERs] = SFC(hoehe_ALT, FF.Ma_CR_ALT, specs.bypass);

FF.mf7 = 1/exp((FF.R_ALT * FF.sfc_CR_ALT_1PERs * (1/(Endwerte_Iteration.CA_CW_Clean))/(FF.v_ALT)));

FF.mf8 = 1;

%% Holding

FF.Ma_HLD = specs.v_HLD / ISA.a(hoehe_HLD);

[FF.sfc_HLD_daNh, FF.sfc_HLD_1PERh, FF.sfc_HLD_1PERs] = SFC(hoehe_HLD, FF.Ma_HLD, specs.bypass);

FF.mf9 = 1/(exp(specs.t_HLD * FF.sfc_HLD_1PERs * (1/Endwerte_Iteration.CA_CW_Clean)));


%% Fuel Fraction ges


% Nach Roskam
FF.mf0 = 0.992; 
FF.mf1 = 0.996;
FF.mf2 = 0.996;
FF.mf10 = 1; % hat kristof gesagt

%%%%%%%%%%%%%%%%%%%%%%%%  mfi muss noch mal überprüft werden, ich habe
%%%%%%%%%%%%%%%%%%%%%%%%  jetzt mit alles Massenanteilen gerechnet, von 0
%%%%%%%%%%%%%%%%%%%%%%%%  bis 10 nicht von 2 bis 10. Muss eventuell nochmal
%%%%%%%%%%%%%%%%%%%%%%%%  angepasst werden
% Matrix für Produkt auf mf0 bis mf10
FF.mfi = [FF.mf0, FF.mf1, FF.mf2, FF.mf3,...
    FF.mf4, FF.mf5, FF.mf6, FF.mf7,...
    FF.mf8, FF.mf9, FF.mf10]; 

FF.Mff = prod(FF.mfi(3:11)); % Faktorprodukt  
FF.mf_oC = (1-FF.Mff) * Ergebnis_basis_m.m_To; %%%%%%%%%%%%%%%%%%%%% Muss verändert werden für iteration

FF.Mff_2_10 = prod(FF.mfi(3:11)); % Gesamtanteil Mff ab TO
FF.Mff_2_5 = prod(FF.mfi(3:6)); % Nur Reisefluganteil ohne DIV

% Kraftstoffmassenfaktor neu
FF.Kappa_ges = 1- FF.Mff_2_10 + 0.05 * (1 - FF.Mff_2_5); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Massenabschaetzung nach Torenbeek
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rumpfmasse mit EAS

% bestimmung der umspuelten Rumpfflaeche Sg_12
% Achtung Approximation der Flaeche, keine genaue Angabe!!!!!!!!!!
% Veraenderliche Werte je nach Geometrie
NR_Rumpf_Geometrie.l_mittelteil = 52; % [m]
NR_Rumpf_Geometrie.l_front = 8.08; % [m]
NR_Rumpf_Geometrie.r1_front = 0.33; % [m]
NR_Rumpf_Geometrie.h_spitze = 0.094; 
NR_Rumpf_Geometrie.l_coan = 15.01; 
NR_Rumpf_Geometrie.r1_coan = 0.96;

% Nebrnrechnungen
NR_Rumpf_Geometrie.s_front = sqrt((specs.R_rumpf - NR_Rumpf_Geometrie.r1_front)^2 + NR_Rumpf_Geometrie.l_front^2);
NR_Rumpf_Geometrie.s_coan = sqrt((specs.R_rumpf - NR_Rumpf_Geometrie.r1_coan)^2 + NR_Rumpf_Geometrie.l_coan^2);


% Rechnungen der Einzelflaechen
NR_Rumpf_Geometrie.Mantel_rumpf = 2 * pi * specs.R_rumpf * NR_Rumpf_Geometrie.l_mittelteil;
NR_Rumpf_Geometrie.Mantel_front = pi * NR_Rumpf_Geometrie.s_front * (specs.R_rumpf + NR_Rumpf_Geometrie.r1_front);
NR_Rumpf_Geometrie.Kugelkappe = pi * (NR_Rumpf_Geometrie.h_spitze^2 + NR_Rumpf_Geometrie.r1_front);
NR_Rumpf_Geometrie.Mantel_coan = pi * NR_Rumpf_Geometrie.s_coan * (specs.R_rumpf + NR_Rumpf_Geometrie.r1_coan);
NR_Rumpf_Geometrie.Endkappe = pi * NR_Rumpf_Geometrie.r1_coan^2;
% Gesamtflaeche
M_Rumpf.Sg_12 = NR_Rumpf_Geometrie.Mantel_rumpf + NR_Rumpf_Geometrie.Mantel_front + NR_Rumpf_Geometrie.Kugelkappe + NR_Rumpf_Geometrie.Mantel_coan + NR_Rumpf_Geometrie.Endkappe;

% Bestimmung v_D_EAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Müssen wir wirklich mit 20000ft rechen der
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mit 34000ft (v_D_2)
Ma_D = specs.Ma_CR + 0.05;
v_D_1_TAS = Ma_D * ISA.a(hoehe_ALT);
v_D_2_TAS = (specs.Ma_CR * ISA.a(hoehe_ALT) + convvel(60, 'kts','m/s')) ; 
v_D_2_TAS_kts = convvel(v_D_2_TAS,'m/s', 'kts') ;
v_D_1_TAS_kts = convvel(v_D_1_TAS,'m/s', 'kts') ;
M_Rumpf.v_D_EAS = v_D_1_TAS * sqrt(ISA.rho(hoehe_ALT)/ISA.rho_0) ;        %%%%%%%%%%%%%%%%%%%%%%%%% Nicht sicher !!!!!!!!!!!!!!!!!!!!!!!!!




% Abstand der phi25 linien von Fluegel und HLW an der Wurzel l_t
NR_l_t.versatz_NP_HLW_VK = HLW.kabinenversatz - HLW.x_NP_ges/2;
NR_l_t.versatz_phi25_HLW_VK = NR_l_t.versatz_NP_HLW_VK - (HLW.l_i_R * 0.25);

NR_l_t.versatz_NP_DTF_VK = NP.versatz_HK + DT.l_i_R - NP.x_NP_ges;
NR_l_t.versatz_phi25_DTF_VK = NR_l_t.versatz_NP_DTF_VK - (DT.l_i_R * 0.25);
% der Leitwerkshebelarm zwischen Fluegel und HLW wurde in FE1 in der PS9
% berechnet und ist unter HLE.r zu finden.
M_Rumpf.l_t = NR_l_t.versatz_phi25_DTF_VK + HLW.r - NR_l_t.versatz_phi25_HLW_VK;

% Berechnung von W_Rumpf



%%%%%%%%%%%%% ??????? sehr niedrig
% W_Rumpf = Rumpfmasse !!!!!!!!!!!!!!!!!!!!!
M_Rumpf.W_Rumpf = 0.23 * sqrt(M_Rumpf.v_D_EAS * (M_Rumpf.l_t)/(specs.D_rumpf + specs.h_rumpf)) * (M_Rumpf.Sg_12)^1.2; 





%% Fluegelmasse ueber Bruchlastvielfaches

NR_M_Fluegel.const = 4.58 * 10^(-3);
NR_M_Fluegel.b_ref = 1.905; % [m]
NR_M_Fluegel.b_s = sqrt((Ergebnisse_Fluegel.b/2)^2 + ( (NP.versatz_HK-(DT.l_a/2)) + (DT.l_i_R/2) )^2);
NR_M_Fluegel.k_no = 1 + sqrt((NR_M_Fluegel.b_ref)/(NR_M_Fluegel.b_s));

NR_M_Fluegel.k = (1 + Ergebnisse_Fluegel.lambda)^(0.4);

NR_M_Fluegel.k_e = 0.95; % two wing mounted engines

NR_M_Fluegel.k_uc = 1; % for wing mounted undercarriage or 0.95 for not Wingmounted undercarriage %%%%%%%% Nicht sicher !!!!!

NR_M_Fluegel.Lambda_50 = tan((DT.s_A)/(NP.versatz_HK + DT.l_i_A/2)); % Annahme Pfeilung von Außenfluegel


NR_M_Fluegel.W_des = Ergebnis_basis_m.m_OE + Ergebnis_basis_m.m_payload;
NR_M_Fluegel.k_st = 1 + (9.06*10^(-4)) * (((Ergebnisse_Fluegel.b * cos(Ergebnisse_Fluegel.phi_VK_max))^3)/(NR_M_Fluegel.W_des)) *...
    (((M_Rumpf.v_D_EAS/100)/(0.13))^2) * cos(NR_M_Fluegel.Lambda_50); % Annahme (t/c)_r = 0.13  W_des = unbekannt

NR_M_Fluegel.k_b = 1; % for catilever wings otherwise k_b = 1 - nue_s^2 || neu_s ditance trut to wing root

NR_M_Fluegel.d_l_root = 0.13;

M_Fluegel.n_max = 2.1 + (10900)/(4540 + Ergebnis_basis_m.m_To); % Achtung hier steht m_To drin, muss füt Iteration veraendert werden

M_Fluegel.n_ult = M_Fluegel.n_max * 1.5;

NR_M_Fluegel.W_F_initial = Ergebnis_basis_m.m_OE * 0.3; %% Initiale Annahme, dass Winggroup weight ca 30% M_OE sind


% Berechnung W_F_basic unter initialen Annahmen
M_Fluegel.W_F_basic = NR_M_Fluegel.const * NR_M_Fluegel.k_no * NR_M_Fluegel.k *...
    NR_M_Fluegel.k_e * NR_M_Fluegel.k_uc * NR_M_Fluegel.k_st *...
    (NR_M_Fluegel.k_b * M_Fluegel.n_ult * (NR_M_Fluegel.W_des - 0.8* NR_M_Fluegel.W_F_initial))^(0.55) * ...
    Ergebnisse_Fluegel.b^(1.675) * NR_M_Fluegel.d_l_root^(-0.45) * cos(NR_M_Fluegel.Lambda_50)^(-1.325); % Annahme d_l_root = (t/c)_r = 0.13






pauschaler_Wert_VK = 0.9; % Wert angenommen als flaeche, die nicht in klappen/ slpats umgewandet werden kann
pauschaler_Wert_HK = 0.15;
ende_querruder = 0.3;


% Rechnung basiert auf Visualisierung Fluegel Plot Tank / Holme
area1 = trapz([0, Ergebnisse_Fluegel.b/2-specs.R_rumpf], [DT.l_a, NP.versatz_HK+DT.l_a_R]); % Flaeche unter gesamten Vorderkante
area2 = trapz([0, Ergebnisse_Fluegel.b/2-specs.R_rumpf], [DT.l_a * (1-0.15), NP.versatz_HK+DT.l_a_R * (1-0.15)]); % Flaeche unter erstem Holm
area_VK_lift_div = (area1 - area2) * (pauschaler_Wert_VK); % Fläche für Vorderkanten High lift divices

area3 = trapz([(Ergebnisse_Fluegel.b/2)*ende_querruder, DT.s_A],...
    [(1-0.65)*Ergebnisse_Fluegel.Fluegeltiefen_eta(1,(100-ende_querruder*100))+(tan(DT.phi_HK_dt)*ende_querruder*(Ergebnisse_Fluegel.b/2)),...
    NP.versatz_HK+DT.l_i_A * (1- 0.65)]); % beginn bei 30% halbspannweite bis kink flaeche unter Holm
area4 = trapz([DT.s_A, DT.s_I + DT.s_A],...
    [NP.versatz_HK+DT.l_i_A * (1- 0.65), NP.versatz_HK+DT.l_i_I * (1- 0.65)]); % Innenflaeche unter Holm

area5 = trapz([(Ergebnisse_Fluegel.b/2)*ende_querruder, DT.s_A],...
    [(tan(DT.phi_HK_dt)*ende_querruder*(Ergebnisse_Fluegel.b/2)), NP.versatz_HK]); % beginn bei 30% halbspannweite bis kink flaeche unter außenfluegel(Hinterkante)

area6 = DT.s_I * NP.versatz_HK;

area_HK_lift_div = ((area3 + area4) - (area5 + area6)) * (1-pauschaler_Wert_HK);


% % masse trailing edge flap und leading edge flap
NR_M_Fluegel.S_f_tef = area_HK_lift_div;
NR_M_Fluegel.S_f_lef = area_VK_lift_div;
NR_M_Fluegel.W_tef = Data_High_lift.Spec_Weight_I_slot_splitflap * NR_M_Fluegel.S_f_tef;
NR_M_Fluegel.W_lef = Data_High_lift.Spec_Weight_Slats * NR_M_Fluegel.S_f_lef;

% Formel Torenbeek S 454 c-8
M_Fluegel.W_high_lift_div = NR_M_Fluegel.W_tef + NR_M_Fluegel.W_lef;
% Formel Torenbeek S 454
M_Fluegel.W_SP = 0.015 * NR_M_Fluegel.W_F_initial; % oder 12.2[kg/m^2]
% 
% 
% Initiales Gewicht der Fluegelgruppe W_W im Torenbeek oder W_F
% Formel Torenbeek S 455 c-11
W_F = M_Fluegel.W_F_basic + 1.2 * (M_Fluegel.W_high_lift_div + M_Fluegel.W_SP); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ergebnis sehr gering
 


% W undercarriage Torenbeek S282 A,B,C,D S283
k_uc = 1;
A = [18.1; 9.1]; % X(1,1) = main Gear, X(2,1) = nose gear beide retractable
B = [0.131; 0.082];
C = [0.019; 0];
D = [2.23*10^(-5); 2.97*10^(-6)];
% Formel (8-17) S.282
W_uc_teilergebnis = (k_uc .* (A + B .* Ergebnis_basis_m.m_To.^(3/4) + C .* Ergebnis_basis_m.m_To + D .* Ergebnis_basis_m.m_To.^(3/2)));
W_uc = sum(W_uc_teilergebnis); %%%%%%%%%%%%%%%%%%%% Ergebnis Plausibel



% W_sc Surface control group S.283

k_sc = [0.23; 0.44; 0.64;]; % [light AC ohne duplicated sys.controls; Transport-AC/Trainer manually controlled; Transportairpanes with powercontrol]
% Formel (8-18) S.283
W_sc_initial = (k_sc(3,1) .* Ergebnis_basis_m.m_To.^(2/3)) .* 0.768; % 0.768 sind fuer umrechnungsfaktor zu kg
% +20% masse fuer leading edge flap +15% masse lift tumper controls
m_cockpitcontrols = 50; % [kg] Formel (8-19)
m_autopilot = 9 * Ergebnis_basis_m.m_To.^(1/5); % [kg] Formel (8-20)

% tabelle 8-7 S 284
m_maneuver_control = 0.773 * Ergebnis_basis_m.m_To^(0.60);
deflaction_angle_flaps = deg2rad(8); %% Wert einfach angenommen keine Ahnung ob richtig
% einer der beiden werte verwenden m_trailing_edge_zylinder_actuation oder m_trailing_edge_fouler_flap
m_trailing_edge_zylinder_actuation = 5.569 * (area_HK_lift_div * sin(deflaction_angle_flaps))^0.92; %% Achtung nur eine der beiden für rechnung verwenden
m_trailing_edge_fouler_flap = 11.02 * (area_HK_lift_div * sin(deflaction_angle_flaps))^0.92; %% Achtung nur eine der beiden für rechnung verwenden
m_trailing_edge = [m_trailing_edge_zylinder_actuation; m_trailing_edge_fouler_flap];
m_slat_control = 11.23 * (area_VK_lift_div)^0.82;
% Masse horizontal tailplane acuator
k_hc = 0.44; % dual actuator otherwise (0.31 single actuator)
S_he = HLW.F - (2 * HLW.F_R) ;    % Flaeche Horizontal tail plane (exposed)
V_max = ISA.a(hoehe_CR) * specs.Ma_MO; % nicht sicher on Ma_MO oder Ma_CR ????
deflaction_Horizontal_tailplane = deg2rad(35);
m_HLW_control = k_hc * (S_he * (V_max^0.5) * sin(deflaction_Horizontal_tailplane))^(0.88); %%%%%%% sehr unsicher
% speedbrake control
s_speed_brake = 20; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Annahme
m_speed_brake = 40.4 * s_speed_brake^0.92;
% Lift dumper controls
deflaction_lift_dampner = deg2rad(45);
s_lift_dampner = 20; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Annahme
m_lift_dampner = 20.0 * (s_lift_dampner * sin(deflaction_lift_dampner))^0.92;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Achtung hier Werte
% %aus Tabelle 8-7 S 284 ueberpruefen !!!!!!!!!!!!!!!!!!!
W_sc = W_sc_initial + W_sc_initial * 0.2 + W_sc_initial *0.15 + m_cockpitcontrols + m_autopilot +...
    m_maneuver_control + m_trailing_edge(1,1) + m_slat_control + m_HLW_control + m_speed_brake + m_lift_dampner; %%%%%%%%% Ergebnis in einem plausibelem bereich kein plan ob richtig (etwas hoch im vergleich 747-400)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Engine section nacelle group





