%% PS 1 2 Fuel Fraction und Erweiterte MAssenabschaetzung nach Toerenbeck

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

% Laden der oben erstellten daten zur schnelleren benutzung muss angepasst 
% werden wenn genaueres zu den Hochauftriebshilfen bekannt ist 
load Data_PS1_High_lift_divice_Torenbeek.mat;   



%% Anfangswerte festlegen

M_TO_initial = Ergebnis_basis_m.m_To;
M_OE_initial = Ergebnis_basis_m.m_OE;
M_del_empty_initial = 140000; % Annahme die sehr nahe am ersten Schleifenduchlauf liegt
M_Zero_Fuel_initial = Ergebnis_basis_m.m_OE + specs.m_pax_all_eco + specs.m_cargo; % Annahme für all Eco Version nicht ACHTUNG
% Kappa_initial = 0.3064; % initiale annahme aus erstem Schleifenduchlauf
Zaehlvariabele = 0;
delta_M_to = 100;       % Anfangswert

% Unter der annahme, dass 50% der Strukturmasse aus CFK gefertigt werden und CFK 40% leichter ist als ALU
Technologiefaktor_ALU_CFK = 0.5 * 0.4;



while abs(delta_M_to) > 0.0001
    
   
    
    %% Berechnungen Fuel Fraction %%
    
    
    Flughoehe_CR = specs.flight_level * 10^2 ;     % in ft
    
    hoehe_CR = round(unitsratio('m','ft')*Flughoehe_CR);     % in m
    hoehe_CL = round(unitsratio('m','ft')*Flughoehe_CR*(2/3) );     % in m
    hoehe_CL_ALT = round(unitsratio('m','ft')*specs.flight_level_ALT * 10^2*(2/3) );     % in m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ?????
    hoehe_ALT = round(unitsratio('m','ft')*specs.flight_level_ALT* 10^2);     % in m
    hoehe_HLD = round(unitsratio('m', 'ft')*1500);
    
    %% Climbn_ult
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
    FF.mf_oC = (1-FF.Mff) * M_TO_initial; %%%%%%%%%%%%%%%%%%%%% Muss verändert werden für iteration
    


    FF.Mff_2_10 = prod(FF.mfi(3:11)); % Gesamtanteil Mff ab TO
    FF.Mff_2_5 = prod(FF.mfi(3:6)); % Nur Reisefluganteil ohne DIV
    FF.Mff_6_10 = prod(FF.mfi(7:11));

    % Kraftstoffmassenfaktor neu
    FF.Kappa_ges = 1- FF.Mff_2_10 + 0.05 * (1 - FF.Mff_2_5); 
    
    M_take_off_initial.M_fuel = M_TO_initial * FF.Kappa_ges;
    
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
    M_Rumpf.v_D_EAS = v_D_1_TAS * sqrt(ISA.rho(hoehe_ALT)/ISA.rho_0) ;       
%     M_Rumpf.v_D_TAS = v_D_1_TAS;
    
    
    
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
    M_Rumpf.W_Rumpf = (0.23 * sqrt(M_Rumpf.v_D_EAS * (M_Rumpf.l_t)/(specs.D_rumpf + specs.h_rumpf)) * (M_Rumpf.Sg_12)^1.2) * (1- Technologiefaktor_ALU_CFK) ; 
    
    
    M_Airframe_Structur.Fuselage_Group = M_Rumpf.W_Rumpf;
    
    
    %% Fluegelmasse ueber Bruchlastvielfaches
    
    NR_M_Fluegel.const = 4.58 * 10^(-3);
    NR_M_Fluegel.b_ref = 1.905; % [m]
    NR_M_Fluegel.b_s = 2 * sqrt((Ergebnisse_Fluegel.b/2)^2 + ( (NP.versatz_HK-(DT.l_a/2)) + (DT.l_i_R/2) )^2);
    NR_M_Fluegel.k_no = 1 + sqrt((NR_M_Fluegel.b_ref)/(NR_M_Fluegel.b_s));
    
    NR_M_Fluegel.k = (1 + Ergebnisse_Fluegel.lambda)^(0.4);
    
    NR_M_Fluegel.k_e = 0.95; % two wing mounted engines
    
    NR_M_Fluegel.k_uc = 1; % for wing mounted undercarriage or 0.95 for not Wingmounted undercarriage %%%%%%%% Nicht sicher !!!!!
    
    NR_M_Fluegel.Lambda_50 = tan((DT.s_A)/(NP.versatz_HK + DT.l_i_A/2)); % Annahme Pfeilung von Außenfluegel
    
    
    NR_M_Fluegel.W_des = M_Zero_Fuel_initial;
    NR_M_Fluegel.k_st = 1 + (9.06*10^(-4)) * (((Ergebnisse_Fluegel.b * cos(Ergebnisse_Fluegel.phi_VK_max))^3)/(NR_M_Fluegel.W_des)) *...
        (((M_Rumpf.v_D_EAS/100)/(0.13))^2) * cos(NR_M_Fluegel.Lambda_50); % Annahme (t/c)_r = 0.13  W_des = unbekannt
    
    NR_M_Fluegel.k_b = 1; % for catilever wings otherwise k_b = 1 - nue_s^2 || neu_s ditance trut to wing root
    
    NR_M_Fluegel.d_l_root = 0.13;
    
    M_Fluegel.n_max = 2.1 + (10900)/(4540 + M_TO_initial); % Achtung hier steht m_To drin, muss füt Iteration veraendert werden
    
    M_Fluegel.n_ult = M_Fluegel.n_max * 1.5;
    
    if Zaehlvariabele == 0;

        NR_M_Fluegel.W_F_initial = M_OE_initial * 0.3; %% Initiale Annahme, dass Winggroup weight ca 30% M_OE sind
    
    else Zaehlvariabele > 0;
        NR_M_Fluegel.W_F_initial = M_Airframe_Structur.Wing_Group;

    end
        
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
    
    % beginn bei 30% halbspannweite bis kink flaeche unter Holm
    area3 = trapz([(Ergebnisse_Fluegel.b/2)*ende_querruder, DT.s_A],...
        [(1-0.65)*Ergebnisse_Fluegel.Fluegeltiefen_eta(1,(100-ende_querruder*100))+(tan(DT.phi_HK_dt)*ende_querruder*(Ergebnisse_Fluegel.b/2)),...
        NP.versatz_HK+DT.l_i_A * (1- 0.65)]); 
    % Innenflaeche unter Holm
    area4 = trapz([DT.s_A, DT.s_I + DT.s_A],...
        [NP.versatz_HK+DT.l_i_A * (1- 0.65), NP.versatz_HK+DT.l_i_I * (1- 0.65)]); 
    % beginn bei 30% halbspannweite bis kink flaeche unter außenfluegel(Hinterkante)
    area5 = trapz([(Ergebnisse_Fluegel.b/2)*ende_querruder, DT.s_A],...
        [(tan(DT.phi_HK_dt)*ende_querruder*(Ergebnisse_Fluegel.b/2)), NP.versatz_HK]); 
    % Fläche unter Fluegel (rechteck zwischen kink und rumpf)
    area6 = DT.s_I * NP.versatz_HK;
    
    area_HK_lift_div = ((area3 + area4) - (area5 + area6)) * (1-pauschaler_Wert_HK);
    
    %%%%%%%%%%%%%%%%%%%
    % Neue Berechnung 
    load Ergebnisse_Hochauftrieb_2.mat
    load Ergebnisse_Start_Landeanforderungen.mat
    %NR_M_Fluegel.W_tef = 

    kf1 = 1.3; % Double Slotted Fowler
    kf2 = 1.25; % Double Slotted flaps with variable Geometry

    k_f = kf1 * kf2;

    
    S_f = F_FOWLER;%(spannweite_flaps * (Ergebnisse_Fluegel.b/2)) * Ergebnisse_Fluegel.l_m;  % Fläche Klappen über mittlere Flügeltiefe
    b_fs = spannweite_flaps * (Ergebnisse_Fluegel.b); % Spannweite Klappen
    V_lf = landeanvorderung.v_50;
    beta = deg2rad(45);%° % Vorgegeben als Optimum PS 6
    pfeilung = (Ergebnisse_Fluegel.phi_25_max); 
    t_c = specs.d_l;
    %rad oder degree?
    formel = 2.706 * k_f * ((S_f * b_fs)^(3/16)) * ((V_lf/100)^2 * ((sin(beta) * cos(pfeilung))/(t_c)))^(3/4) * S_f;

    
    
    %%%%%%%%%%%%%%%%%%%
    
    % % masse trailing edge flap und leading edge flap
    NR_M_Fluegel.S_f_tef = area_HK_lift_div;
    NR_M_Fluegel.S_f_lef = area_VK_lift_div;
    NR_M_Fluegel.W_tef = formel; %Data_High_lift.Spec_Weight_I_slot_splitflap * NR_M_Fluegel.S_f_tef;
    NR_M_Fluegel.W_lef = Data_High_lift.Spec_Weight_Slats * NR_M_Fluegel.S_f_lef;
    
    % Formel Torenbeek S 454 c-8
    M_Fluegel.W_high_lift_div = NR_M_Fluegel.W_tef + NR_M_Fluegel.W_lef;
    % Formel Torenbeek S 454
    M_Fluegel.W_SP = 0.015 * NR_M_Fluegel.W_F_initial; % oder 12.2[kg/m^2]
    % 
    

    
    % Initiales Gewicht der Fluegelgruppe W_W im Torenbeek oder W_F
    % Formel Torenbeek S 455 c-11
    M_Fluegel.W_F = (M_Fluegel.W_F_basic + 1.2 * (M_Fluegel.W_high_lift_div + M_Fluegel.W_SP))...
        * (1 - Technologiefaktor_ALU_CFK); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ergebnis sehr gering M_Fluegel.W_F
    M_Airframe_Structur.Wing_Group = M_Fluegel.W_F;
    




    % Masse TailGroup S.281 annahme aus Text leicht opimistischer (technologiefaktor)
    W_tail = 0.02 * M_TO_initial;
    
    M_Airframe_Structur.Tail_Group = W_tail;
    
%% Höhenleitwerk Massenberechnung
    
    HLW_Mass.const = 4.58 * 10^(-3);
    HLW_Mass.b_ref = 1.905; % [m]
    HLW_Mass.b_s = HLW.l_phi50;
    HLW_Mass.k_no = 1 + sqrt((HLW_Mass.b_ref)/(HLW_Mass.b_s));
    
    HLW_Mass.k = (1 + HLW.phi_VK)^(0.4);
    
    HLW_Mass.k_e = 1; % no engines
    
    HLW_Mass.k_uc = 0.95; % for wing mounted undercarriage or 0.95 for not Wingmounted undercarriage %%%%%%%% Nicht sicher !!!!!
    
    HLW_Mass.Lambda_50 = HLW.phi_50; % Annahme Pfeilung von Außenfluegel
    HLW_Mass.k_st = 1 + 9.06*10^(-4)*(((HLW.b*cos(HLW.phi_VK))^3/NR_M_Fluegel.W_des)*(M_Rumpf.v_D_EAS/(100*specs.HLW_d2l))^2*cos(HLW_Mass.Lambda_50));
  
    
    HLW_Mass.k_b = 1; % for catilever wings otherwise k_b = 1 - nue_s^2 || neu_s ditance trut to wing root
    
    if Zaehlvariabele == 0;

       HLW_Mass.W_F_initial = 3000; %% Initiale Annahme für HLW based on Uebung Schwerpunkt
    
    else Zaehlvariabele > 0;
        HLW_Mass.W_F_initial = M_Airframe_Structur.HLW;

    end
        
    %Berechnung W_HLW_basic
    M_HLW.W_HLW_basic = HLW_Mass.const * HLW_Mass.k_no * HLW_Mass.k *...
        HLW_Mass.k_e * HLW_Mass.k_uc * HLW_Mass.k_st *...
        (HLW_Mass.k_b * M_Fluegel.n_ult * (NR_M_Fluegel.W_des - 0.8* HLW_Mass.W_F_initial))^(0.55) * ...
        HLW.b^(1.675) * specs.HLW_d2l^(-0.45) * cos(HLW_Mass.Lambda_50)^(-1.325);

    M_Airframe_Structur.HLW = M_HLW.W_HLW_basic;

%% Seitenleitwerk Massenberechnung
    
    SLW_Mass.const = 4.58 * 10^(-3);
    SLW_Mass.b_ref = 1.905; % [m]
    SLW_Mass.b_s = SLW.l_phi50;
    SLW_Mass.k_no = 1 + sqrt((SLW_Mass.b_ref)/(SLW_Mass.b_s));
    
    SLW_Mass.k = (1 + SLW.phi_VK)^(0.4);
    
    SLW_Mass.k_e = 1; % no engines
    
    SLW_Mass.k_uc = 0.95; % for wing mounted undercarriage or 0.95 for not Wingmounted undercarriage %%%%%%%% Nicht sicher !!!!!
    
    SLW_Mass.Lambda_50 = SLW.phi_50; % Annahme Pfeilung von Außenfluegel
    SLW_Mass.k_st = 1 + 9.06*10^(-4)*(((SLW.b*cos(SLW.phi_VK))^3/NR_M_Fluegel.W_des)*(M_Rumpf.v_D_EAS/(100*specs.HLW_d2l))^2*cos(SLW_Mass.Lambda_50));
  
    
    SLW_Mass.k_b = 1; % for catilever wings otherwise k_b = 1 - nue_s^2 || neu_s ditance trut to wing root
    
    if Zaehlvariabele == 0;

       SLW_Mass.W_F_initial = 3000; %% Initiale Annahme für SLW based on Uebung Schwerpunkt
    
    else Zaehlvariabele > 0;
        SLW_Mass.W_F_initial = M_Airframe_Structur.SLW;

    end
        
    %Berechnung W_SLW_basic
    M_SLW.W_SLW_basic = SLW_Mass.const * SLW_Mass.k_no * SLW_Mass.k *...
        SLW_Mass.k_e * SLW_Mass.k_uc * SLW_Mass.k_st *...
        (SLW_Mass.k_b * M_Fluegel.n_ult * (NR_M_Fluegel.W_des - 0.8* SLW_Mass.W_F_initial))^(0.55) * ...
        SLW.b^(1.675) * specs.HLW_d2l^(-0.45) * cos(SLW_Mass.Lambda_50)^(-1.325);

    M_Airframe_Structur.SLW = M_SLW.W_SLW_basic;

    



    % W undercarriage Torenbeek S282 A,B,C,D S283
    k_uc = 1;
    A = [18.1; 9.1]; % X(1,1) = main Gear, X(2,1) = nose gear beide retractable
    B = [0.131; 0.082];
    C = [0.019; 0];
    D = [2.23*10^(-5); 2.97*10^(-6)];
    % Formel (8-17) S.282
    W_uc_teilergebnis = (k_uc .* (A + B .* M_TO_initial.^(3/4) + C .* M_TO_initial + D .* M_TO_initial.^(3/2)));
    W_uc = sum(W_uc_teilergebnis); %%%%%%%%%%%%%%%%%%%% Ergebnis Plausibel
    
    M_Airframe_Structur.Undercarriage_Group = W_uc;


    % MAIN GEAR
    % Von https://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_10_Mass.pdf
    W_MainGear = 18.1 + 0.131*M_TO_initial^(3/4) + 0.019*M_TO_initial + 2.23*10^(-5)*M_TO_initial^(3/2);
    
    M_Airframe_Structur.MainGear = W_MainGear;

    % FRONT GEAR
    % Von https://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_10_Mass.pdf
    W_FrontGear = 9.1 + 0.082*M_TO_initial^(3/4)+ 2.97*10^(-6)*M_TO_initial^(3/2);
    
    M_Airframe_Structur.FrontGear = W_FrontGear;
    
    
    % W_sc Surface control group S.283
    
    NR_W_sc.k_sc = [0.23; 0.44; 0.64;]; % [light AC ohne duplicated sys.controls; Transport-AC/Trainer manually controlled; Transportairpanes with powercontrol]
    % Formel (8-18) S.283
    NR_W_sc.W_sc_initial = (NR_W_sc.k_sc(3,1) .* M_TO_initial.^(2/3)) .* 0.768; % 0.768 sind fuer umrechnungsfaktor zu kg
    % +20% masse fuer leading edge flap +15% masse lift tumper controls
    NR_W_sc.m_cockpitcontrols = 50; % [kg] Formel (8-19)
    NR_W_sc.m_autopilot = 9 * M_TO_initial.^(1/5); % [kg] Formel (8-20)
    
    % tabelle 8-7 S 284
    NR_W_sc.m_maneuver_control = 0.773 * M_TO_initial^(0.60);
    deflaction_angle_flaps = deg2rad(8); %% Wert einfach angenommen keine Ahnung ob richtig
    % einer der beiden werte verwenden m_trailing_edge_zylinder_actuation oder m_trailing_edge_fouler_flap
    NR_W_sc.m_trailing_edge_zylinder_actuation = 5.569 * (area_HK_lift_div * sin(deflaction_angle_flaps))^0.92; %% Achtung nur eine der beiden für rechnung verwenden
    NR_W_sc.m_trailing_edge_fouler_flap = 11.02 * (area_HK_lift_div * sin(deflaction_angle_flaps))^0.92; %% Achtung nur eine der beiden für rechnung verwenden
    NR_W_sc.m_trailing_edge = [NR_W_sc.m_trailing_edge_zylinder_actuation; NR_W_sc.m_trailing_edge_fouler_flap];
    NR_W_sc.m_slat_control = 11.23 * (area_VK_lift_div)^0.82;
    % Masse horizontal tailplane acuator
    NR_W_sc.k_hc = 0.44; % dual actuator otherwise (0.31 single actuator)
    NR_W_sc.S_he = HLW.F - (2 * HLW.F_R) ;    % Flaeche Horizontal tail plane (exposed)
    V_max = ISA.a(hoehe_CR) * specs.Ma_MO; % nicht sicher on Ma_MO oder Ma_CR ????
    deflaction_Horizontal_tailplane = deg2rad(35);
    NR_W_sc.m_HLW_control = NR_W_sc.k_hc * (NR_W_sc.S_he * (V_max^0.5) * sin(deflaction_Horizontal_tailplane))^(0.88); %%%%%%% sehr unsicher
    % speedbrake control
    NR_W_sc.s_speed_brake = 20; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Annahme
    NR_W_sc.m_speed_brake = 40.4 * NR_W_sc.s_speed_brake^0.92;
    % Lift dumper controls
    deflaction_lift_dampner = deg2rad(45);
    NR_W_sc.s_lift_dampner = 20; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Annahme
    NR_W_sc.m_lift_dampner = 20.0 * (NR_W_sc.s_lift_dampner * sin(deflaction_lift_dampner))^0.92;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Achtung hier Werte
    % %aus Tabelle 8-7 S 284 ueberpruefen !!!!!!!!!!!!!!!!!!!
    W_sc = NR_W_sc.W_sc_initial + NR_W_sc.W_sc_initial * 0.2 + NR_W_sc.W_sc_initial *0.15 + NR_W_sc.m_cockpitcontrols + NR_W_sc.m_autopilot +...
        NR_W_sc.m_maneuver_control + NR_W_sc.m_trailing_edge(1,1) + NR_W_sc.m_slat_control + NR_W_sc.m_HLW_control + NR_W_sc.m_speed_brake + NR_W_sc.m_lift_dampner; %%%%%%%%% Ergebnis in einem plausibelem bereich kein plan ob richtig (etwas hoch im vergleich 747-400)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    M_Airframe_Structur.Surface_Control_Group = W_sc;
    
    %% Engine section nacelle group (Ole)
    % nach verfahren Appendix B P. 449 zur bestimmung der umspuehlten
    % Oberflaeche
    
    % ???????????????????????????????? Quelle der Werte ???????????????????
    % deklarationen
    feet = 3.28084;     %% Umrechnungsfactor
    ln = specs.ln_TW * feet;        %%länge fan cowling
    lh = specs.lh_TW * feet;      %%länge zur größten Stelle an der urbine
    Dn = specs.Dn_TW * feet;      %%größter Durchmesser Engine
    Dh = specs.Dh_TW * feet;        %%Durchmesser Eingang Engine
    Gesamtlaenge = 4.49 * feet;
    beta_c = ln / Gesamtlaenge;   %%Gesamtlänge des Triebwerkes
    
    l_nacell_pylon = (ln - lh)/feet; % Annahme das das die laenge ist die der pylon mit der nacell interagiert

    % (B-13)
    A_cowling_sqFeet = ln * Dn * (2 + 0.35 * beta_c * (Dh / Dn) + 1.15 * (1 - beta_c));
    A_cowling = A_cowling_sqFeet * 0.092903; % umrechnung in m^2
    
    % Gas Generator section
    lg = 0.3 * feet;     %%Länge section
    Dg = 1.3 * feet;     %%größter Durchmesser section
    Deg = 0.9 * feet;    %%KLeinster Durchmesser section
    
    % Formel (B-14)
    A_generator_sqFeet = pi * lg * Dg * (1 -(1 / 3) * (1 -(Deg / Dg))) * 1 - 0.18 *((Dg / lg)^(5/3));
    A_generator = A_generator_sqFeet * 0.092903; % umrechnung in m^2
    
    
    %Plug
    lp = 0.9 * feet;    %%Länge section
    Dp = 0.5 * feet;     %%Anfangs Durchmesser
    % Formel (B-15)
    A_plug_sqFeet = 0.7 * pi * lp * Dp;
    A_plug = A_plug_sqFeet * 0.092903; % umrechnung in m^2
    
    %Alles zsm
    A_Turbine_sqFeet = (A_generator_sqFeet + A_cowling_sqFeet + A_plug_sqFeet); % Achtung Ergebnis in Feet
    
    A_Turbine = A_Turbine_sqFeet * 0.092903; % umrechnung in m^2
    
    % %Gewicht berechnen nach:
    % % Commercinal airplane design principals
    % %  Refined Weight and Balance Estimate  s.317
    % 
    % W_e = 2.20462 * specs.m_TW;       %%Gewicht Triebwerk
    % S_p = 10.3 * feet^2;              %%umspülte fläche des pilon
    % n_ult = M_Fluegel.n_ult;                    %%Danach soll itteriert werden sagt jasper
    % 
    % 
    % W_n = 35.45 * 2 *((2.33 * (1.1 * W_e) * A_Turbine_sqFeet)/(10000))^(0.59);    %%Ausrechnen in pounds
    % W_p = 24.11 * 2 * S_p^(0.381) * ((1.46 * n_ult * 1.1 * W_e * ln * Dn)/(10^6 * cos(deg2rad(80)))^(0.952));% 
    
    % Gewicht = W_n / 2.20462
    % Gewicht = W_p / 2.20462
    
    % Gewicht nach Torenbeek Tabelle 8-8 P.284
    
    m_nacell_stuktur = 0.405 * sqrt(M_Rumpf.v_D_EAS) * A_Turbine^(1.3);
    
    m_gas_gen_couling = 14.6 * (A_plug + A_generator);
    
    m_noisesuppression = 1.71 * A_cowling;
    
    M_Airframe_Structur.Engine_saction_nacelle_group = (m_nacell_stuktur + m_gas_gen_couling + m_noisesuppression) * specs.n_TW;
    
    
    %%%%%%%% Wo wird die Masse der Pylons Eingerechnet ??????????????????????????
    M_Airframe_Structur.Zwischensumme = M_Airframe_Structur.Wing_Group +...
        M_Airframe_Structur.Tail_Group + M_Airframe_Structur.Fuselage_Group +...
        M_Airframe_Structur.Undercarriage_Group + M_Airframe_Structur.Surface_Control_Group +...
        M_Airframe_Structur.Engine_saction_nacelle_group;
    
    
    
    
    %% Propulsion group
    
    M_Propulsion_Group = specs.n_TW * specs.m_TW;
    
    %% Airframe Service and equipment (Mac)
    
    M_Airframe_service_and_equipment.M_APU = specs.m_APU;
    
    % MAssen Intrumente Navigation equip und elektronik gruppe
    % Instrumente S.289
    NR_Airframe_service.W_ieg_alt = 0.347 * M_del_empty_initial^(5/9) * specs.max_range_basis_km^(1/4);
    % ????????????????????????????????????????????????????????????????????
    Technologiefaktor_W_ieg = 0.8;
    NR_Airframe_service.W_ieg = NR_Airframe_service.W_ieg_alt * Technologiefaktor_W_ieg;
    
    M_Airframe_service_and_equipment.M_intruments_nav_electr = NR_Airframe_service.W_ieg;
    
    
    % Massen Hydraulik Pneumatik und Elektronik gruppe
    % Absolut nicht sicher !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % Electronics P.290
    NR_Airframe_service.v_cabine = 62.1 * pi * (specs.R_rumpf^2) /2;
    
    % (8-43)
    NR_Airframe_service.P_el_APU = 3.64 * NR_Airframe_service.v_cabine^0.7;
    % (8-41)
    NR_Airframe_service.W_el_AC = 16.3 * NR_Airframe_service.P_el_APU * (1 - 0.033 * sqrt(NR_Airframe_service.P_el_APU));
    
    % (8-40)
    NR_Airframe_service.W_el_DC = 0.02 * M_TO_initial + 181;
    NR_Airframe_service.W_el_ges = NR_Airframe_service.W_el_AC + NR_Airframe_service.W_el_DC;
    
    
    % Hydraulic / Pneumatic S.290 (8-39)
    NR_Airframe_service.W_hyd = 0.015 * M_del_empty_initial + 272;
    
    M_Airframe_service_and_equipment.M_hydraulic_electric_group = NR_Airframe_service.W_hyd +NR_Airframe_service.W_el_ges;
    
    
    
    % Massen Furnishing and equipment group
    % Flight Deck Accomodation S.291
    NR_Airframe_service.W_flight_deck_acc = 16.5 * M_del_empty_initial^(0.285);
    
    % PAX Cabin Accomodities S.291
    
    % Pax Sitze gemäß Torenbeek Table 3-2 S.76
    % Spalten: Single, Double, Triple
    % Zeilen: Eco, Business, First Class
    NR_Airframe_service.seats = [[8 16 24]; 
        [13.6 25.4 35.4].*0.9; 
        [21.3 65 0]]; % [[10.9 21.3 29.9]; [13.6 25.4 35.4]; [21.3 31.8 0]]; orginal
    % Quelle fuer ECO: https://www.expliseat.com/products-seat-lines/ ||
    %       zum Sitz TISEAT E2 X-LINE wurden noch 1.5 kg addiert um
    %       Elektronik zu beruecksichtigen

    % Quelle fuer Bussiness: Annahme von einer Gewichtsverbesserug von 10%
    
    % Quelle fuer First: https://wingdesign.com/shop/flugzeugmoebel/flugzeugsitze/flugzeugsitz-doppelsitzbank-businessclass-leder-schwarz/
    %       Doppelsitz First Class 10 kg leichter als Modernisierungsfaktor

   
    
    % ATT Seats
    NR_Airframe_service.n_flight_att = 9;           %% Bei All Eco ist es 9
    NR_Airframe_service.W_Flight_att = 8.2 * NR_Airframe_service.n_flight_att;     
    
    
    % Galley/Pantry
    Galley_main = 113.4;
    Galley_medium = 45.3;
    Galley_small = 29.5;
    
    NR_Airframe_service.W_Galley = Galley_main * 2 + Galley_medium * 4 + Galley_small * 5;
    
    %Lavatory/ Toiletts
    lavatory = 136;
    NR_Airframe_service.W_lavatory = lavatory * 7 + lavatory * 0.2 * 1; 
    
    %Floor Covering
    NR_Airframe_service.W_floor = 1.25 * (5.85 * 57.34)^(1.15);
    
    % Sound Proofing and Insulation ??????????????????????????????????????
    V_pax = 5.85 * 57.34 * 1.676;
    V_ch  = 1.68 * 4.15 * 57.34;
    NR_Airframe_service.W_soundProof = 6.17 * (V_pax + V_ch)^(1.07);
    
    % Cargo Accomodations S.291
    % Cargo Retraint + Handling Provisions
    NR_Airframe_service.W_cargo_restraint = 1.28 * V_ch;
    NR_Airframe_service.W_cargo_palletprov = 13.67 * V_ch;
    
    % Standard Emergency Equipment S.291
    % Fixed Oxygen System
    NR_Airframe_service.W_oxygenSys = 18.1 + 1.09 * specs.n_pax; % Pax 
    
    % Fire Detection
    W_to = M_TO_initial;
    NR_Airframe_service.W_fireDetect = 0.0012 * W_to;
    
    % Evacuation
    NR_Airframe_service.W_evac = 0.453 * specs.n_pax_all_eco;
    
    M_Airframe_service_and_equipment.M_Furnishing_equip = NR_Airframe_service.W_flight_deck_acc + ...
        NR_Airframe_service.W_Flight_att + NR_Airframe_service.W_Galley + NR_Airframe_service.W_lavatory +...
        NR_Airframe_service.W_floor + NR_Airframe_service.W_soundProof + NR_Airframe_service.W_cargo_restraint + ...
        NR_Airframe_service.W_cargo_palletprov + NR_Airframe_service.W_oxygenSys + NR_Airframe_service.W_fireDetect + ...
        NR_Airframe_service.W_evac;
    
    
    % Aircon
    % (8-45) P.293
    NR_Airframe_service.W_aircon = 14 * 57.34^1.28; % Annahme l_pc =57.34m
    
    M_Airframe_service_and_equipment.M_Aircon_antifreeze = NR_Airframe_service.W_aircon;
    
    % Misscellaneous Weight
       % Wird weggelassen als Technologiefaktor, aufgrund seines geringen
       % Außmaßes von 1% Delivery empty mass
    M_Airframe_service_and_equipment.M_Miscellaneous = 0.005 * M_TO_initial;
    
    M_Airframe_service_and_equipment.Zusammen = M_Airframe_service_and_equipment.M_Miscellaneous + M_Airframe_service_and_equipment.M_Aircon_antifreeze + ...
        M_Airframe_service_and_equipment.M_Furnishing_equip + M_Airframe_service_and_equipment.M_hydraulic_electric_group + ...
        M_Airframe_service_and_equipment.M_intruments_nav_electr + M_Airframe_service_and_equipment.M_APU;
    
    
    %% Opperational Items Tabelle (8-13) S.292
    
    % Crew Provisions
    m_opp_items.m_prov_crew = 93 * specs.n_flight + 68 * specs.n_crew; 
    
    % Passenger Cabin suppys
    m_opp_items.m_computer = 0.453 * specs.n_pax;
    m_opp_items.m_snacks = 2.27 * specs.n_pax;
    m_opp_items.m_main_meal = 8.62 * specs.n_pax;
    m_opp_items.passenger_cabin_supp = m_opp_items.m_computer + m_opp_items.m_snacks + m_opp_items.m_main_meal;
    % Portable Wather / Toilet Chemicals
    m_opp_items.m_port_water = 90.7 * specs.n_Toilette;
    
    % Safy equip
    m_opp_items.m_safty_equip = 3.4 * specs.n_pax;
    
    % Residual Fuel 
    m_opp_items.m_residual_fuel = 0.151 * Tank.V_Tank^(2/3);
    
    % Seating 
    % PAX Seats                 %% Bei All ECO ist es 432 
    % m_opp_items.W_PAX_seats_ECO = (specs.n_pax / 3) * NR_Airframe_service.seats(1,3); 
    m_opp_items.W_PAX_seats_Basic = (30/2) * NR_Airframe_service.seats(3,2) + 10 * NR_Airframe_service.seats(2,3) + 16 * NR_Airframe_service.seats(2,2) + (219) *NR_Airframe_service.seats(1,1); 
                                     %First Class                           %business Class                                                              % Economy 
    m_opp_items.W_PAX_seats_AllEco = 136 * NR_Airframe_service.seats(1,3) + 12 * NR_Airframe_service.seats(1,2);
    m_opp_items.SeatDeltaBasicAllEco = m_opp_items.W_PAX_seats_AllEco - m_opp_items.W_PAX_seats_Basic;

    m_opp_items.Zusammen = m_opp_items.m_prov_crew + m_opp_items.m_computer +...
        m_opp_items.m_snacks + m_opp_items.m_main_meal + m_opp_items.m_port_water +...
        m_opp_items.m_safty_equip + m_opp_items.m_residual_fuel + m_opp_items.W_PAX_seats_Basic;
    
    
    Masse_delivery_empty = M_Airframe_service_and_equipment.Zusammen +...
        M_Propulsion_Group + M_Airframe_Structur.Zwischensumme;
    Masse_opperating_empty = m_opp_items.Zusammen + M_Airframe_service_and_equipment.Zusammen +...
        M_Propulsion_Group + M_Airframe_Structur.Zwischensumme;
    % Speicher Zwischenergebnis stand 21.5.23 nicht itteriert
    % Masse_opperating_empty = 1.480585 *e+05
    
    %% Payload Vergleiche FE1 PS.4
    
   
    M_Payload.M_payload_basis = specs.m_pax + specs.m_cargo;
    
    %% Zero Fuel Mass 
    M_Zero_Fuel.M_ZF = Masse_opperating_empty + M_Payload.M_payload_basis;
    
    
    %% takeoff Masse
    
    M_take_off_initial.test_mto = M_Zero_Fuel.M_ZF/(1-FF.Kappa_ges);
    M_take_off_initial.M_TO = M_Zero_Fuel.M_ZF+ M_take_off_initial.M_fuel;
    


    %% Check wie sich die Daten verändert haben

    delta_M_to = M_TO_initial - M_take_off_initial.M_TO;
    delta_M_to_test = M_TO_initial - M_take_off_initial.test_mto;


    %% Kopie der Werte um weiter rechnen zu können
    
    M_TO_initial = M_take_off_initial.M_TO;
    M_OE_initial = Masse_opperating_empty;
    M_del_empty_initial = Masse_delivery_empty;
    M_Zero_Fuel_initial = M_Zero_Fuel.M_ZF;
    %M_Airframe_Structur.HLW
   
    Zaehlvariabele = Zaehlvariabele + 1; % test
end     % Ende der Iteration

%% Speichern der Daten

Ergebnisse_Massen_FE2.M_TO = M_TO_initial;
Ergebnisse_Massen_FE2.M_OE = M_OE_initial;
Ergebnisse_Massen_FE2.M_DE = M_del_empty_initial;
Ergebnisse_Massen_FE2.M_ZF = M_Zero_Fuel_initial;
Ergebnisse_Massen_FE2.M_F = M_TO_initial * FF.Kappa_ges;




% Airplane Structure
Anteile_einzel_Massen_FE2.Airplane_Structure.Wing_group = M_Airframe_Structur.Wing_Group;
Anteile_einzel_Massen_FE2.Airplane_Structure.Tail_group = M_Airframe_Structur.Tail_Group;
Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.M = M_Airframe_Structur.Fuselage_Group;
Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.M_Rumpf = M_Rumpf; 
Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.Sg_12 = M_Rumpf.Sg_12;
Anteile_einzel_Massen_FE2.Airplane_Structure.Undercarriage_group = M_Airframe_Structur.Undercarriage_Group;
Anteile_einzel_Massen_FE2.Airplane_Structure.Surface_control_group = M_Airframe_Structur.Surface_Control_Group;
Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Masse = M_Airframe_Structur.Engine_saction_nacelle_group;
Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Turbine_Area = A_Turbine;
Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Turbine_Diameter = Dn / feet;
Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Turbine_length = Gesamtlaenge / feet;
Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.l_nacell_pylon = l_nacell_pylon;
Anteile_einzel_Massen_FE2.Airplane_Structure.Zusammen = M_Airframe_Structur.Zwischensumme;
Anteile_einzel_Massen_FE2.Airplane_Structure.NR.Rumpf = NR_Rumpf_Geometrie;
Anteile_einzel_Massen_FE2.Airplane_Structure.NR.FLuegel = NR_M_Fluegel;
Anteile_einzel_Massen_FE2.Airplane_Structure.MainGear = M_Airframe_Structur.MainGear;
Anteile_einzel_Massen_FE2.Airplane_Structure.FrontGear = M_Airframe_Structur.FrontGear;
% Engine Group
Anteile_einzel_Massen_FE2.Propulsion.Propulsion_group = M_Propulsion_Group;
% Service and Equipment
Anteile_einzel_Massen_FE2.Airframe_Service_equipment.APU = M_Airframe_service_and_equipment.M_APU;
Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Intruments_Nav_Electr = M_Airframe_service_and_equipment.M_intruments_nav_electr;
Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Hydraulics_Electric = M_Airframe_service_and_equipment.M_hydraulic_electric_group;
Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Furnishing_equipment = M_Airframe_service_and_equipment.M_Furnishing_equip;
Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Aircon_AntiIce = M_Airframe_service_and_equipment.M_Aircon_antifreeze;
Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Miscellaneous = M_Airframe_service_and_equipment.M_Miscellaneous;
Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Zusammen = M_Airframe_service_and_equipment.Zusammen;
Anteile_einzel_Massen_FE2.Airframe_Service_equipment.NR_Airframe_service = NR_Airframe_service;
% Opperational Items
Anteile_einzel_Massen_FE2.Opperational_Items.Crew_provi = m_opp_items.m_prov_crew;
Anteile_einzel_Massen_FE2.Opperational_Items.Passenger_cabin_supp = m_opp_items.passenger_cabin_supp;
Anteile_einzel_Massen_FE2.Opperational_Items.Liquids = m_opp_items.m_port_water;
Anteile_einzel_Massen_FE2.Opperational_Items.Safty_equip = m_opp_items.m_safty_equip;
Anteile_einzel_Massen_FE2.Opperational_Items.Seating = m_opp_items.W_PAX_seats_Basic;
Anteile_einzel_Massen_FE2.Opperational_Items.Residual_Fuel = m_opp_items.m_residual_fuel;
Anteile_einzel_Massen_FE2.Opperational_Items.Zusammen = m_opp_items.Zusammen;
Anteile_einzel_Massen_FE2.Opperational_Items.M_opp_items = m_opp_items;

% Treibstoff
Anteile_einzel_Massen_FE2.Fuel_fractions.Kappa = FF.Kappa_ges;
Anteile_einzel_Massen_FE2.Fuel_fractions.M_Fuel = Ergebnisse_Massen_FE2.M_F;
Anteile_einzel_Massen_FE2.Fuel_fractions.FF = FF;

save Ergebnisse_Massen_FE2.mat Ergebnisse_Massen_FE2 Anteile_einzel_Massen_FE2 M_HLW M_SLW FF



% sichern aller Berechnungen 
save Zwischenergebnisse_FE2_PS1.mat





