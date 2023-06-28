%% PS2 Schwerpunktberechnung

clc
clear all
close all

%% Laden von Werten
load Projekt_specs.mat;
load Ergebnisse_Massen_FE2.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Auftrieb_Momente.mat;
load Ergebnisse_Leitwerke.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;

%% Berechnung Rumpfschwerpunkt

% Annahme von Prozentwerten für x_SP bezogen auf Gesamtlänge
Rumpf_SP_Faktoren.xSP_Rumpf = 0.43;
Rumpf_SP_Faktoren.xSP_HLW = 0.96;
Rumpf_SP_Faktoren.xSP_SLW = 0.96;
Rumpf_SP_Faktoren.xSP_Bugfahrwerk = 0.05;

Rumpf_SP_Faktoren.xSP_APU = 0.98;
Rumpf_SP_Faktoren.xSP_CockpitInstruments = 0.07;
Rumpf_SP_Faktoren.xSP_HydraulicsElectricalWing = 0.5;
Rumpf_SP_Faktoren.xSP_HydraulicsElectricalTail = 0.9;
Rumpf_SP_Faktoren.xSP_AC_AntiIce = 0.5;
Rumpf_SP_Faktoren.xSP_Misc = 0.5;

Rumpf_SP_Faktoren.xSP_CrewProvisions = 0.49;
Rumpf_SP_Faktoren.xSP_PassengerCabinSupplies = 0.55;
Rumpf_SP_Faktoren.xSP_WaterToiletChem = 0.8;
Rumpf_SP_Faktoren.xSP_SafetyEq = 0.5;
Rumpf_SP_Faktoren.xSP_Seating = 0.58;
%Rumpf_SP_Faktoren.xSP_ResFuel = Flügelschwerpunkt;

% Annahme von Prozentwerten für z_SP in Metern
Rumpf_SP_Faktoren.zSP_Rumpf = 0;
Rumpf_SP_Faktoren.zSP_HLW = 1;
Rumpf_SP_Faktoren.zSP_SLW = 2.5;
Rumpf_SP_Faktoren.zSP_Bugfahrwerk = -4.5;

Rumpf_SP_Faktoren.zSP_APU = 0;
Rumpf_SP_Faktoren.zSP_CockpitInstruments = -0.5;
Rumpf_SP_Faktoren.zSP_HydraulicsElectricalWing = -2;
Rumpf_SP_Faktoren.zSP_HydraulicsElectricalTail = 0;
Rumpf_SP_Faktoren.zSP_AC_AntiIce = -2;
Rumpf_SP_Faktoren.zSP_Misc = 0;

Rumpf_SP_Faktoren.zSP_CrewProvisions = 0;
Rumpf_SP_Faktoren.zSP_PassengerCabinSupplies = 0;
Rumpf_SP_Faktoren.zSP_WaterToiletChem = -1;
Rumpf_SP_Faktoren.zSP_SafetyEq = 0.5;
Rumpf_SP_Faktoren.zSP_Seating = 0;
Rumpf_SP_Faktoren.zSP_ResFuel =-2.5;

% Berechnung Momente und CGX, CGZ
% [Masse,x,y,z]

CG_Data.Rumpf = [Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.M + Anteile_einzel_Massen_FE2.Airplane_Structure.Tail_group - M_HLW.W_HLW_basic - M_SLW.W_SLW_basic,0.43,0,0];
CG_Data.HLW = [M_HLW.W_HLW_basic, 0.96, 0, 1];
CG_Data.SLW = [M_SLW.W_SLW_basic, 0.96, 0, 5];
CG_Data.Bugfahrwerk = [Anteile_einzel_Massen_FE2.Airplane_Structure.FrontGear, 0.05, 0, -4.5];
CG_Data.APU = [specs.m_APU, 0.97, 0, 0];
CG_Data.CockpitInstruments = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Intruments_Nav_Electr, 0.07,0,0];
CG_Data.HydraulicsElectricalWing = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Hydraulics_Electric*0.6,0.45, 0, -2];
CG_Data.HydraulicsElectricalTail = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Hydraulics_Electric*0.4,0.9, 0, 0];
CG_Data.Furnishing = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Furnishing_equipment, 0.5, 0, 0];
CG_Data.AC_AntiIce = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Aircon_AntiIce, 0.5,0,-2];
CG_Data.Misc = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Miscellaneous, 0.5,0,0];
CG_Data.CrewProvisions = [Anteile_einzel_Massen_FE2.Opperational_Items.Crew_provi, 0.49,0,0];
CG_Data.PassengerCabinSupplies = [Anteile_einzel_Massen_FE2.Opperational_Items.Passenger_cabin_supp, 0.55,0,0];
CG_Data.WaterToiletChem = [Anteile_einzel_Massen_FE2.Opperational_Items.Liquids, 0.80,0,-1];
CG_Data.SafetyEq = [Anteile_einzel_Massen_FE2.Opperational_Items.Safty_equip, 0.5,0,0.5];
CG_Data.Seating = [Anteile_einzel_Massen_FE2.Opperational_Items.Seating, 0.58,0,0];
CG_Data.ResFuel = [Anteile_einzel_Massen_FE2.Opperational_Items.Residual_Fuel, 0.5, 0, -2.5];

CG_DataMatrix=[CG_Data.Rumpf;CG_Data.HLW;CG_Data.SLW;CG_Data.Bugfahrwerk;CG_Data.APU;CG_Data.CockpitInstruments;...
    CG_Data.HydraulicsElectricalWing;CG_Data.HydraulicsElectricalTail;CG_Data.Furnishing;CG_Data.AC_AntiIce;CG_Data.Misc;...
    CG_Data.CrewProvisions;CG_Data.PassengerCabinSupplies;CG_Data.WaterToiletChem;CG_Data.SafetyEq;CG_Data.Seating;CG_Data.ResFuel];
%ResFuel Annahme SP bei x=0.5

CG_Moment_X=0;
CG_M=0;
for C=1:length(CG_DataMatrix)       %%Der Schwerpunkt wird ausgrechnet mit der Schwerpunkt Formel
    CG_Moment_X=CG_Moment_X+(CG_DataMatrix(C,1)*(CG_DataMatrix(C,2)*specs.l_rumpf));
    CG_M=CG_M+CG_DataMatrix(C,1);
end
CG_Rumpf_X=CG_Moment_X/CG_M;
CG_Rumpf_X_Prozent = CG_Rumpf_X / specs.l_rumpf;

CG_Moment_Z=0;
CG_MZ=0;
for C=1:length(CG_DataMatrix)       %%Der Schwerpunkt wird ausgrechnet mit der Schwerpunkt Formel
    CG_Moment_Z=CG_Moment_Z+(CG_DataMatrix(C,1)*(CG_DataMatrix(C,4)));
    CG_MZ=CG_MZ+CG_DataMatrix(C,1);
end

CG_Rumpf_Z=CG_Moment_Z/CG_MZ; 


%% Berechnung Flügelschwerpunkt
CG_Data_Wing.Fluegel = [Anteile_einzel_Massen_FE2.Airplane_Structure.Wing_group, (NP.versatz_HK + DT.l_i_R + tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)-NP.x_SP_ges, 0, 1]; 
CG_Data_Wing.MainGear = [Anteile_einzel_Massen_FE2.Airplane_Structure.MainGear, 8, 0, -3]; 
CG_Data_Wing.SurfaceControls = [Anteile_einzel_Massen_FE2.Airplane_Structure.Surface_control_group, (NP.versatz_HK + DT.l_i_R + tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)-NP.x_SP_ges, 0, 1]; 
CG_Data_Wing.EngineSection = [Anteile_einzel_Massen_FE2.Propulsion.Propulsion_group, 4.5, 0, -2.5]; 
CG_Data_Wing.Nacelle = [Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Masse, 5.5, 0, 0.5]; 
 
CG_DataMatrix_Wing = [CG_Data_Wing.Fluegel;CG_Data_Wing.MainGear;CG_Data_Wing.SurfaceControls;CG_Data_Wing.EngineSection;CG_Data_Wing.Nacelle]; 
 
CG_Wing_Moment_X=0; 
CG_Wing_M=0; 
for C=1:length(CG_DataMatrix_Wing)       %%Der Schwerpunkt wird ausgrechnet mit der Schwerpunkt Formel 
    CG_Wing_Moment_X=CG_Wing_Moment_X+(CG_DataMatrix_Wing(C,1)*(CG_DataMatrix_Wing(C,2))); 
    CG_Wing_M=CG_Wing_M+CG_DataMatrix_Wing(C,1); 
end 
CG_Wing_X=CG_Wing_Moment_X/CG_Wing_M; 

CG_Wing_Moment_Z=0; 
CG_Wing_MZ=0; 
for C=1:length(CG_DataMatrix_Wing)       %%Der Schwerpunkt wird ausgrechnet mit der Schwerpunkt Formel 
    CG_Wing_Moment_Z=CG_Wing_Moment_Z+(CG_DataMatrix_Wing(C,1)*(CG_DataMatrix_Wing(C,4))); 
    CG_Wing_MZ=CG_MZ+CG_DataMatrix_Wing(C,1); 
end 
CG_Wing_Z=CG_Wing_Moment_Z/CG_Wing_MZ;

%% Bestimmung von X_MAC

Wing_MAC.xSP_MAC_lmue = 0.23; % siehe Übung: Wert zwischen 20% und 25%
Wing_MAC.xSP_MAC_FG = CG_Wing_X - (CG_Data_Wing.Fluegel(2)-NP.l_mue_ges);

Wing_MAC.XMAC = CG_Rumpf_X + Wing_MAC.xSP_MAC_FG*(CG_Wing_M/CG_M) - Wing_MAC.xSP_MAC_lmue*(1+(CG_Wing_M/CG_M))*NP.l_mue_ges;

%% Bestimmung des Gesamtschwerpunktes

Rumpf_MAC.xSP_MAC_RG = -Wing_MAC.XMAC + CG_Rumpf_X;

CG_Gesamt_x = (Rumpf_MAC.xSP_MAC_RG*CG_M + Wing_MAC.xSP_MAC_FG*CG_Wing_M)/(CG_M + CG_Wing_M);