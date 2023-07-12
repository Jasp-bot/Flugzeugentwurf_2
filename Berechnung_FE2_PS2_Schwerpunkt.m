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
% [Masse,x(%),y,z]

CG_Data.Rumpf = [Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.M + Anteile_einzel_Massen_FE2.Airplane_Structure.Tail_group - M_HLW.W_HLW_basic - M_SLW.W_SLW_basic,0.42,0,0];
CG_Data.HLW = [M_HLW.W_HLW_basic, 0.93, 0, 1];
CG_Data.SLW = [M_SLW.W_SLW_basic, 0.95, 0, 5];
CG_Data.Bugfahrwerk = [Anteile_einzel_Massen_FE2.Airplane_Structure.FrontGear, 0.06, 0, -4.5];
CG_Data.APU = [specs.m_APU, 0.97, 0, 1];
CG_Data.CockpitInstruments = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Intruments_Nav_Electr, 0.035,0,-1.5];
CG_Data.HydraulicsElectricalWing = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Hydraulics_Electric*0.6,0.3, 0, -2];
CG_Data.HydraulicsElectricalTail = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Hydraulics_Electric*0.4,0.7, 0, -1];
CG_Data.Furnishing = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Furnishing_equipment, 0.5, 0, 2];
CG_Data.AC_AntiIce = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Aircon_AntiIce, 0.5,0,-2];
CG_Data.Misc = [Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Miscellaneous, 0.5,0,0];
CG_Data.CrewProvisions = [Anteile_einzel_Massen_FE2.Opperational_Items.Crew_provi, 0.48,0,0];
CG_Data.PassengerCabinSupplies = [Anteile_einzel_Massen_FE2.Opperational_Items.Passenger_cabin_supp, 0.52,0,0];
CG_Data.WaterToiletChem = [Anteile_einzel_Massen_FE2.Opperational_Items.Liquids, 0.80,0,-1];
CG_Data.SafetyEq = [Anteile_einzel_Massen_FE2.Opperational_Items.Safty_equip, 0.5,0,0.5];
CG_Data.Seating = [Anteile_einzel_Massen_FE2.Opperational_Items.Seating, 0.51,0,0];
CG_Data.ResFuel = [Anteile_einzel_Massen_FE2.Opperational_Items.Residual_Fuel, 0.49, 0, -2.5];

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
% Achtung: Alles im Flügelkoordinatensystem
CG_Data_Wing.Fluegel = [Anteile_einzel_Massen_FE2.Airplane_Structure.Wing_group, (NP.versatz_HK + DT.l_i_R + tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)-NP.x_SP_ges-3, 0, 1]; 
CG_Data_Wing.MainGear = [Anteile_einzel_Massen_FE2.Airplane_Structure.MainGear, DT.s_R*tan(DT.phi_VK_max)+0.82*DT.l_i_R, 0, -3]; 
CG_Data_Wing.SurfaceControls = [Anteile_einzel_Massen_FE2.Airplane_Structure.Surface_control_group, (NP.versatz_HK + DT.l_i_R + tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)-NP.x_SP_ges, 0, 0.5]; 
CG_Data_Wing.EngineSection = [Anteile_einzel_Massen_FE2.Propulsion.Propulsion_group, 4, 0, -1.5]; 
CG_Data_Wing.Nacelle = [Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Masse, 4.5, 0, -0.5]; 
 
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
    CG_Wing_MZ=CG_Wing_MZ+CG_DataMatrix_Wing(C,1); 
end 
CG_Wing_Z_FG=CG_Wing_Moment_Z/CG_Wing_MZ;
CG_Wing_Z_RG=CG_Wing_Z_FG-2.19;

%% Bestimmung von X_MAC

Wing_MAC.xSP_MAC_lmue = 0.35; % siehe Übung: Wert zwischen 20% und 25%

Wing_MAC.xSP_MAC_FG = 0.5*NP.l_mue_ges -(CG_Data_Wing.Fluegel(2)-CG_Wing_X);

Wing_MAC.XMAC = CG_Rumpf_X + Wing_MAC.xSP_MAC_FG*(CG_Wing_M/CG_M) - Wing_MAC.xSP_MAC_lmue*(1+(CG_Wing_M/CG_M))*NP.l_mue_ges;

%% Bestimmung des Gesamtschwerpunktes

Rumpf_MAC.xSP_MAC_RG = -Wing_MAC.XMAC + CG_Rumpf_X;

CG_Gesamt_x = (Rumpf_MAC.xSP_MAC_RG*CG_M + Wing_MAC.xSP_MAC_FG*CG_Wing_M)/(CG_M + CG_Wing_M);
CG_Gesamt_z = (CG_Rumpf_Z*CG_MZ+CG_Wing_Z_RG*CG_Wing_MZ)/(CG_MZ+CG_Wing_MZ);

%% Flügellage
% Von Flugzeugnase zu imaginärer Spitze im Rumpf
Wing_Position1 = Wing_MAC.XMAC - CG_Data_Wing.Fluegel(2) +0.5*NP.l_mue_ges;
% Von Flugzeugnase zur geraden Sektion im Rumpf
Wing_Position2 = Wing_MAC.XMAC - CG_Data_Wing.Fluegel(2) + DT.s_R*tan(Ergebnisse_Fluegel.phi_VK_max)+0.5*NP.l_mue_ges;

%% Bestimmung Schwerpunkt Tank

% Schwerpunkt Außentrapez
CG_Fuel_X.Epsilon_Aussen = Tank.c_aussen*tan(Ergebnisse_Fluegel.phi_VK_max) + Tank.b1_A;
CG_Fuel_X.Aussentrapez = (Tank.b1_I^2 - Tank.b1_A^2 + CG_Fuel_X.Epsilon_Aussen*(Tank.b1_I + 2*Tank.b1_A))/(3*(Tank.b1_A+Tank.b1_I));
CG_Fuel_X.Flaeche_Aussentrapez = 0.5*(Tank.b1_I + Tank.b1_A)*Tank.c_aussen;

% Schwerpunkt Innentrapez
CG_Fuel_X.Epsilon_Innen = Tank.c_innen*tan(Ergebnisse_Fluegel.phi_VK_max) + Tank.b1_I;
CG_Fuel_X.Innentrapez = (Tank.b1_R^2 - Tank.b1_I^2 + CG_Fuel_X.Epsilon_Innen*(Tank.b1_R + 2*Tank.b1_I))/(3*(Tank.b1_R+Tank.b1_I));
CG_Fuel_X.Flaeche_Innentrapez = 0.5*(Tank.b1_I + Tank.b1_R)*Tank.c_innen;

%Schwerpunkt Rumpfteil
CG_Fuel_X.Rumpf = 0.5*Tank.b1_R;
CG_Fuel_X.Flaeche_Rumpf = Tank.b1_R*Tank.c_rumpf;

% Gesamtschwerpunkt Tank im Flügelkoordinatensystem
CG_Fuel_X.DeltaTankWing = DT.l_i_R*0.15 + DT.s_R*tan(Ergebnisse_Fluegel.phi_VK_max);
CG_Fuel_X.Aussentrapez_FG = CG_Fuel_X.Aussentrapez + Tank.c_innen*tan(Ergebnisse_Fluegel.phi_VK_max) + CG_Fuel_X.DeltaTankWing;
CG_Fuel_X.Innentrapez_FG = CG_Fuel_X.Innentrapez + CG_Fuel_X.DeltaTankWing;
CG_Fuel_X.Rumpf_FG = CG_Fuel_X.Rumpf + CG_Fuel_X.DeltaTankWing;

CG_Fuel_X.CG_Tank_FG = (CG_Fuel_X.Flaeche_Aussentrapez*CG_Fuel_X.Aussentrapez_FG + CG_Fuel_X.Flaeche_Innentrapez*CG_Fuel_X.Innentrapez_FG + CG_Fuel_X.Flaeche_Rumpf*CG_Fuel_X.Rumpf_FG)/(CG_Fuel_X.Flaeche_Rumpf + CG_Fuel_X.Flaeche_Innentrapez + CG_Fuel_X.Flaeche_Aussentrapez);

% Gesamtschwerpunkt Tank im MAC-Koordinatensystem
CG_Fuel_X.Aussentrapez_MAC = 0.5*NP.l_mue_ges -(CG_Data_Wing.Fluegel(2)-CG_Fuel_X.Aussentrapez_FG);
CG_Fuel_X.Innentrapez_MAC = 0.5*NP.l_mue_ges -(CG_Data_Wing.Fluegel(2)-CG_Fuel_X.Innentrapez_FG);
CG_Fuel_X.Rumpf_MAC = 0.5*NP.l_mue_ges -(CG_Data_Wing.Fluegel(2)-CG_Fuel_X.Rumpf_FG);

CG_Fuel_X.CG_Tank_MAC = 0.5*NP.l_mue_ges -(CG_Data_Wing.Fluegel(2)-CG_Fuel_X.CG_Tank_FG);

%% Bestimmung Schwerpunkt Frachträume
CG_Fracht.X_vorn_Prozent = 0.27;
CG_Fracht.X_hinten_Prozent = 0.63;
CG_Fracht.Masse_vorn_theoretisch = 14000;
CG_Fracht.Masse_hinten_theoretisch = 7344;
CG_Fracht.FrachtVornMac = CG_Fracht.X_vorn_Prozent*specs.l_rumpf - Wing_MAC.XMAC;
CG_Fracht.FrachtHintenMac = CG_Fracht.X_hinten_Prozent*specs.l_rumpf - Wing_MAC.XMAC;


%% Beladung

% FUEL
% Treibstoffmasse
Betankung.Masse_Fuel_Aussen_Theoretisch = 2*Tank.V_OB_A*Tank.ka*Tank.kb*1000*Tank.rho_kerosin;
Betankung.Masse_Fuel_Innen = 2*Tank.V_OB_I*Tank.ka*Tank.kb*1000*Tank.rho_kerosin;
Betankung.Masse_Fuel_Aussen_Praktisch = Ergebnisse_Massen_FE2.M_F-Betankung.Masse_Fuel_Innen;

% Betankung Außentank
Betankung.CG_BetankterInnentank = (CG_Gesamt_x*Ergebnisse_Massen_FE2.M_OE + CG_Fuel_X.Innentrapez_MAC*Betankung.Masse_Fuel_Innen)/(Ergebnisse_Massen_FE2.M_OE + Betankung.Masse_Fuel_Innen);
Betankung.P1 = [CG_Gesamt_x/NP.l_mue_ges; Ergebnisse_Massen_FE2.M_OE];
Betankung.P2 = [Betankung.CG_BetankterInnentank/NP.l_mue_ges; Ergebnisse_Massen_FE2.M_OE + Betankung.Masse_Fuel_Innen];
% Betankung Innentank
Betankung.CG_BetankterAussentank = (Betankung.CG_BetankterInnentank*(Ergebnisse_Massen_FE2.M_OE + Betankung.Masse_Fuel_Innen) + CG_Fuel_X.Aussentrapez_MAC*Betankung.Masse_Fuel_Aussen_Praktisch)/(Betankung.Masse_Fuel_Aussen_Praktisch + Ergebnisse_Massen_FE2.M_OE + Betankung.Masse_Fuel_Innen);
Betankung.P3 = [Betankung.CG_BetankterAussentank/NP.l_mue_ges; (Betankung.Masse_Fuel_Aussen_Praktisch + Ergebnisse_Massen_FE2.M_OE + Betankung.Masse_Fuel_Innen)];


% PAX AUßEN
Pax_Bording = readtable("Bording.xlsx","Sheet","Beladung");
AussenPax = table2array(Pax_Bording(1:38,"x"));
AussenPaxMAC = AussenPax - Wing_MAC.XMAC;
MasseAussenPax = table2array(Pax_Bording(1:38,"m"));
% Front to Back
CG_Shift_Outer = zeros(length(AussenPaxMAC),1);
NewMassCounter = zeros(length(AussenPaxMAC),1);
Mass_Shift_Outer = Ergebnisse_Massen_FE2.M_OE + Ergebnisse_Massen_FE2.M_F;
CG_Startposition = Betankung.CG_BetankterAussentank;
for i = 1:length(AussenPaxMAC)
    newCGvalue = (CG_Startposition*Mass_Shift_Outer + AussenPaxMAC(i)*MasseAussenPax(i))/(Mass_Shift_Outer + MasseAussenPax(i)); % Hier können Sie den neuen Wert je nach Bedarf berechnen
    CG_Shift_Outer(i) = newCGvalue; % Wert an der entsprechenden Stelle im Vektor zuweisen
    NewMassCounter(i) = Mass_Shift_Outer + MasseAussenPax(i);
    CG_Startposition = newCGvalue;
    Mass_Shift_Outer = Mass_Shift_Outer + MasseAussenPax(i);
end
% Back to Front
BackwardsAussenPaxMAC = flipud(AussenPaxMAC);
BackwardsMasseAussenPax = flipud(MasseAussenPax);
BackwardsCG_Shift_Outer = zeros(length(AussenPaxMAC),1);
BackwardsNewMassCounter = zeros(length(AussenPaxMAC),1);
BackwardsMass_Shift_Outer = Ergebnisse_Massen_FE2.M_OE + Ergebnisse_Massen_FE2.M_F;
BackwardsCG_Startposition = Betankung.CG_BetankterAussentank;
for i = 1:length(BackwardsAussenPaxMAC)
    BackwardsnewCGvalue = (BackwardsCG_Startposition*BackwardsMass_Shift_Outer + BackwardsAussenPaxMAC(i)*BackwardsMasseAussenPax(i))/(BackwardsMass_Shift_Outer + BackwardsMasseAussenPax(i)); % Hier können Sie den neuen Wert je nach Bedarf berechnen
    BackwardsCG_Shift_Outer(i) = BackwardsnewCGvalue; % Wert an der entsprechenden Stelle im Vektor zuweisen
    BackwardsNewMassCounter(i) = BackwardsMass_Shift_Outer + BackwardsMasseAussenPax(i);
    BackwardsCG_Startposition = BackwardsnewCGvalue;
    BackwardsMass_Shift_Outer = BackwardsMass_Shift_Outer + BackwardsMasseAussenPax(i);
end

% PAX INNEN
InnenPax = table2array(Pax_Bording(1:42,"x_1"));
InnenPaxMAC = InnenPax - Wing_MAC.XMAC;
MasseInnenPax = table2array(Pax_Bording(1:42,"m_1"));
% Front to Back
CG_Shift_Inner = zeros(length(InnenPaxMAC),1);
InnerMassCounter = zeros(length(InnenPaxMAC),1);
Mass_Shift_Inner = BackwardsMass_Shift_Outer;
CG_Startposition_Innen = BackwardsCG_Startposition;
for i = 1:length(InnenPaxMAC)
    newCGvalue1 = (CG_Startposition_Innen*Mass_Shift_Inner + InnenPaxMAC(i)*MasseInnenPax(i))/(Mass_Shift_Inner + MasseInnenPax(i)); % Hier können Sie den neuen Wert je nach Bedarf berechnen
    CG_Shift_Inner(i) = newCGvalue1; % Wert an der entsprechenden Stelle im Vektor zuweisen
    InnerMassCounter(i) = Mass_Shift_Inner + MasseInnenPax(i);
    CG_Startposition_Innen = newCGvalue1;
    Mass_Shift_Inner = Mass_Shift_Inner + MasseInnenPax(i);
end
% Back to Front
BackwardsInnenPaxMAC = flipud(InnenPaxMAC);
BackwardsMasseInnenPax = flipud(MasseInnenPax);
BackwardsCG_Shift_Inner = zeros(length(InnenPaxMAC),1);
BackwardsInnerMassCounter = zeros(length(InnenPaxMAC),1);
BackwardsMass_Shift_Inner = BackwardsMass_Shift_Outer;
BackwardsCG_Startposition_Innen = BackwardsCG_Startposition;
for i = 1:length(BackwardsInnenPaxMAC)
    BackwardsnewCGvalue1 = (BackwardsCG_Startposition_Innen*BackwardsMass_Shift_Inner + BackwardsInnenPaxMAC(i)*BackwardsMasseInnenPax(i))/(BackwardsMass_Shift_Inner + BackwardsMasseInnenPax(i)); % Hier können Sie den neuen Wert je nach Bedarf berechnen
    BackwardsCG_Shift_Inner(i) = BackwardsnewCGvalue1; % Wert an der entsprechenden Stelle im Vektor zuweisen
    BackwardsInnerMassCounter(i) = BackwardsMass_Shift_Inner + BackwardsMasseInnenPax(i);
    BackwardsCG_Startposition_Innen = BackwardsnewCGvalue1;
    BackwardsMass_Shift_Inner = BackwardsMass_Shift_Inner + BackwardsMasseInnenPax(i);
end

% FRACHT
% Front to back
CG_Fracht.FrachtmasseTotal = Ergebnisse_Massen_FE2.M_ZF-Ergebnisse_Massen_FE2.M_OE-(BackwardsMass_Shift_Inner-Betankung.P3(2));
CG_Fracht.CG_BeladenVorn = (BackwardsCG_Startposition_Innen*BackwardsMass_Shift_Inner + CG_Fracht.FrachtVornMac*CG_Fracht.Masse_vorn_theoretisch)/(BackwardsMass_Shift_Inner+CG_Fracht.Masse_vorn_theoretisch);
CG_Fracht.Masse_BeladenVorn = BackwardsMass_Shift_Inner+CG_Fracht.Masse_vorn_theoretisch;

CG_Fracht.Masse_hinten_real = CG_Fracht.FrachtmasseTotal - CG_Fracht.Masse_vorn_theoretisch;
CG_Fracht.CG_BeladenHinten = (CG_Fracht.CG_BeladenVorn*CG_Fracht.Masse_BeladenVorn + CG_Fracht.FrachtHintenMac*CG_Fracht.Masse_hinten_real)/(CG_Fracht.Masse_BeladenVorn+CG_Fracht.Masse_hinten_real);
CG_Fracht.Masse_BeladenHinten = CG_Fracht.Masse_BeladenVorn+CG_Fracht.Masse_hinten_real;

% Back to front
CG_Fracht.CG_BeladenHinten1 = (BackwardsCG_Startposition_Innen*BackwardsMass_Shift_Inner + CG_Fracht.FrachtHintenMac*CG_Fracht.Masse_hinten_theoretisch)/(BackwardsMass_Shift_Inner+CG_Fracht.Masse_hinten_theoretisch);
CG_Fracht.Masse_BeladenHinten1 = BackwardsMass_Shift_Inner+CG_Fracht.Masse_hinten_theoretisch;

CG_Fracht.Masse_vorn_real = CG_Fracht.FrachtmasseTotal - CG_Fracht.Masse_hinten_theoretisch;
CG_Fracht.CG_BeladenVorn1 = (CG_Fracht.CG_BeladenHinten1*CG_Fracht.Masse_BeladenHinten1 + CG_Fracht.FrachtVornMac*CG_Fracht.Masse_vorn_real)/(CG_Fracht.Masse_BeladenHinten1+CG_Fracht.Masse_vorn_real);
CG_Fracht.Masse_BeladenVorn1 = CG_Fracht.Masse_BeladenHinten1+CG_Fracht.Masse_vorn_real;

%% Berechnung Grenzen
% LÄNGSSTABILITÄT AM BODEN
LS.x_MainGear_MAC = (0.5*NP.l_mue_ges -(CG_Data_Wing.Fluegel(2)-CG_Data_Wing.MainGear(2)))/NP.l_mue_ges; %[Prozent l_mue]
LS.l_MainGear = 2.5;
LS.delta = deg2rad(15); % Min 15deg
% LS in Prozent
LS.Laengsstabilitaet = LS.x_MainGear_MAC -(LS.l_MainGear+0.5*specs.D_rumpf-CG_Gesamt_z)*(tan(LS.delta)/NP.l_mue_ges);

% MINIMALE BUGFAHRWERKSLAST
BFWL.x_CG_BFW_Min_MAC = 0.06*(CG_Data.Bugfahrwerk(2)*specs.l_rumpf - Wing_MAC.XMAC + ((1/0.06)-1)*LS.x_MainGear_MAC*NP.l_mue_ges);
BFWL.x_CG_BFW_Min_MAC_Prozent = BFWL.x_CG_BFW_Min_MAC/NP.l_mue_ges;

% MAXIMALE BUGFAHRWERKSLAST
BFWL.m_to_max = Ergebnisse_Massen_FE2.M_TO;
BFWL.MomentanMasse = linspace(Ergebnisse_Massen_FE2.M_TO,Ergebnisse_Massen_FE2.M_ZF,100)';
% ???? Warum negativ
BFWL.x_CG_BFW_Max_MAC = (BFWL.m_to_max./BFWL.MomentanMasse).*0.2.*(CG_Data.Bugfahrwerk(2)*specs.l_rumpf - Wing_MAC.XMAC + ((BFWL.MomentanMasse./BFWL.m_to_max).*((1/0.2)-1).*LS.x_MainGear_MAC.*NP.l_mue_ges));
BFWL.x_CG_BFW_Max_MAC_Prozent = BFWL.x_CG_BFW_Max_MAC./NP.l_mue_ges;

% KIPPSTABILITÄT
KS.s = 12.7/2;
KS.psi = deg2rad(60);
KS.h_sp = (311*80*0.3 + CG_Gesamt_z*Ergebnisse_Massen_FE2.M_OE)/(Ergebnisse_Massen_FE2.M_OE+311*80) + LS.l_MainGear + 0.5*specs.D_rumpf;
KS.arctan = atan(KS.s/(Wing_MAC.XMAC-(CG_Data.Bugfahrwerk(2)*specs.l_rumpf)+(LS.x_MainGear_MAC*NP.l_mue_ges)));
KS.Kippstabilitaet = (KS.h_sp/(tan(KS.psi)*sin(KS.arctan))) - Wing_MAC.XMAC + (CG_Data.Bugfahrwerk(2)*specs.l_rumpf);
KS.Kippstabilitaet_Prozent = KS.Kippstabilitaet/NP.l_mue_ges;

% STATISCHE STABILITÄT
StatStab.c_A_alpha_F = VWA.c_AF_anstieg; 
Umrechnung_Z_FG_RG = 2.19;
StatStab.z_abstand = -CG_Data_Wing.Fluegel(4)+Umrechnung_Z_FG_RG + CG_Data.HLW(4); % Abstand zwischen Profilsehnen angenommen vergleiche Torenbeek s480

l_fn=Wing_Position2;            %%Abstand Flugzeugnase zu Fluegel Schnittpunkt mit Rumpf NICHT MAC!!! Muss noch ver�ndert werden
Durchmesser_Flugzeug_wo_Fluegel_durchgeht=2*sqrt((0.5*specs.D_rumpf)^2-Umrechnung_Z_FG_RG^2); 
r_H = specs.l_rumpf-Wing_Position2-specs.coanlaenge+specs.HLW_beginn;
Abwindfaktor = 1.75 * (StatStab.c_A_alpha_F/(pi * Ergebnisse_Fluegel.streckung_phi25_max *...
    (Ergebnisse_Fluegel.lambda * (r_H/(Ergebnisse_Fluegel.b/2)))^0.25 *...
    (1+ (abs(StatStab.z_abstand/(Ergebnisse_Fluegel.b/2))))));


c_A_alpha_H = (pi * HLW.streckung_phi25)/(1 + sqrt(1 + 0.25 * HLW.streckung_phi25^2 * (tan(HLW.phi_50)^2 + (1 - specs.Ma_CR^2))));

c_A_alpha = StatStab.c_A_alpha_F * (1+ ((c_A_alpha_H)/(StatStab.c_A_alpha_F)) *...
    (HLW.F/Ergebnisse_stat_Flaechenbelastung.F) * 0.85 * (1 - Abwindfaktor) );



d_F1_X_NP_durch_l_mue=(-1.8/StatStab.c_A_alpha_F)*(specs.D_rumpf*specs.h_rumpf *l_fn)/(Ergebnisse_Fluegel.F*Ergebnisse_Fluegel.l_mue);   %%Einfluss des rumpfes vor und hinter dem Fluegel Formel 28

%-------- Ab hier noch Kontrollieren
d_F2_X_NP_durch_l_mue= tan(DT.phi_VK_max)*(0.273/(1+Ergebnisse_Fluegel.lambda))*((specs.D_rumpf*Ergebnisse_Fluegel.l_m*(specs.D_rumpf+Ergebnisse_Fluegel.b))/((Ergebnisse_Fluegel.l_mue^2)*(Ergebnisse_Fluegel.b+2.15*specs.D_rumpf)));

d_TW_X_NP_durch_l_mue=-4*(specs.Dn_TW^2 *specs.l_TW)/(Ergebnisse_Fluegel.F*Ergebnisse_Fluegel.l_mue*StatStab.c_A_alpha_F); %Einfluss Triebwerk auf NP Formel 30 PS02

l_i_Mitte= DT.l_i_I+(tan(Ergebnisse_Fluegel.phi_VK_max)*DT.s_R)    %Ausrechnen von tiefe des Fluges imagin�r IM rumpf     

Flaeche_im_Rumpf_oberes_dreieck=tan(Ergebnisse_Fluegel.phi_VK_max)*DT.s_R*DT.s_R;


%Gesamte Fluegelflaeche mit dem Dreieck im Rumpf
F_ges_Fluegel_MAC=Flaeche_im_Rumpf_oberes_dreieck+Ergebnisse_Fluegel.F;
%gesamte streckung mit teil im rumpf (dreieck)
Streckung_ges=(Ergebnisse_Fluegel.b^2)/F_ges_Fluegel_MAC;

%%NP des Fl�gels berechnen Formel 24 PS02
X_NP_F = tan(DT.phi_VK_max)*(DT.s_A+DT.s_I+DT.s_R) + DT.l_a - NP.x_NP_ges;
X_NP_F_MAC = 0.5*NP.l_mue_ges -(CG_Data_Wing.Fluegel(2)-X_NP_F)

X_NP_OH_durch_l_mue = (X_NP_F_MAC/Ergebnisse_Fluegel.l_mue) + d_F1_X_NP_durch_l_mue + d_F2_X_NP_durch_l_mue + d_TW_X_NP_durch_l_mue;


X_NP_durch_l_mue = X_NP_OH_durch_l_mue +((r_H/Ergebnisse_Fluegel.l_mue)*0.85*(HLW.F/Ergebnisse_Fluegel.F)*(c_A_alpha_H/c_A_alpha)*(1-Abwindfaktor))


Neutralpunkt=X_NP_durch_l_mue*Ergebnisse_Fluegel.l_mue+l_fn



% NEUTRALPUNKTLAGE
StatStab.CG_sigma_x = X_NP_durch_l_mue -0.05;

%% Plotten

figure(1)
hold on 
grid on
xlim([20 60])
ylim([Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000])

plot([Betankung.P1(1)*100 Betankung.P2(1)*100], [Betankung.P1(2) Betankung.P2(2)],"b-")
plot([Betankung.P2(1)*100 Betankung.P3(1)*100], [Betankung.P2(2) Betankung.P3(2)],"bo-")
plot(CG_Shift_Outer*100/NP.l_mue_ges, NewMassCounter,"rx-")
plot(BackwardsCG_Shift_Outer*100/NP.l_mue_ges, BackwardsNewMassCounter,"mx-")
plot(CG_Shift_Inner*100/NP.l_mue_ges, InnerMassCounter,"rx-")
plot(BackwardsCG_Shift_Inner*100/NP.l_mue_ges, BackwardsInnerMassCounter,"mx-")
plot([BackwardsCG_Startposition_Innen*100/NP.l_mue_ges, CG_Fracht.CG_BeladenVorn*100/NP.l_mue_ges],[BackwardsMass_Shift_Inner, CG_Fracht.Masse_BeladenVorn],"g-")
plot([CG_Fracht.CG_BeladenVorn*100/NP.l_mue_ges, CG_Fracht.CG_BeladenHinten*100/NP.l_mue_ges],[CG_Fracht.Masse_BeladenVorn, CG_Fracht.Masse_BeladenHinten],"go-")
plot([BackwardsCG_Startposition_Innen*100/NP.l_mue_ges, CG_Fracht.CG_BeladenHinten1*100/NP.l_mue_ges],[BackwardsMass_Shift_Inner, CG_Fracht.Masse_BeladenHinten1],"g--")
plot([CG_Fracht.CG_BeladenHinten1*100/NP.l_mue_ges, CG_Fracht.CG_BeladenVorn1*100/NP.l_mue_ges],[CG_Fracht.Masse_BeladenHinten1, CG_Fracht.Masse_BeladenVorn1],"go--")

%Grenzen
plot([BFWL.x_CG_BFW_Min_MAC_Prozent*100,BFWL.x_CG_BFW_Min_MAC_Prozent*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"-k")
plot(BFWL.x_CG_BFW_Max_MAC_Prozent*100,BFWL.MomentanMasse,"--k")
plot([LS.Laengsstabilitaet*100,LS.Laengsstabilitaet*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"-b")
plot([StatStab.CG_sigma_x*100,StatStab.CG_sigma_x*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],":b")

title('Beladung 3-Klassenbestuhlung','FontSize',20);
xlabel('X^{MAC}_{SP}/l_{\mu}','FontSize',16)
ylabel('kg','FontSize',16)
legend('','','','','','','','','','','Minimale BFWL','Maximale BFWL','Längsstabilität','NP -5%')

hold off

%-------------------------------------
% Plot Schwerpunktlage
figure(2)
hold on
grid on
% RUMPF
% Mittelpunkt der Ellipse
center_x = specs.l_rumpf/2;
center_y = 0;
% Erstellen der Ellipse
theta = linspace(0, 2*pi, 100); % Winkel für die Ellipsenpunkte
x = center_x + (specs.l_rumpf/2) * cos(theta);
y = center_y + (specs.D_rumpf/2) * sin(theta);

% ERSATZFLÜGEL
x1 = Wing_MAC.XMAC;
y1 = -24;
x2 = Wing_MAC.XMAC + NP.l_mue_ges;
y2 = 24;


% PLOTTEN
% Rumpf
plot(x, y, 'r', 'LineWidth', 2);
% Ersatzflügel
rectangle('Position', [x1, y1, x2-x1, y2-y1], 'EdgeColor', 'red', 'LineWidth', 0.5);
% Flügel
plot([Wing_Position1, Wing_Position1 + (DT.s_A+DT.s_I+DT.s_R)*tan(DT.phi_VK_max)], [0, (DT.s_A+DT.s_I+DT.s_R)], 'r', 'LineWidth', 2);
plot([Wing_Position1, Wing_Position1 + (DT.s_A+DT.s_I+DT.s_R)*tan(DT.phi_VK_max)], [0, -(DT.s_A+DT.s_I+DT.s_R)], 'r', 'LineWidth', 2);
plot([Wing_Position1+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max), Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max)],[DT.s_R+DT.s_I, DT.s_R+DT.s_I], 'r', 'LineWidth', 2);
plot([Wing_Position1+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max), Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max)],[-(DT.s_R+DT.s_I), -(DT.s_R+DT.s_I)], 'r', 'LineWidth', 2);
plot([Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max),DT.l_a+Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max)],[-(DT.s_R+DT.s_I+DT.s_A),-(DT.s_R+DT.s_I+DT.s_A)], 'r', 'LineWidth', 2);
plot([Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max),DT.l_a+Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max)],[DT.s_R+DT.s_I+DT.s_A,DT.s_R+DT.s_I+DT.s_A], 'r', 'LineWidth', 2);
plot([Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max),DT.l_a+Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max)],[-(DT.s_R+DT.s_I),-(DT.s_R+DT.s_I+DT.s_A)], 'r', 'LineWidth', 2);
plot([Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max),DT.l_a+Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max)],[DT.s_R+DT.s_I,DT.s_R+DT.s_I+DT.s_A], 'r', 'LineWidth', 2);
plot([Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max),Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max)],[DT.s_R+DT.s_I,-(DT.s_R+DT.s_I)], 'r', 'LineWidth', 2);
% Schwerpunkte Tank
% plot([CG_Fuel_X.Aussentrapez+(DT.s_I*tan(DT.phi_VK_max))+Wing_Position1+DT.s_R*tan(DT.phi_VK_max)+(0.15*DT.l_i_R),CG_Fuel_X.Aussentrapez+(DT.s_I*tan(DT.phi_VK_max))+Wing_Position1+DT.s_R*tan(DT.phi_VK_max)+(0.15*DT.l_i_R)],[-30,30], '-b');
%plot([Wing_MAC.XMAC+CG_Fuel_X.CG_Tank_MAC,Wing_MAC.XMAC+CG_Fuel_X.CG_Tank_MAC],[-30,30], '--m');
%plot([Wing_Position1+CG_Fuel_X.CG_Tank_FG,Wing_Position1+CG_Fuel_X.CG_Tank_FG],[-30,30], ':g');

%plot([Wing_MAC.XMAC+LS.x_MainGear_MAC*NP.l_mue_ges,Wing_MAC.XMAC+LS.x_MainGear_MAC*NP.l_mue_ges],[-30,30], '--m');
%plot([Wing_Position1+CG_Data_Wing.MainGear(2),Wing_Position1+CG_Data_Wing.MainGear(2)],[-30,30], ':g');

% CGs
plot(CG_Rumpf_X,0,'xb');
plot(Wing_Position1+CG_Wing_X,0,'ob');
plot(LS.x_MainGear_MAC*NP.l_mue_ges +Wing_MAC.XMAC, 0, 'xk');
plot(CG_Data.Bugfahrwerk(2)*specs.l_rumpf,0,'*k');
plot(CG_Gesamt_x + Wing_MAC.XMAC,0,'og','LineWidth', 3)
plot([Wing_Position1+CG_Data_Wing.Fluegel(2),Wing_Position1+CG_Data_Wing.Fluegel(2)],[-30,30],'-k')
plot([Wing_MAC.XMAC+(X_NP_durch_l_mue*NP.l_mue_ges),Wing_MAC.XMAC+(X_NP_durch_l_mue*NP.l_mue_ges)],[-30,30],'--g');

axis equal;
legend('','','','','','','','','','','CG Rumpf','CG Flügelgruppe','CG HFW', 'CG BFW', 'CG Gesamt')

% Achsenbeschriftungen und Titel
xlabel('x [m]');
ylabel('y [m]');
title('Schwerpunktlagen');
