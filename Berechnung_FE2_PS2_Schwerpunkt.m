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

Wing_MAC.xSP_MAC_lmue = 0.21; % siehe Übung: Wert zwischen 20% und 25%
Wing_MAC.xSP_MAC_FG = CG_Wing_X - (CG_Data_Wing.Fluegel(2)-NP.l_mue_ges);

Wing_MAC.XMAC = CG_Rumpf_X + Wing_MAC.xSP_MAC_FG*(CG_Wing_M/CG_M) - Wing_MAC.xSP_MAC_lmue*(1+(CG_Wing_M/CG_M))*NP.l_mue_ges;

%% Bestimmung des Gesamtschwerpunktes

Rumpf_MAC.xSP_MAC_RG = -Wing_MAC.XMAC + CG_Rumpf_X;

CG_Gesamt_x = (Rumpf_MAC.xSP_MAC_RG*CG_M + Wing_MAC.xSP_MAC_FG*CG_Wing_M)/(CG_M + CG_Wing_M);

%% Flügellage
% Von Flugzeugnase zu imaginärer Spitze im Rumpf
Wing_Position1 = Wing_MAC.XMAC - CG_Data_Wing.Fluegel(2);
% Von Flugzeugnase zur geraden Sektion im Rumpf
Wing_Position2 = Wing_MAC.XMAC - CG_Data_Wing.Fluegel(2) + DT.s_R*tan(Ergebnisse_Fluegel.phi_VK_max);

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
CG_Fuel_X.Aussentrapez_MAC = CG_Fuel_X.Aussentrapez_FG - (CG_Data_Wing.Fluegel(2)-NP.l_mue_ges);
CG_Fuel_X.Innentrapez_MAC = CG_Fuel_X.Innentrapez_FG - (CG_Data_Wing.Fluegel(2)-NP.l_mue_ges);
CG_Fuel_X.Rumpf_MAC = CG_Fuel_X.Rumpf_FG - (CG_Data_Wing.Fluegel(2)-NP.l_mue_ges);

CG_Fuel_X.CG_Tank_MAC = CG_Fuel_X.CG_Tank_FG - (CG_Data_Wing.Fluegel(2)-NP.l_mue_ges);

%% Bestimmung Schwerpunkt Frachträume

% CAD???

%% Beladung

% FUEL
% Treibstoffmasse
Betankung.Masse_Fuel_Aussen = 2*Tank.V_OB_A*Tank.ka*Tank.kb*1000*Tank.rho_kerosin;
Betankung.Masse_Fuel_Innen_Theoretisch = 2*Tank.V_OB_I*Tank.ka*Tank.kb*1000*Tank.rho_kerosin;
Betankung.Masse_Fuel_Innen_Praktisch = Ergebnisse_Massen_FE2.M_F-Betankung.Masse_Fuel_Aussen;
% Betankung Außentank
Betankung.CG_BetankterAussentank = (CG_Gesamt_x*Ergebnisse_Massen_FE2.M_OE + CG_Fuel_X.Aussentrapez_MAC*Betankung.Masse_Fuel_Aussen)/(Ergebnisse_Massen_FE2.M_OE + Betankung.Masse_Fuel_Aussen);
Betankung.P1 = [CG_Gesamt_x/NP.l_mue_ges; Ergebnisse_Massen_FE2.M_OE];
Betankung.P2 = [Betankung.CG_BetankterAussentank/NP.l_mue_ges; Ergebnisse_Massen_FE2.M_OE + Betankung.Masse_Fuel_Aussen];
% Betankung Innentank
Betankung.CG_BetankterInnentank = (Betankung.CG_BetankterAussentank*(Ergebnisse_Massen_FE2.M_OE + Betankung.Masse_Fuel_Aussen) + CG_Fuel_X.Innentrapez_MAC*Betankung.Masse_Fuel_Innen_Praktisch)/(Betankung.Masse_Fuel_Innen_Praktisch + Ergebnisse_Massen_FE2.M_OE + Betankung.Masse_Fuel_Aussen);
Betankung.P3 = [Betankung.CG_BetankterInnentank/NP.l_mue_ges; (Betankung.Masse_Fuel_Innen_Praktisch + Ergebnisse_Massen_FE2.M_OE + Betankung.Masse_Fuel_Aussen)];


% PAX Außen
Pax_Bording = readtable("Bording.xlsx","Sheet","Beladung");
AussenPax = table2array(Pax_Bording(1:49,"x"));
AussenPaxMAC = AussenPax - Wing_MAC.XMAC;
MasseAussenPax = table2array(Pax_Bording(1:49,"m"));
% Front to Back
CG_Shift_Outer = zeros(length(AussenPaxMAC),1);
NewMassCounter = zeros(length(AussenPaxMAC),1);
Mass_Shift_Outer = Ergebnisse_Massen_FE2.M_OE + Ergebnisse_Massen_FE2.M_F;
CG_Startposition = Betankung.CG_BetankterInnentank;
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
BackwardsCG_Startposition = Betankung.CG_BetankterInnentank;
for i = 1:length(BackwardsAussenPaxMAC)
    BackwardsnewCGvalue = (BackwardsCG_Startposition*BackwardsMass_Shift_Outer + BackwardsAussenPaxMAC(i)*BackwardsMasseAussenPax(i))/(BackwardsMass_Shift_Outer + BackwardsMasseAussenPax(i)); % Hier können Sie den neuen Wert je nach Bedarf berechnen
    BackwardsCG_Shift_Outer(i) = BackwardsnewCGvalue; % Wert an der entsprechenden Stelle im Vektor zuweisen
    BackwardsNewMassCounter(i) = BackwardsMass_Shift_Outer + BackwardsMasseAussenPax(i);
    BackwardsCG_Startposition = BackwardsnewCGvalue;
    BackwardsMass_Shift_Outer = BackwardsMass_Shift_Outer + BackwardsMasseAussenPax(i);
end

%% Plotten

figure(1)
hold on 
grid on
xlim([0 60])
ylim([Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000])

plot([Betankung.P1(1)*100 Betankung.P2(1)*100], [Betankung.P1(2) Betankung.P2(2)],"b-")
plot([Betankung.P2(1)*100 Betankung.P3(1)*100], [Betankung.P2(2) Betankung.P3(2)],"bo-")
plot(CG_Shift_Outer*100/NP.l_mue_ges, NewMassCounter,"rx-")
plot(BackwardsCG_Shift_Outer*100/NP.l_mue_ges, BackwardsNewMassCounter,"mx-")

title('Beladung 3-Klassenbestuhlung','FontSize',20);
xlabel('X^{MAC}_{SP}/l_{\mu}','FontSize',16)
ylabel('kg','FontSize',16)