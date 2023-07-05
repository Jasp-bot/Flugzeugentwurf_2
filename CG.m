%% CG Berechnung PS02
clc
clear all
close all


load Projekt_specs.mat;
load Ergebnisse_Massen_FE2.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Auftrieb_Momente.mat;
load Ergebnisse_Leitwerke.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;






%% getroffene Annahmen
c_A_alpha_F = GRA.c_a_anstieg;      % Annahme bitte ueberpruefen !!!!!!!!!!!!!!!!!!!!
z_abstand = 3; % Abstand zwischen Profilsehnen angenommen vergleiche Torenbeek s480

l_fn=40;            %%Abstand Flugzeugnase zu Fluegel Schnittpunkt mit Rumpf NICHT MAC!!! Muss noch ver�ndert werden
Durchmesser_Flugzeug_wo_Fluegel_durchgeht=specs.D_rumpf;      % Erste annahme!!!! Der fluegel geht nicht in der mitte durch den Rumpf sondern etwas weiter unten daher ist der rumpf da etwas schmaler



CG_Rumpf_Prozent=[Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.M,0.43,0,0];
CG_Hoenleitwerk_Prozent=0.96;
CG_Seitenleitwerk_Prozent=0.96;
CG_Bugfahrwerk_Prozent=0.5;
CG_APU_Prozent=[specs.m_APU,0.98,0,0];
CG_Instruments_and_co_Prozent=[Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Intruments_Nav_Electr, 0.7,0,0];
CG_Hydraulics_Wing_Prozent=0.50;
CG_Hydraulics_Tail_Prozent=0.90;
CG_AC_AntiICE_Prozent=[Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Aircon_AntiIce, 0.5,0,-2];
CG_Miscellaneous_Prozent=[Anteile_einzel_Massen_FE2.Airframe_Service_equipment.Miscellaneous, 0.5,0,0];
CG_Crew_Provisions_Prozent=[Anteile_einzel_Massen_FE2.Opperational_Items.Crew_provi, 0.49,0,0];
CG_Passenger_Supplies_Prozent=[Anteile_einzel_Massen_FE2.Opperational_Items.Crew_provi, 0.55,0,0];
CG_Water_Toilet_Chemicals_Prozent=[Anteile_einzel_Massen_FE2.Opperational_Items.Liquids, 0.80,0,-1];
CG_Safety_equipment_Prozent=[Anteile_einzel_Massen_FE2.Opperational_Items.Safty_equip, 0.5,0,0.5];
CG_Seatings_prozent=[Anteile_einzel_Massen_FE2.Opperational_Items.Seating, 0.58,0,0];




%%CG_Prozent_residual_fuel_Prozent=[Anteile_einzel_Massen_FE2.Opperational_Items.Residual_Fuel,
%%0,0,-2.5];      muss im Schwerpunkt liegen daher noch kein Wert



% [Masse, x, y, z];
% seats_MAT=[Anteile_einzel_Massen_FE2.Opperational_Items.Seating,0.55*specs.l_rumpf,0,0];%% Z h�he liegt direkt im Rumpfmittelpunkt Y=0 weil symetrisch
% cockpit_MAT=[4500,0.05*specs.l_rumpf,0,0];
% tail_MAT=[Anteile_einzel_Massen_FE2.Airplane_Structure.Tail_group];
% wing_MAT=[Anteile_einzel_Massen_FE2.Airplane_Structure.Wing_group];

CG_gesamt_Mat=[CG_Rumpf_Prozent;CG_APU_Prozent;CG_Instruments_and_co_Prozent;CG_AC_AntiICE_Prozent;CG_Miscellaneous_Prozent;
    CG_Crew_Provisions_Prozent;CG_Passenger_Supplies_Prozent;CG_Water_Toilet_Chemicals_Prozent;CG_Safety_equipment_Prozent;
    CG_Seatings_prozent];   %%Es sind noch nicht alle einzelteile des Flugzeuges enthalten
CG_M_und_x=0;
CG_M=0;


for C=1:length(CG_gesamt_Mat)       %%Der Schwerpunkt wird ausgrechnet mit der Schwerpunkt Formel
    CG_M_und_x=CG_M_und_x+(CG_gesamt_Mat(C,1)*(CG_gesamt_Mat(C,2)*specs.l_rumpf));
    CG_M=CG_M+CG_gesamt_Mat(C,1);
end
CG_ges=CG_M_und_x/CG_M;              %% Ges Schwerpunkt von der Flugzeugnase aus in Meter   
%%specs.l_rumpf

F_BFW_min=specs.m_k_min*specs.g*CG_M;   %% min Bugfahrwerkslast zur Bestimmung des Schwerpunktes

%%Aus berechnung von Jasper oder PS2 (Bestimmung des Neutralpunktes)

Abwindfaktor = 1.75 * ((c_A_alpha_F)/(pi * Ergebnisse_Fluegel.streckung_phi25_max *...
    (Ergebnisse_Fluegel.lambda * (HLW.r/(Ergebnisse_Fluegel.b/2))^0.25 *...
    (1+ (abs(z_abstand))/((Ergebnisse_Fluegel.b/2))) )));


c_A_alpha_H = (pi * HLW.streckung_phi25)/(1 + sqrt(1 + 0.25 * HLW.streckung_phi25^2 * (tan(HLW.phi_50)^2 + (1 - specs.Ma_CR^2))));

c_A_alpha = c_A_alpha_F * (1+ ((c_A_alpha_H)/(c_A_alpha_F)) *...
    (HLW.F/Ergebnisse_stat_Flaechenbelastung.F) * (1 - Abwindfaktor) );



d_F1_X_NP_durch_l_mue=(1.8/c_A_alpha_F)*(specs.D_rumpf^2 *l_fn)/(Ergebnisse_Fluegel.F*Ergebnisse_Fluegel.l_mue)    %%Einfluss des rumpfes vor und hinter dem Fluegel Formel 28

d_F2_X_NP_durch_l_mue=(0.273*specs.D_rumpf^2*Ergebnisse_Fluegel.l_m*(Ergebnisse_Fluegel.b+specs.D_rumpf))/((1+Ergebnisse_Fluegel.lambda)*Ergebnisse_Fluegel.l_mue^2*(Ergebnisse_Fluegel.b+2.15*specs.D_rumpf))*tan(Ergebnisse_Fluegel.phi_25_max) %Einfluss des Rumpf Fl�gel �berganges auf den NP FORMEL 29 PS02

d_TW_X_NP_durch_l_mue=(specs.Dn_TW^2 *specs.l_TW)/(Ergebnisse_Fluegel.F*Ergebnisse_Fluegel.l_mue*c_A_alpha_F) %Einfluss Triebwerk auf NP Formel 30 PS02

l_i_Mitte= DT.l_i_I+(tan(Ergebnisse_Fluegel.phi_VK_max)*Durchmesser_Flugzeug_wo_Fluegel_durchgeht*0.5)    %Ausrechnen von tiefe des Fluges imagin�r IM rumpf     

Flaeche_im_Rumpf_oberes_dreieck=(l_i_Mitte-DT.l_i_I)*Durchmesser_Flugzeug_wo_Fluegel_durchgeht*0.5;


%Gesamte Fluegelflaeche mit dem Dreieck im Rumpf
F_ges_Fluegel_MAC=Flaeche_im_Rumpf_oberes_dreieck+Ergebnisse_Fluegel.F;
%gesamte streckung mit teil im rumpf (dreieck)
Streckung_ges=(Ergebnisse_Fluegel.b^2)/F_ges_Fluegel_MAC;

%%NP des Fl�gels berechnen Formel 24 PS02
X_NP_F=l_i_Mitte*(0.25+(Streckung_ges/12)*(1+2*Ergebnisse_Fluegel.lambda*tan(Ergebnisse_Fluegel.phi_25_max)));


X_NP_OH_durch_l_mue=(X_NP_F/Ergebnisse_Fluegel.l_mue)+d_F1_X_NP_durch_l_mue+d_F2_X_NP_durch_l_mue+d_TW_X_NP_durch_l_mue


%X_NP_OH_durch_l_mue*Ergebnisse_Fluegel.l_mue

X_NP_durch_l_mue=X_NP_OH_durch_l_mue+((HLW.r/Ergebnisse_Fluegel.l_mue)*0.85*c_A_alpha_H/c_A_alpha*(1-Abwindfaktor))


Neutralpunkt=X_NP_durch_l_mue*Ergebnisse_Fluegel.l_mue+l_fn



save Ergebnisse_CG.mat Neutralpunkt F_ges_Fluegel_MAC Flaeche_im_Rumpf_oberes_dreieck


