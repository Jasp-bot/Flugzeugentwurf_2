% recovery fuer FE1
clc
clear all
close all


Startwerte_Iteration.CA_CW_LR = 17.3;
Startwerte_Iteration.CA_CW_TO = 10.5;
Startwerte_Iteration.CA_CW_LDG = 7.4;


Berechnung_PS4_basis_stat_Massen(0);
Berechnung_PS5_Flaechenbelastung
Berechnung_PS6_Startschub_Landeanforderung(Startwerte_Iteration, 0);

Berechnung_PS8_Fluegel_Tank;
Berechnung_PS9_Auftrieb_Momente;
Berechnung_PS9_Leitwerke;
Ausgabewert = Berechnungen_PS10_Widerstand(0);

% FE2
Berechnung_FE2_PS1_M_Mf;
    
% FE2 PS2 Schwepunkt
Berechnung_FE2_PS2_Schwerpunkt
        
% FE2 PS4 Widerstand
Berechnung_FE2_PS4_Widerstand;
    
% FE2 PS5 Hochauftrieb 1
Berechnung_FE2_PS5_Hochauftrieb_1;
    
% FE2 PS6 Hochauftrieb 2
 Berechnung_FE2_PS6_Hochauftrieb_2
    
% FE2 PS7 Flugleistung 1
Berechnung_FE2_PS7_Flugleistung1;
    
% FE 2 PS8 Flugleistung 2
Berechnung_FE2_PS8_Flugleistung2_NRD;



