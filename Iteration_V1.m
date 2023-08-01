%% Iteration v1

clc
clear all
close all

x=0;


%% Vorbereitung / Veranschaulichung

% E_CR = Startwerte_Iteration.CA_CW_LR;
% E_TO = Startwerte_Iteration.CA_CW_TO;
% E_LDG = Startwerte_Iteration.CA_CW_LDG;

%% Beginn der Iteration

dX = 1;
zaehlvariabele_test = 0; 

 % Beginn mit FE1 PS 3 Massenabschaetzung
    Berechnung_PS4_basis_stat_Massen(0); % eingabewert iteration entscheidet ob die WERTE  von FE1 oder FE2 verwendet werden
        % Eingabewert_Iteration =1 (Werte FE2) | Eingabewert_Iteration =0 (Werte FE1) 
    
    % FE1 PS4 Flaechenbelastung 
    
    Berechnung_PS5_Flaechenbelastung; 
        % hier kann die Flaechenbelastung und das c_A vom Flugzeug angepasst werden
    
    % FE1 PS5 Familienbildung, hier nicht sicher ob nötig
    Berechnung_PS5_familie_stat_Massen;
Startwerte_Iteration = Berechnungen_PS10_Widerstand(x);
while abs(dX) > 0.0001
    Berechnung_PS6_Startschub_Landeanforderung(Startwerte_Iteration,x);
    Berechnung_PS8_Fluegel_Tank;
    Berechnung_PS9_Auftrieb_Momente;
    Berechnung_PS9_Leitwerke;
    Startwerte_Iteration = Berechnungen_PS10_Widerstand(x);
    Endwerte_Iteration = Berechnungen_PS10_Widerstand(x);
    % Vergleiche
    load Ergebnisse_Start_Landeanforderungen.mat;
    dX_CR = Endwerte_Iteration.CA_CW_LR - schub_CR.Eta;
    dX_TO = Endwerte_Iteration.CA_CW_TO - startschub.Eta_To_inv(3,1);
    dX_LDG = Endwerte_Iteration.CA_CW_LDG - (1/landeanvorderung.Eta_LDG);

    dX = dX_CR + dX_TO + dX_LDG;
    zaehlvariabele_test = zaehlvariabele_test +1;


end

Endwerte_Iteration = Berechnungen_PS10_Widerstand(x);


save Ergebnisse_Endwerte_Iteration_V1.mat Endwerte_Iteration

%HAllo 
