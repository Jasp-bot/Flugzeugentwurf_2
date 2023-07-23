%% Iteration

clc
clear all
close all

Startwerte_Iteration = Berechnungen_PS10_Widerstand;

%% Vorbereitung / Veranschaulichung

% E_CR = Startwerte_Iteration.CA_CW_LR;
% E_TO = Startwerte_Iteration.CA_CW_TO;
% E_LDG = Startwerte_Iteration.CA_CW_LDG;

%% Beginn der Iteration

dX = 1;
%zaehlvariabele_test = 0; 

while abs(dX) > 0.0001
    Berechnung_PS6_Startschub_Landeanforderung(Startwerte_Iteration);
    Berechnung_PS8_Fluegel_Tank;
    Berechnung_PS9_Auftrieb_Momente;
    Berechnung_PS9_Leitwerke;
    Startwerte_Iteration = Berechnungen_PS10_Widerstand;
    Endwerte_Iteration = Berechnungen_PS10_Widerstand;
    % Vergleiche
    load Ergebnisse_Start_Landeanforderungen.mat;
    dX_CR = Endwerte_Iteration.CA_CW_LR - schub_CR.Eta;
    dX_TO = Endwerte_Iteration.CA_CW_TO - startschub.Eta_To_inv(3,1);
    dX_LDG = Endwerte_Iteration.CA_CW_LDG - (1/landeanvorderung.Eta_LDG);

    dX = dX_CR + dX_TO + dX_LDG;
    %zaehlvariabele_test = zaehlvariabele_test +1;


end

Endwerte_Iteration = Berechnungen_PS10_Widerstand;


save Ergebnisse_Endwerte_Iteration_V1.mat Endwerte_Iteration

%HAllo 
