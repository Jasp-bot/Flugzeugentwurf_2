% Iteration V2 Pseudocode

clear all
close all
clc

load Ergebnisse_Widerstand_FE2.mat;

%% Anfangsbedingung

% festlegen einer startbedingung
% noch unsicher wie diese aussieht, da ich mir nicht sicher bein was die
% Endbedingung ist? Stand 18.07.2023 11Uhr



%% Beginn der Iteration

dX = 1; % fuer die while schlefe, als Anfangsbedingung
zaehlvariabele = 1; % um zu wissen, wie oft die Schleife durchlaufen wird

Startwerte_Iteration.CA_CW_LR = 1 ./ Ergebnisse_Widerstand_FE2.cW_cA_off_D;
% Startwerte_Iteration.CA_CW_TOL = ;          % Werte von Mac
% Startwerte_Iteration.CA_CW_LDG = ;         % Werte von Mac



Eingabewert_Iteration = 0; % Startwert

while abs(dX) > 0.0001
% FE1 Iteration
    % Beginn mit FE1 PS 3 Massenabschaetzung
    Berechnung_PS4_basis_stat_Massen(Eingabewert_Iteration); % eingabewert iteration entscheidet ob die WERTE  von FE1 oder FE2 verwendet werden
        % Eingabewert_Iteration =1 (Werte FE2) | Eingabewert_Iteration =0 (Werte FE1) 
    
    % FE1 PS4 Flaechenbelastung 
    
    Berechnung_PS5_Flaechenbelastung; 
        % hier kann die Flaechenbelastung und das c_A vom Flugzeug angepasst werden
    
    % FE1 PS5 Familienbildung, hier nicht sicher ob nötig
    Berechnung_PS5_familie_stat_Massen;

    % FE1 PS6 Schubanforderung
    Berechnung_PS6_Startschub_Landeanforderung(Startwerte_Iteration); % Startwerte_Iteration ist ein Struct

    
    % FE1 PS7 Fluegel/Tank
    Berechnung_PS8_Fluegel_Tank;
    % FE1 PS8 Auftriebsverteilung
    Berechnung_PS9_Auftrieb_Momente
    % FE1 PS8 Leitwerke
    Berechnung_PS9_Leitwerke
    % FE1 PS9 Widerstand
    [Startwerte_Iteration_FE1] = Berechnungen_PS10_Widerstand(Eingabewert_Iteration);
    
% FE2 Ieration
    while dx_FE2 > 0.00001

    % FE2 PS1 Neue Massenabschaetzung
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
    

    load Ergebnisse_Massen_FE2.mat;
    load Ergebnisse_Flugleistung_2.mat;

    dx_FE2 = abs(FF.mf4 - FFneu.mf4);

    end
    




    Eingabewert_Iteration = 1;
end


    % FE2 PS3 Fahrwerk
    Berechnung_FE2_PS3_Fahrwerk;

    


