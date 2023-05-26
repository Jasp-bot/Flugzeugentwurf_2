function Berechnung_PS5_Flaechenbelastung

clc
clear all
close all

%% Inputs

load Projekt_specs.mat;
load Zwischenergebnisse_PS4_basis_stat_Massen.mat;
load Ergebnisse_Basis_stat_m.mat;
load Ergebnisse_ISA_DATA.mat


% Matrix der Fluegeldaten (gegeben) bereits konvertiert in .mat Datei
load importdata.mat

% excel_fluegeldaten.Flugzeugnamen = flaechenbelastungsData(:,1) 
excel_fluegeldaten.m_ToM = flaechenbelastungsData(:,2); % max take off mass
excel_fluegeldaten.F = flaechenbelastungsData(:,3);        % fleuegelflaeche
excel_fluegeldaten.G_To_F = flaechenbelastungsData(:,4);   % Flaechenbelastung
excel_fluegeldaten.Ma_mo = flaechenbelastungsData(:,5);    % Machzahl mo
excel_fluegeldaten.n_E = flaechenbelastungsData(:,6);      % Anzahl Triebwerke



%% Regression zur Bestimmung Fluegelflaeche u. Flaechenbelastung

% Regressionsanalyse Flaeche

zwischenergebnis_ps5.f_F = fit(excel_fluegeldaten.m_ToM, excel_fluegeldaten.F, 'power1');
zwischenergebnis_ps5.coeffizienten_F = coeffvalues(zwischenergebnis_ps5.f_F);
zwischenergebnis_ps5.A_F = zwischenergebnis_ps5.coeffizienten_F(1,1);
zwischenergebnis_ps5.B_F = zwischenergebnis_ps5.coeffizienten_F(1,2);

zwischenergebnis_ps5.Fluegelflaeche_FZ_regression = zwischenergebnis_ps5.A_F .* Ergebnis_basis_m.m_To .^ zwischenergebnis_ps5.B_F; %%%%%%%%%%%%%%%%%% wert ist zu klein um eine gute Faechenbalastung zu bekommen...
% deswegen wird die Fluegelflaeche erhoeht bei constantem Ergebnis_basis_m.m_To

Ergebnisse_stat_Flaechenbelastung.F = 435;



% Regressionanalyse Flaechenbelastung

zwischenergebnis_ps5.f_G_To_F = fit(excel_fluegeldaten.m_ToM, excel_fluegeldaten.G_To_F, 'power1');
zwischenergebnis_ps5.coeffizienten_G_To_F = coeffvalues(zwischenergebnis_ps5.f_G_To_F);
zwischenergebnis_ps5.A_G_To_F = zwischenergebnis_ps5.coeffizienten_G_To_F(1,1);
zwischenergebnis_ps5.B_G_To_F = zwischenergebnis_ps5.coeffizienten_G_To_F(1,2);

zwischenergebnis_ps5.Fleachenbelastung_FZ_regression = zwischenergebnis_ps5.A_G_To_F .* Ergebnis_basis_m.m_To .^ zwischenergebnis_ps5.B_G_To_F;

% Berechung der neuen Flachenbelastung

Ergebnisse_stat_Flaechenbelastung.G_To = Ergebnis_basis_m.m_To * specs.g;
Ergebnisse_stat_Flaechenbelastung.G_To_ICA = Ergebnis_basis_m.m_To * specs.g * 0.98;

Ergebnisse_stat_Flaechenbelastung.Fleachenbelastung = Ergebnisse_stat_Flaechenbelastung.G_To_ICA ...
    / Ergebnisse_stat_Flaechenbelastung.F; 



%% Berechnung Auftrieb ueber Fluegelflaeche

Flughoehe = specs.flight_level * 10^2 ;                         % in ft
hoehe =round(unitsratio('m','ft')*Flughoehe);                     % in m


Ergebnisse_stat_Flaechenbelastung.C_A_CR = Ergebnisse_stat_Flaechenbelastung.Fleachenbelastung * ...
    1/((ISA.rho(hoehe,1)/2) * (ISA.a(hoehe,1)* specs.Ma_CR)^2);



%% Speichern der Daten in .mat
save Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche Ergebnisse_stat_Flaechenbelastung
save Zwischenergebnisse_PS5_Fluegelflaechen excel_fluegeldaten zwischenergebnis_ps5













