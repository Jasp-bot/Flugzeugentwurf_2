%% Visualisierung der Daten zur Flächenbelastung und Fluegelflaeche

clc
clear all 
close all

% Ausfuehrunf der Funktion
Berechnung_PS5_Flaechenbelastung

%% Inputs
load Projekt_specs.mat;
load Ergebnisse_Basis_stat_m.mat;
load Ergebnisse_ISA_DATA.mat
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche
load Zwischenergebnisse_PS5_Fluegelflaechen.mat

%% Visualisierung

figure(1)
% Fluegelflaeche
hold on
grid on
plot(zwischenergebnis_ps5.f_F, excel_fluegeldaten.m_ToM, excel_fluegeldaten.F, '.r');
plot(Ergebnis_basis_m.m_To, zwischenergebnis_ps5.Fluegelflaeche_FZ_regression, 'xb');
plot(Ergebnis_basis_m.m_To, Ergebnisse_stat_Flaechenbelastung.F, '*k');

title(['Fluegelflaeche über Maximales Startgewicht']);
xlabel('Maximales Startgewicht in kg');
ylabel('Fluegelflaeche');
legend('Datenpunkte', 'Regressionskurve', 'Abgeschätzte Flaeche', 'Gewählte Fläche');

figure(2)
% Fleachenbelastung
hold on 
grid on
plot(zwischenergebnis_ps5.f_G_To_F, 'b', excel_fluegeldaten.m_ToM, excel_fluegeldaten.G_To_F); %, m_To_Flugzeug, Fleachenbelastung_FZ, 'x')
plot(Ergebnis_basis_m.m_To, zwischenergebnis_ps5.Fleachenbelastung_FZ_regression, 'xr');
plot(Ergebnis_basis_m.m_To, Ergebnisse_stat_Flaechenbelastung.Fleachenbelastung, '*black');

title(['Flaechnbelastung über Maximales Startgewicht']);
xlabel('Maximales Startgewicht in kg');
ylabel('Flaechenbelastung');
legend('Datenpunkte', 'Regressionskurve', 'Abgeschätzte Flaechenbelastung','Berechnete Flaechenbelastung', 'Location', 'southeast');



