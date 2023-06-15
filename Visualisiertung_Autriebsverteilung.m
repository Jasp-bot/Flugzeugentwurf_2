% Visualisierung Auftriebsverteilung und Fluegelmomente

clc
clear all
close all

%% Ausfuerung des Programms

Berechnung_PS9_Auftrieb_Momente;

%% Laden der Dateien

load Projekt_specs.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Auftrieb_Momente.mat;





%% plot der Zirkulationsanteile

[x1, y1]= max(Ergebnisse_Auftriebsverteilung.c_a_eta);

figure(1)
hold on
grid on
plot(Ergebnisse_Auftriebsverteilung.eta,GRA.gamma_a_eta)
plot(Ergebnisse_Auftriebsverteilung.eta,VWA.gamma_b_eta)
plot(Ergebnisse_Auftriebsverteilung.eta,Ergebnisse_Auftriebsverteilung.gamma_eta_ges,'r')
plot(Ergebnisse_Auftriebsverteilung.eta,Ergebnisse_Auftriebsverteilung.c_a_eta)
plot([0 1], [Ergebnisse_stat_Flaechenbelastung.C_A_CR, Ergebnisse_stat_Flaechenbelastung.C_A_CR],'--k')
plot(y1*10^(-3),x1,'*r')
legend('Grundrissabhaengiger Anteil', 'Verwindungsabhaengiger Anteil', 'Gesamtzirkulation',...
    'Auftriebsbeiwertverteilung', 'Reiseauftriebsbeiwert C_a','Maximaler Auftriebsbeiwert','FontSize',18)

title('Auftriebsbeiwertverteilung','FontSize',20)
xlabel('Dimensionslose Halbspannweite','FontSize',16)

hold off
