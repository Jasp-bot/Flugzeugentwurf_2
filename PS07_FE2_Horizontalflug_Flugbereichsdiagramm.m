%%PS 07 Parameterstudie zum Horizontalflug und Flugbereichsdiagramm

clc
clear all
close all


load Ergebnisse_Massen_FE2.mat
load Projekt_specs.mat
load Ergebnisse_ISA_DATA.mat
load Ergebnisse_Widerstand.mat
load CaVerteilungFE12.mat
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat
load Ergebnisse_Start_Landeanforderungen.mat
load Ergebnisse_Widerstand_FE2.mat



%OD Polare
%Variabeln
G=Ergebnisse_Massen_FE2.M_TO*specs.g;       %Max Abfluggewicht * g (9.80665)
Flughoehe_CR = specs.flight_level * 10^2 ;     % in ft
hoehe_CR = round(unitsratio('m','ft')*Flughoehe_CR);    %umrechnung in m
T=ISA.T;    %Temperatur
a=ISA.a;    %MachzahlISA
Druckdifferenz=ISA.p(hoehe_CR)/ISA.p0;

roh_zu_roh0=ISA.rho(hoehe_CR)/ISA.rho_0;    %Dichteverhältnis
CA_j=Ergebnisse_Widerstand_FE2.n_iteration_vec                     

V_EAS=sqrt((2)./(ISA.rho_0*CA_j)*Ergebnisse_stat_Flaechenbelastung.Fleachenbelastung);
Machzahl_j=(V_EAS/sqrt(roh_zu_roh0))/a(hoehe_CR);     %Machzahl Formel 5

%%Schub Gewichts verhaeltnis ausrechnen über die hoehe

Schubkennfeld_S_zu_S0_KF_j=specs.Drosselgrad(2)*(roh_zu_roh0)*exp(-0.35*Machzahl_j*Druckdifferenz*sqrt(specs.bypass));
S_zu_S0_E=1-(1.3+0.25*specs.bypass)*roh_zu_roh0;        
%Standschubverhaeltnis
S_zu_S0_j=Schubkennfeld_S_zu_S0_KF_j*S_zu_S0_E                          %Schub zu Startschub verhältnis


