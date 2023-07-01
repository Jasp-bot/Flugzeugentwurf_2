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


%OD Polare
%Variabeln
G=Ergebnisse_Massen_FE2.M_TO*specs.g;       %Max Abfluggewicht * g (9.80665)
Flughoehe_CR = specs.flight_level * 10^2 ;     % in ft
hoehe_CR = round(unitsratio('m','ft')*Flughoehe_CR);    %umrechnung in m
T=ISA.T;    %Temperatur
a=ISA.a;    %MachzahlISA
roh_zu_roh0=ISA.rho(hoehe_CR)/ISA.rho_0;    %Dichteverhältnis
CA_j=CaVerteilungFE12;           %%FALSCH

V_EAS=sqrt((2*G)./(ISA.rho_0*CA_j*Ergebnisse_stat_Flaechenbelastung.F));
Machzahl_j=V_EAS/sqrt(roh_zu_roh0);     %Machzahl Formel 5



