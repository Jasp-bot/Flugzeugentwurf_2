%% CG Berechnung
clc
clear all
close all

load Ergebnisse_Massen_FE2.mat;
load Projekt_specs.mat;

cg_seats=43.85;

% [Masse, x, y, z];
seats_MAT=[Anteile_einzel_Massen_FE2.Opperational_Items.Seating,0.55*specs.l_rumpf,0,0];%% Z h�he liegt direkt im Rumpfmittelpunkt Y=0 weil symetrisch
cockpit_MAT=[4500,0.05*specs.l_rumpf,0,0];
tail_MAT=[Anteile_einzel_Massen_FE2.Airplane_Structure.Tail_group];



