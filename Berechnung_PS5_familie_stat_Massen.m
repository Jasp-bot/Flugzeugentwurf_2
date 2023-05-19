function Berechnung_PS5_familie_stat_Massen

clc
clear all
close all

%% Inputs

load Projekt_specs.mat;
load Zwischenergebnisse_PS4_basis_stat_Massen.mat;
load Ergebnisse_Basis_stat_m.mat;

%% Berechnung der MAssenanteile Shrink-Version


% Berechnung voruebergehende Cargomasse

Ergebnis_shrink_m.m_cargo_vorl = Ergebnis_basis_m.m_payload - ((0.2)/(1)) * specs.n_pax * specs.m_bag;

% payloadmasse shink

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% nochmal überprüfen lassen
Ergebnis_shrink_m.m_payload = (1+(specs.delta_m_p)/1) * specs.n_pax + specs.m_cargo;      %Ergebnis_shrink_m.m_cargo_vorl;
% nicht 100% sicher ob man hier mit den 14t rechnet oder mit m_cargo_vorl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% m_OE shink

Ergebnis_shrink_m.delta_m_p_shink_kg = Ergebnis_shrink_m.m_payload - Ergebnis_basis_m.m_payload;

Ergebnis_shrink_m.m_OE = Ergebnis_basis_m.m_OE * ...
    (1+(specs.k_vOE * (Ergebnis_shrink_m.delta_m_p_shink_kg)/(Ergebnis_basis_m.m_payload)));

% m_k_shrink

Ergebnis_shrink_m.m_fuel = Ergebnis_basis_m.m_To * Ergebnis_basis_m.kappa_DP_Shrink;

% m_to Shrink

Ergebnis_shrink_m.m_To_vorl = Ergebnis_shrink_m.m_fuel + Ergebnis_shrink_m.m_OE + Ergebnis_shrink_m.m_payload;

% Nach vorgabe wird die ideale Familie erreicht, indem die Fracht im Rumpf
% substituiert wird

Ergebnis_shrink_m.delta_cargo = Ergebnis_shrink_m.m_To_vorl - Ergebnis_basis_m.m_To;
Ergebnis_shrink_m.m_cargo = specs.m_cargo - Ergebnis_shrink_m.delta_cargo;

Ergebnis_shrink_m.m_To = Ergebnis_shrink_m.m_To_vorl - Ergebnis_shrink_m.delta_cargo;

save Ergebnisse_Shrink_massen Ergebnis_shrink_m


