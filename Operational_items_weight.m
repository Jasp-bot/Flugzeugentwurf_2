clc
clear 
close
load Projekt_specs.mat
load Ergebnisse_Basis_stat_m.mat;

W_DeliveryEmpty = 140000;

%% Instrumente S.289
W_ieg_alt = 0.347 * W_DeliveryEmpty^(5/9) * specs.max_range_basis_km^(1/4);
% ????????????????????????????????????????????????????????????????????
Technologiefaktor_W_ieg = 0.8;
W_ieg = W_ieg_alt * Technologiefaktor_W_ieg;


% Absolut nicht sicher !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % Electronics P.290
% v_cabine = 62.1 * pi * (specs.R_rumpf^2) /2;
% 
% % (8-43)
% P_el_APU = 3.64 * v_cabine^0.7;
% % (8-41)
% W_el_AC = 16.3 * P_el_APU * (1 - 0.033 * sqrt(P_el_APU));
% 
% % (8-40)
% W_el_DC = 0.02 * Ergebnis_basis_m.m_To + 181;


%% Hydraulic / Pneumatic S.290 (8-39)
W_hyd = 0.015 * W_DeliveryEmpty + 272;


%% Flight Deck Accomodation S.291
W_flight_deck_acc = 16.5 * W_DeliveryEmpty^(0.285);

%% PAX Cabin Accomodities S.291

% Pax Sitze gemäß Torenbeek Table 3-2 S.76
seats = [[10.9 21.3 29.9]; [13.6 25.4 35.4]; [21.3 31.8 0]]
% Spalten Eco bis Firts- Jewils single bis triple Seat

% ATT Seats
n_flight_att = 9;           %% Bei All Eco ist es 9
W_Flight_att = 8.2 * n_flight_att;     

% PAX Seats                 %% Bei All ECO ist es 432 
W_PAX_seats_ECO = (specs.n_pax_all_eco / 3) * seats(1,3); 
%W_PAX_seats_Basic =

% Galley/Pantry
Galley_main = 113.4;
Galley_medium = 45.3;
Galley_small = 29.5;

W_Galley = Galley_main * 2 + Galley_medium * 4 + Galley_small * 5;

%Lavatory/ Toiletts
lavatory = 136;
W_lavatory = lavatory * 7 + lavatory * 0.2 * 1; 

%Floor Covering
W_floor = 1.25 * (5.85 * 57.34)^(1.15);

% Sound Proofing and Insulation ??????????????????????????????????????
V_pax = 5.85 * 57.34 * 1.676;
V_ch  = 1.68 * 4.15 * 57.34;
W_soundProof = 6.17 * (V_pax + V_ch)^(1.07);

%% Cargo Accomodations S.291
% Cargo Retraint + Handling Provisions
W_cargo_restraint = 1.28 * V_ch;
W_cargo_palletprov = 13.67 * V_ch;

%% Standard Emergency Equipment S.291
% Fixed Oxygen System
W_oxygenSys = 18.1 + 1.09 * specs.n_pax_all_eco; % Pax 

% Fire Detection
W_to = Ergebnis_basis_m.m_To;
W_fireDetect = 0.0012 * W_to;

% Evacuation
W_evac = 0.453 * specs.n_pax_all_eco;

%% Aircon
% (8-45) P.293
W_aircon = 14 * 57.34^1.28; % Annahme l_pc =57.34m

%% Misscellaneous Weight
   % Wird weggelassen als Technologiefaktor, aufgrund seines geringen
   % Außmaßes von 1% Delivery empty mass


M_Airframe_service_and_equipment = W_hyd + W_flight_deck_acc + W_Flight_att + W_PAX_seats_ECO + W_Galley + W_lavatory + W_floor +...
    W_soundProof + W_cargo_restraint + W_cargo_palletprov + W_oxygenSys + W_fireDetect +W_evac + W_aircon;
