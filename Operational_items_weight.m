clc
clear 
close
load Projekt_specs.mat

W_DeliveryEmpty = 150000;

%% Instrumente S.289
W_ieg = 0.347 * W_DeliveryEmpty^(5/9) * 9352.6;


%% Hydraulic / Pneumatic S.290
W_hyd = 0.015 * W_DeliveryEmpty + 272;


%% Flight Deck Accomodation S.291
W_flight_deck_acc = 16.5 * W_DeliveryEmpty;

%% PAX Cabin Accomodities S.291
% Pax Sitze gemäß Torenbeek Table 3-2 S.76

seats = [[10.9 21.3 29.9]; [13.6 25.4 35.4]; [21.3 31.8 0]]
% Spalten Eco bis Firts- Jewils single bis triple Seat

% ATT Seats
n_flight_att = 8;           %% Bei All Eco ist es 9
W_Flight_att = 8.2 * 8;     

% PAX Seats                 %% Bei All ECO ist es 432 
W_PAX_seats_ECO = (432 / 3) * seats(1,3); 
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
W_floor = 1.25 * (5.85 * 57.34);

% Sound Proofing and Insulation
V_pax = 5.85 * 57.34 * 1.676;
V_ch  = 1.68 * 4.15 * 57.34;
W_soundProof = 6.17 * (V_pax + V_ch)^(1.07);

%% Cargo Accomodations S.291
% Cargo Retraint + Handling Provisions
W_cargo_restraint = 1.28 * V_ch;
W_cargo_palletprov = 13.67 * V_ch;

%% Standard Emergency Equipment S.291
% Fixed Oxygen System
W_oxygenSys = 18.1 + 1.09 * 432; % Pax 

% Fire Detection
W_to = 300000;
W_fireDetect = 0.0012 * W_to;

% Evacuation
W_evac = 0.453 * 432;