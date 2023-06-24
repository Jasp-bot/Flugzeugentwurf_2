%% Projekt Spezifikationen

specs.g = 9.80665;                                  % [m/(s^2)]

% Reichweite
specs.max_range_basis_nm = 5050;
specs.max_range_basis_km = unitsratio('km','nm')*specs.max_range_basis_nm;     % umrechnen von nm in km
specs.max_range_Shrink_nm = 9200;
specs.max_range_Shrink_km = unitsratio('km','nm')*specs.max_range_Shrink_nm;


% massenabschaetzung Ausweich-Reichweite
specs.v_HLD = convvel(250, 'kts','m/s') ;           % [m/s]
specs.t_HLD = 45*60;                                % [s] 
specs.R_HLD = (specs.v_HLD * specs.t_HLD);          % [m]
specs.R_HLD_km = specs.R_HLD / 1000;                % [km]

specs.R_ALT = unitsratio('m','nm')*200;               % [m]
specs.R_ALT_km = specs.R_ALT / 1000;                % [km]

specs.H_HLD = unitsratio('m', 'ft')*1500;
specs.flight_level_ALT = 200;


% Kabine und Massenanteile
specs.n_pax = 311;                                  % Paxanzahl
specs.m_pax = 97 * specs.n_pax;                     % m_pax = (Gewicht Passagier(80kg)) + (Gewicht gepäck des Passagiers(17kg))
specs.m_cargo = 14000;                              % Masse Cargo in[kg]
specs.m_bag = 17;                                   % Masse Cargo pro pax [kg]
specs.n_pax_shrink = 249;                           % Paxanzahl shrink
specs.n_flight = 2;                                 % Anzahl Piloten
specs.n_crew = 9;                                   % Anzahl Crew
specs.n_pax_all_eco = 432;                          % Anzahl passagiere All eco
specs.m_pax_all_eco = 97 * specs.n_pax_all_eco;     % Gewicht in kg
specs.n_Toilette = 8;                               % Anzahl Toiletten



specs.delta_m_p = -0.2;                             % Verringerung von Basis nach shrink (20% kleiner) [1]
specs.k_vOE = 0.5;                                  % Rüstmassenfaktor [1]

% Aerodynamik

specs.kappa_luft = 1.4;
specs.Ma_CR = 0.82;
specs.Ma_MO = 0.84;
specs.flight_level = 340;

% Triebwerksdaten 

specs.bypass = 10;                                  % Bypassverhaeltnis 
specs.n_TW = 2;      % Anzahl Triebwerke
m_TW = 7549;   % Quelle: https://www.easa.europa.eu/sites/default/files/dfu/EASA%20E111%20TCDS%20issue%2012_Trent%20XWB.pdf
specs.m_TW = m_TW * 0.10 + m_TW; % Masse eines TW in kg bsp terd XWB

% Landespezifikationen

specs.s_LDG_specs = 2600;  % max Landestrecke bis zum Stillstand
specs.s_LDG_safty = 0.6 * specs.s_LDG_specs; % Vorgeschriebene Strecke bis zum Halt des FZ (0.6 Sicherheitswert für error)

% Fluegel/ Geometrie

specs.D_rumpf = 6.21;
specs.R_rumpf = specs.D_rumpf/2;
specs.l_rumpf = 74.92;
specs.h_rumpf = 6.21;
specs.d_l = 0.13;
% Leitwerke

specs.coanlaenge = 15;
specs.HLW_beginn = 6.6;
specs.HLW_d2l = 0.12;

% Massenanteile 

specs.m_APU = 335; % [kg] HGT1700 Auxiliary Power Unit




save Projekt_specs.mat specs;