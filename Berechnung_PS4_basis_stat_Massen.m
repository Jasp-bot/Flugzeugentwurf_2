function Berechnung_PS4_basis_stat_Massen(Eingabewert_Iteration)

if Eingabewert_Iteration == 0
    % Statistische Massenberechnung nonlinear regressionsanalyse
    % clc
    % clear all
    % close all
    
    %% Input Excel
    
    AC_Data.aircraftData = importfile_org('aircraftData_org.xlsx');
    % aircraftData(:,1) = type Aircraft				/
    % aircraftData(:,2) = entry int service			/
    % aircraftData(:,3) = Pax designpoint				pax_design
    % aircraftData(:,4) = wing loading [N/m^2]			pressure_wing
    % aircraftData(:,5) = wing area [m^2]				A_wing
    % aircraftData(:,6) = m_ToM [kg]				    m_To
    % aircraftData(:,7) = max Landing Mass [kg]			m_La
    % aircraftData(:,8) = max zero fuel mass [kg]		m_max_zfuel
    % aircraftData(:,9) = operating empty mass [kg]		m_OE
    % aircraftData(:,10) = cruise speed [km/h]			V_cruise
    % aircraftData(:,11) = phys range max PAYLOAD [km]	R_1
    % aircraftData(:,12) = max payload mass [kg]		m_p1
    % aircraftData(:,13) = phys range max FUEL [km]		R_2
    % aircraftData(:,14) = payload mass max fuel [kg]	m_p2
    
    
    AC_Data.pax_ACdata = AC_Data.aircraftData(:,3);      % Pax AC Data in 1				
    AC_Data.pressure_wing = AC_Data.aircraftData(:,4);   % wing loading in N/m^2			
    AC_Data.A_wing = AC_Data.aircraftData(:,5);          % wing area [m^2]				
    AC_Data.m_To = AC_Data.aircraftData(:,6);            %  maximum take off mass m_ToM [kg]
    AC_Data.m_La_max = AC_Data.aircraftData(:,7);        % max landing mass
    AC_Data.m_max_zfuel = AC_Data.aircraftData(:,8);     % max zero fuel mass [kg]
    AC_Data.m_OE = AC_Data.aircraftData(:,9);            % operating empty mass [kg]	
    AC_Data.V_cruise = AC_Data.aircraftData(:,10);       % cruise speed [km/h]		
    AC_Data.R_1 = AC_Data.aircraftData(:,11) * 1000;     % phys range max PAYLOAD [m]	
    AC_Data.m_p1 = AC_Data.aircraftData(:,12);           % max payload mass [kg]		
    AC_Data.R_2 = AC_Data.aircraftData(:,13) * 1000;     % phys range max FUEL [m]
    AC_Data.m_p2 = AC_Data.aircraftData(:,14);           % payload mass max fuel [kg]	
    
    %% Import Daten
    
    load Projekt_specs.mat;              % Projekt spezofikationen im Struct [specs]
    
    
    %% ab hier beginnt die richtige Rechnung Parameterstudie Statistische Massenbestimmung
    %% Berechnungen aller Parameter
    
    
    % a.) nue (zur bestimmung des sfc)
    
    % Transportguete
    % nue_12 = R_12/ (log((m_To)/(m_OE + m_p12))); 
    
    basis_stat_m.nue_1 = AC_Data.R_1 ./ (log((AC_Data.m_To) ...
        ./ (AC_Data.m_OE + AC_Data.m_p1))); % [m]
    
    basis_stat_m.nue_2 = AC_Data.R_2 ./ (log((AC_Data.m_To) ...
        ./ (AC_Data.m_OE + AC_Data.m_p2))); % [m]
    
    
    
    % b.) Bestimmung der Reichweite für Ausweichflug und Warteschleife
    % v_HLD und t_HLD us Aufgabenstelung
    % Konvertieren der Einheiten
    % Vergleiche Projektspezifikationen [specs]
    
    
    % c.) Reservereichweite (masse)
    % Vergleiche Projektspezifikationen [specs]
    % berechnung des Reservetreibstoffs
    % m_fuel_res = m_RF_ALT + m_RF_HLD = (m_OE + m_p12) *...
    % exp((R_HLD+R_ALT)/(nue_12)) - (m_OE + m_p12)
    
    basis_stat_m.m_fuel_res_1 = ((AC_Data.m_OE + AC_Data.m_p1) ...
        .* exp((specs.R_HLD + specs.R_ALT)./(basis_stat_m.nue_1))) - (AC_Data.m_OE + AC_Data.m_p1);
    
    basis_stat_m.m_fuel_res_2 = ((AC_Data.m_OE + AC_Data.m_p2)...
        .* exp((specs.R_HLD + specs.R_ALT)./(basis_stat_m.nue_2))) - (AC_Data.m_OE + AC_Data.m_p2);
    
    % d.) Bestimmung der opperationellen Reichweite
    
    % m_To = m_OE + m_p12 + m_FT + f_RFCR * m_FT + m_fuel_res_12
    % Reserveaufschlag Reisekraftstoff
    basis_stat_m.f_RFCR = 0.05;                % 5 Prozent von m_FT
    
    % Berechnung masse tripfuel
    
    basis_stat_m.m_FT1 = (AC_Data.m_To - (AC_Data.m_OE + AC_Data.m_p1)...
        .* exp((specs.R_HLD + specs.R_ALT) ./ (basis_stat_m.nue_1))) ./ (1 + basis_stat_m.f_RFCR); 
    
    basis_stat_m.m_FT2 = (AC_Data.m_To - (AC_Data.m_OE + AC_Data.m_p2)...
        .* exp((specs.R_HLD + specs.R_ALT) ./ (basis_stat_m.nue_2))) ./ (1 + basis_stat_m.f_RFCR);
    
    % op-Range
    
    basis_stat_m.R_OPS1 = basis_stat_m.nue_1 .* log(AC_Data.m_To ./ (AC_Data.m_To - basis_stat_m.m_FT1));
    basis_stat_m.R_OPS2 = basis_stat_m.nue_2 .* log(AC_Data.m_To ./ (AC_Data.m_To - basis_stat_m.m_FT2));
    
    
    %% Berechnung der Mittelwerte
    
    basis_stat_m.R_OPS_mittel = ((basis_stat_m.R_OPS1 + basis_stat_m.R_OPS2)/2)/1000;  % [km]
    basis_stat_m.m_p_DP = 0.5 * (AC_Data.m_p1 + AC_Data.m_p2);
    basis_stat_m.m_FT_DP = 0.5 * (basis_stat_m.m_FT1 + basis_stat_m.m_FT2);
    basis_stat_m.m_fuel_res_DP = 0.5 * (basis_stat_m.m_fuel_res_1 + basis_stat_m.m_fuel_res_2);
    
    % Berechnung von gamma
    basis_stat_m.gamma_DP = basis_stat_m.m_p_DP./ AC_Data.m_To;
    
    
    % Berechnug von Kappa
    basis_stat_m.kappa_DP = (basis_stat_m.m_FT_DP + basis_stat_m.m_fuel_res_DP) ./ AC_Data.m_To;
    
    basis_stat_m.f_gamma = fit(basis_stat_m.R_OPS_mittel, basis_stat_m.gamma_DP, 'power1');
    basis_stat_m.f_kappa = fit(basis_stat_m.R_OPS_mittel ,basis_stat_m.kappa_DP, 'power1');
    
    
    
    %% Berechnung für Massenabschätzung eigenes FZ
    
    % x_test = max_range_km; 
    % a_test = 0.8758;
    % b_test = -0.1867;
    %  
    % y_test = a_test*x_test^b_test;
    
    % Berechnung gamma Designpunkt eigenes Flugzeug
    coeffizienten_gamma = coeffvalues(basis_stat_m.f_gamma); % ermoeglicht die koeffizienten der POWER Function aus fit zu extrahieren
    A_gamma = coeffizienten_gamma(1,1);
    B_gamma = coeffizienten_gamma(1,2);
    
    Ergebnis_basis_m.gamma_DP = A_gamma .* specs.max_range_basis_km .^ B_gamma;
    
    % Berechnung kappa Designpunkt eigenes Flugzeug
    
    coeffizienten_kappa = coeffvalues(basis_stat_m.f_kappa);
    A_kappa = coeffizienten_kappa(1,1);
    B_kappa= coeffizienten_kappa(1,2);
    
    Ergebnis_basis_m.kappa_DP = A_kappa .* specs.max_range_basis_km .^ B_kappa;
    
    % Berechnung eigenes m_To fuer Flugzeug
    
    
    Ergebnis_basis_m.m_payload = specs.m_pax + specs.m_cargo;
    
    Ergebnis_basis_m.m_To = Ergebnis_basis_m.m_payload / Ergebnis_basis_m.gamma_DP;
    
    % Berechnung der Treibstoffmasse
    
    Ergebnis_basis_m.m_fuel = Ergebnis_basis_m.kappa_DP * Ergebnis_basis_m.m_To;
    
    % Berechnung der Opperativen Leermasse
    
    Ergebnis_basis_m.m_OE = Ergebnis_basis_m.m_To - Ergebnis_basis_m.m_fuel - Ergebnis_basis_m.m_payload;
    
    Ergebnis_basis_m.nue_FZ = specs.max_range_basis_km * 1000 ...
        / (log(Ergebnis_basis_m.m_To / (Ergebnis_basis_m.m_OE + Ergebnis_basis_m.m_payload))); % [m]
    
    % berechnung reservekraftstoff
    Ergebnis_basis_m.m_fuel_res_basis = ((Ergebnis_basis_m.m_OE + Ergebnis_basis_m.m_payload)...
        .* exp((specs.R_HLD + specs.R_ALT)./(Ergebnis_basis_m.nue_FZ)))...
        - (Ergebnis_basis_m.m_OE + Ergebnis_basis_m.m_payload);
    
    Ergebnis_basis_m.m_TF_basis = (Ergebnis_basis_m.m_To - (Ergebnis_basis_m.m_OE + Ergebnis_basis_m.m_payload)...
        .* exp((specs.R_HLD + specs.R_ALT) ./ (Ergebnis_basis_m.nue_FZ))) ./ (1 +basis_stat_m.f_RFCR); 
    
    %% Berechnung Kappa_Shrink
    
    Ergebnis_basis_m.kappa_DP_Shrink = A_kappa .* specs.max_range_Shrink_km .^ B_kappa;
    Ergebnis_basis_m.gamma_DP_Shrink = A_gamma .* specs.max_range_Shrink_km .^ B_gamma;

    save Zwischenergebnisse_PS4_basis_stat_Massen.mat AC_Data basis_stat_m;
%% FE2

elseif Eingabewert_Iteration == 1
    load Ergebnisse_Massen_FE2.mat;
    load Ergebnisse_Basis_stat_m.mat

    Ergebnis_basis_m.gamma_DP = (specs.m_pax + specs.m_cargo)./ Ergebnisse_Massen_FE2.M_TO;
    Ergebnis_basis_m.kapap_DP = Ergebnisse_Massen_FE2.M_F / Ergebnisse_Massen_FE2.M_TO;
    Ergebnis_basis_m.m_payload = (specs.m_pax + specs.m_cargo);
    Ergebnis_basis_m.m_To = Ergebnisse_Massen_FE2.M_TO;
    Ergebnis_basis_m.m_fuel = Ergebnisse_Massen_FE2.M_F;
    Ergebnis_basis_m.m_OE = Ergebnisse_Massen_FE2.M_OE;
    Ergebnis_basis_m.nue_FZ = Ergebnis_basis_m.nue_FZ; %% kein plan
    Ergebnis_basis_m.m_fuel_res_basis = Ergebnisse_Massen_FE2.M_F_ALT;
    Ergebnis_basis_m.m_TF_basis = Ergebnisse_Massen_FE2.M_TO - Ergebnisse_Massen_FE2.M_Z_Tripfuel;
    Ergebnis_basis_m.kappa_DP_Shrink = Ergebnis_basis_m.kappa_DP_Shrink
    Ergebnis_basis_m.gamma_DP_Shrink = Ergebnis_basis_m.gamma_DP_Shrink

end

%% Safe File Ergebnisse in .mat

save Ergebnisse_Basis_stat_m.mat Ergebnis_basis_m;






