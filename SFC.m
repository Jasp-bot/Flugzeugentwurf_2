function [sfc_daNh, sfc_1PERh, sfc_1PERs] = SFC (Height,MachNumber,BPR)
%% SPECIFIC FUEL CONSUMPTION TEMPLATE
% Vorlage für die Berechnung des spezifischen Kraftstoffverbrauches
%
%   author:     Andreas Gobbin / edited by Tobias Beelitz
%   date:       19.05.2022
%
%   literature: Torenbeek/ FE Skript/ ISA
%   
%   INPUT:      - Height:       Flughöhe in m
%               - MachNumber:   Flug-Machzahl ohne Einheit
%               - BPR:          Bypass-Ratio ohne Einheit
%
%   OUTPUT:     - sfc_daNH:     SFC in kg/(daNh)
%               - sfc_1PERh:    SFC in 1/h
%               - sfc_1PERs:    SFC in 1/s
%
% clc
% clear all
% close all
load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat;

% Flughoehe = specs.flight_level * 10^2 ;                         % in ft
% hoehe_CR = round(distdim(Flughoehe ,'ft','m'));     % in m
% hoehe_CL = round(distdim(Flughoehe*(2/3) ,'ft','m'));     % in m
% 
% Height = hoehe_CL;
% MachNumber = specs.Ma_CR * 2/3;
% BPR = specs.bypass;

%% Constants and Variables to be defined 
    % Gravity constant in m/s²
    g0 = 9.80665; 

    % Overall engine pressure ratio, from engine type certificate
    OAPR = 60;

    % Turbine entrance temperature, from engine type cert. or from
    % literature 
    TET = 1800;        % [K]

%% Atmospheric Parameters
    % ISA plus condition at Mean Sea Level, from Project Specification
    ISAPlus = 15;
    % atmosdata need to be calculated from your own FE1 atmosphere function (0. Parameterstudie)
    %       (my function is called "f_atmosphere" and returns all important
    %       atmospheric values taking into account the predefined ISAPlus
    %       Offset)
    %           - a:                sonicspeed at Height in m/s
    %           - p:                pressure at Height in kg/m²
    %           - rho:              density at Height in kg/³
    %           - T:                temperature at Height in K
    %           - relativeRho:      (rho/rho_0) 
    %           - relativePressure: (p/p_0)
   % [a, p, rho, T, relativeRho, relativePressure] = f_atmosphere(Height, ISAPlus);    
    
    p = ISA.p(Height);
    a = ISA.a(Height);
    rho = ISA.rho(Height);
    T = ISA.T(Height);
    relativeRho = ISA.rel_rho(Height);
    relativePressure = p/ISA.p0;

    % static temperature for gas generator calculations
    t0  = 288.15;      % temperature at MSL in K
    R = 287.05287;     % gas constant in J/(mol * K)
    tH = p / R / rho;  % static temperature in K

%% Engine Parameters
    % Inlet pressure loss 
    eta_d = 1-(1.3+0.25*BPR)*0.02;

    % Gas generator inlet efficiency
    eta_i = 1-0.7*MachNumber^2*(1-eta_d)/(1+0.2*MachNumber^2);

    % Compressor efficiency (0.84-0.86)
    eta_c = 0.87;

    % Turbine efficiency (0.87-0.89)
    eta_t = 0.90;

    % Nozzle efficiency (0.96-0.98)
    eta_n = 0.98;

    % Fan efficiency (0.85-0.87 Takeoff; 0.82-0.85 Cruise)
    if MachNumber<0.3
        eta_f = 0.89;
    else
        eta_f = 0.87;
    end
    % Gas generator function and parameters
    kappa = 1.4;
    THETA = 1+(kappa-1)*0.5*MachNumber^2;
    PHI   = TET/tH;
    XSI   = THETA*(OAPR^((kappa-1)/kappa)-1);
    G     = (PHI-XSI/eta_c)*(1-1.01/(eta_i^((kappa-1)/kappa)*(XSI+THETA)*...
            (1-XSI/(PHI*eta_c*eta_t))));

%% Specific Fuel Consumption 
    % SFC in kg/(daNh)
    sfc_daNh = 0.697 * sqrt(tH/t0) * (PHI-THETA-XSI/eta_c)/(sqrt(5 * eta_n * ...
                (1 + eta_f * eta_t * BPR) * (G + 0.2 * MachNumber^2 * BPR * eta_d/(eta_f * eta_t)))...
                    - MachNumber * (1 + BPR));
    % SFC in 1/h
    sfc_1PERh = sfc_daNh * g0/10;

    % SFC in 1/s
    sfc_1PERs = sfc_1PERh / 3600;

%end
