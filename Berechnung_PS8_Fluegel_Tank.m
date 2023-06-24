%% PS8 Fluegelgeometrie und Tankauslwgung

function Berechnung_PS8_Fluegel_Tank

clc
clear all
close all

%% Laden der Werte

load Projekt_specs.mat;
load Ergebnisse_Basis_stat_m.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Start_Landeanforderungen.mat;
load Ergebnisse_Shrink_massen.mat;


% Deklaration
Flughoehe = specs.flight_level * 10^2 ;                         % in ft
hoehe =round(unitsratio('m','ft')*Flughoehe);                     % in m


p_CR = ISA.p(hoehe,1);
rho_CR = ISA.rho(hoehe,1);
a_CR = ISA.a(hoehe,1);                              % Schallgeschw. in CR-hoehe


a_MO = a_CR;                                      % Schallgeschw. in CR-hoehe = Schallgeschw. max Opperating



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Berechnungen %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Buffet Onset  (BO)

BO.v_MO = specs.Ma_MO * a_MO; % max opperating speed
BO.k_R = 0.98;

% Auftriebsbeiwert bei Ma_MO und CR altitude
BO.c_A_ICA_MO = 2/(rho_CR * BO.v_MO^2) * ( Ergebnisse_stat_Flaechenbelastung.G_To / Ergebnisse_stat_Flaechenbelastung.F) * BO.k_R; %

% Berechnungen für Auftriebsbeiwert buffet onset
BO.k_trim = (1.05 : 0.01 : 1.1).';
BO.k_VA = (1.1 : 0.01 : 1.2);
BO.k_n = 1.3;

BO.c_A_Bo = BO.c_A_ICA_MO .* BO.k_n .* BO.k_trim .* BO.k_VA; % Matrix x =  k_trim , y= k_VA => 6x11 matrix
BO.c_A_Bo_plot = BO.c_A_Bo(1,1);

% Berechnung Ma_BO

BO.wert = linspace(0.65, 0.85);        % Wert in dem sich ca Ma_BO bewegt
BO.f = @(x) (1459.9 .* x.^4 - 4203.5 .* x.^3 + 4483.4 .* x.^2 - 2100.9 .*x + 366.16);      % Funktion zu erstellung der BO-Kurve

BO.vec_1 = [BO.wert; BO.f(BO.wert)];
BO.vec_2 = [BO.wert; repmat(BO.c_A_Bo_plot, 1, 100)];          % repmat füllt eine MAtrix mit werten repmat(wiederholter wert, x, y)
schnittpunkt_BO_Diagramm = InterX(BO.vec_1, BO.vec_2);

BO.Ma_BO = schnittpunkt_BO_Diagramm(1,1);


%% Pfeilung

% unter der annahme, dass Ma_eff = Ma_Bo und Ma = Ma_MO
%phi_25 = acos((Ma_eff/Ma)^2);
Ergebnisse_Fluegel.phi_25_max = (acos((BO.Ma_BO / specs.Ma_MO)^2)); 

%% Streckung

% Formel aus Anhang PS8 "Pitch up limit"
Ergebnisse_Fluegel.f_streckung = @(x) (-8.892 .* log(x) + 40.285);
Ergebnisse_Fluegel.streckung_phi25_max = Ergebnisse_Fluegel.f_streckung(rad2deg(Ergebnisse_Fluegel.phi_25_max));

%% Zuspitzung

% lambda soll mit 0.3 bis 0.4 angenommen werden

Ergebnisse_Fluegel.lambda = 0.3;

%% Spannweite
% b = sqrt(F*strechung_phi_25_max)
Ergebnisse_Fluegel.b = sqrt(Ergebnisse_stat_Flaechenbelastung.F * Ergebnisse_Fluegel.streckung_phi25_max);


%% Vorderkantenpfeilung
% GL4 PS8
% (0-0.25) kommt von den koordinaten (x/l)_vk - (x/l)_25
Ergebnisse_Fluegel.phi_VK_max = atan( tan(Ergebnisse_Fluegel.phi_25_max) - (4/Ergebnisse_Fluegel.streckung_phi25_max)*(-0.25)...
    * ((1-Ergebnisse_Fluegel.lambda)/(1+Ergebnisse_Fluegel.lambda)));
Ergebnisse_Fluegel.phi_VK_max_grad = rad2deg(Ergebnisse_Fluegel.phi_VK_max);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fluegelgeometrie %%%%%%%%%%%%%%%%%%%%%
%% Einfach-Trapezfluegel (ET)

ET.l_i = (2 * Ergebnisse_stat_Flaechenbelastung.F) / (Ergebnisse_Fluegel.b * (Ergebnisse_Fluegel.lambda+1)); % Fluegelbreite innen
ET.l_a = Ergebnisse_Fluegel.lambda * ET.l_i; % Fluegelbreite außen

% Berechnung feur außentrapez

ET.phi_HK_a = atan( tan(Ergebnisse_Fluegel.phi_25_max) - (4/Ergebnisse_Fluegel.streckung_phi25_max)*(-0.25+1)...
    * ((1-Ergebnisse_Fluegel.lambda)/(1+Ergebnisse_Fluegel.lambda)));
ET.phi_HK_a_grad = rad2deg(ET.phi_HK_a); 

ET.A_1 = (((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2)*(Ergebnisse_Fluegel.b/2))/2) + (Ergebnisse_Fluegel.b/2) * ET.l_a;
ET.A_2 = ((tan(ET.phi_HK_a)*Ergebnisse_Fluegel.b/2)*(Ergebnisse_Fluegel.b/2))/2;
ET.A_trapezfluegel = ET.A_1 - ET.A_2;
ET.A_Fluegel_halb = Ergebnisse_stat_Flaechenbelastung.F/2;

%% Doppel-Trapezfluegel

DT.l_a = ET.l_a;

% Itteration zur Bestimmung des Hinterkantenwinkels und der Flügelflaeche

phi_HK_dt_low = deg2rad(ET.phi_HK_a_grad);
phi_HK_dt_high = deg2rad(ET.phi_HK_a_grad + 15);
dF = 1;

while abs(dF) > 0.0001
%for zaehlvar = 1:1000
    phi_HK_dt_try = (phi_HK_dt_low + phi_HK_dt_high)/(2);

   % phi_HK_dt =deg2rad(25.8573) ; % variabele in grad 

    y_pos_kink = 0.25:0.01:0.35;% 0.3
    y_pos_kink_vector = 6; % y_pos_kink_vector = wert;  wert von 1 bis 11 für Ansatzpunkt des Kinks

% Berechnung von laengenabschnitten von b/2

    DT.s_A = (Ergebnisse_Fluegel.b/2) -((Ergebnisse_Fluegel.b/2) .* y_pos_kink(1,y_pos_kink_vector));
    DT.s_I = ((Ergebnisse_Fluegel.b/2) .* y_pos_kink(1,y_pos_kink_vector))- specs.R_rumpf;
    DT.s_R = specs.R_rumpf;

% Berechnung F_A
    
    DT.F_A = DT.l_a .* DT.s_A + (((tan(Ergebnisse_Fluegel.phi_VK_max).* DT.s_A).* DT.s_A)/2)...
        - (((tan(phi_HK_dt_try).* DT.s_A).* DT.s_A)/2);

    DT.l_i_A = (tan(Ergebnisse_Fluegel.phi_VK_max).* DT.s_A) + DT.l_a - (tan(phi_HK_dt_try).* DT.s_A);

% Berechnung F_I
    DT.l_a_I = DT.l_i_A;
    DT.F_I = DT.l_a_I .* DT.s_I + (((tan(Ergebnisse_Fluegel.phi_VK_max).* DT.s_I).* DT.s_I)/2);

    DT.l_i_I = DT.l_a_I + (tan(Ergebnisse_Fluegel.phi_VK_max).* DT.s_I);

% berechnungen F_R

    DT.l_i_R = DT.l_i_I;

    DT.F_R = specs.R_rumpf .* DT.l_i_R;

% F_ges 

% Ergebniss_fleache = solve((l_a_dt .* s_A + (((tan(phi_VK_max).*s_A).*s_A)/2) - (((tan(phi_HK_dt).*s_A).*s_A)/2)) + (l_a_I .* s_I + (((tan(phi_VK_max).*s_I).*s_I)/2)) + F_R == Fluegelflaeche_FZ/2, phi_HK_dt)

    DT.F_check = DT.F_A + DT.F_I + DT.F_R;
    dF = (Ergebnisse_stat_Flaechenbelastung.F /2) - DT.F_check;
    
    if dF < 0
        phi_HK_dt_low = phi_HK_dt_try;
    else
        phi_HK_dt_high = phi_HK_dt_try;
    
    end

end
DT.phi_HK_dt = phi_HK_dt_try;
DT.y_kink = y_pos_kink(1,y_pos_kink_vector);
DT.phi_VK_max = Ergebnisse_Fluegel.phi_VK_max;


% rad2deg(DT.phi_HK_dt)

%% Neutralpunke (Doppeltrapez)


% Aussentrapez
NP.lambda_A = DT.l_a / DT.l_i_A;

NP.x_phi_25_A = (0.75 * DT.l_i_A) + (tan(DT.phi_HK_dt)* DT.s_A) - (DT.l_a * 0.75); % 
NP.phi_25_A = tan(NP.x_phi_25_A / DT.s_A);
NP.streckung_A = 2 * (DT.s_A)^2 / DT.F_A;

NP.l_mue_A = DT.l_i_A * ((2/3) *((1+ NP.lambda_A + NP.lambda_A^2)/(1 + NP.lambda_A)) );
NP.y_SP_A = DT.s_A * ((1/3)* ((1 + 2 * NP.lambda_A)/(1 + NP.lambda_A)));
NP.x_NP_A = DT.l_i_A * ((1/4) + (NP.streckung_A/12) * (1 + NP.lambda_A * 2) * tan(NP.phi_25_A));
% x_NP_A = l_i_A * ((1/4) + (streckung_phi25_max/12) * (1 + lambda_A * 2) * tan(phi_25_max));

NP.v50_A = [0,  ((Ergebnisse_Fluegel.b/2)*(1 - y_pos_kink(1,y_pos_kink_vector)));...
    (DT.l_a/2), ((tan(DT.phi_HK_dt) * DT.s_A) + DT.l_i_A/2)];
NP.y_SP_A_vector = [NP.y_SP_A, NP.y_SP_A;...
                0, 30];

NP.schnitt_A_SP = InterX(NP.v50_A, NP.y_SP_A_vector);
NP.x_SP_A = NP.schnitt_A_SP(2); 


% Innentrapez

NP.lambda_I = DT.l_a_I / DT.l_i_I;

NP.x_phi_25_I = (0.75 * DT.l_i_I)  - (DT.l_a_I * 0.75);
NP.phi_25_I = tan(NP.x_phi_25_I / DT.s_I);
NP.streckung_I = 2*(DT.s_I)^2 / DT.F_I;


NP.l_mue_I = DT.l_i_I * ((2/3) *((1+ NP.lambda_I + NP.lambda_I^2)/(1+ NP.lambda_I)) );
NP.y_SP_I = DT.s_I * ((1/3)* ((1+2* NP.lambda_I)/(1+ NP.lambda_I)));
NP.x_NP_I = DT.l_i_I * ((1/4) + (NP.streckung_I / 12) * (1 + NP.lambda_I * 2) * tan(NP.phi_25_I));
% x_NP_I = l_i_I * ((1/4) + (streckung_phi25_max/12) * (1 + lambda_I * 2) * tan(phi_25_max));

NP.v50_I = [DT.s_A, (Ergebnisse_Fluegel.b/2)- specs.R_rumpf;...
    ((tan(DT.phi_HK_dt) * DT.s_A) + DT.l_i_A/2), ((tan(DT.phi_HK_dt) * DT.s_A)+ DT.l_i_I/2)];
NP.y_SP_I_vector = [DT.s_A+ NP.y_SP_I, DT.s_A+ NP.y_SP_I;...
                0, 50];

NP.schnitt_I_SP = InterX(NP.v50_I, NP.y_SP_I_vector);
NP.x_SP_I = NP.schnitt_I_SP(2);

% Rumpfteil
DT.l_i_R = DT.l_i_I;
DT.l_a_R = DT.l_i_R;
NP.lambda_R = DT.l_a_R / DT.l_i_R;
NP.phi_25_R = deg2rad(0);

NP.streckung_R = 2*(DT.s_R)^2 / DT.F_R;

NP.l_mue_R = DT.l_i_R * ((2/3) *((1+ NP.lambda_R + NP.lambda_R^2)/(1+NP.lambda_R)) );
NP.y_SP_R = DT.s_R * ((1/3)* ((1+2*NP.lambda_R)/(1+NP.lambda_R)));
NP.x_NP_R = DT.l_i_R * ((1/4) + (NP.streckung_R/12) * (1 + NP.lambda_I * 2) * tan(NP.phi_25_R));
% x_NP_R = l_i_R * ((1/4) + (streckung_phi25_max/12) * (1 + lambda_I * 2) * tan(phi_25_max));

NP.x_SP_R = (DT.l_i_R/2)+(tan(DT.phi_HK_dt) * DT.s_A);

%% Gesamtschwerpunkt

NP.versatz_HK = tan(DT.phi_HK_dt)*DT.s_A;

NP.l_mue_ges = (NP.l_mue_A*DT.F_A + NP.l_mue_I*DT.F_I + NP.l_mue_R*DT.F_R)...
    /(DT.F_A+DT.F_I+DT.F_R);

NP.y_SP_ges = ((NP.y_SP_A+specs.R_rumpf+DT.s_I)*DT.F_A ...
    +(NP.y_SP_I+specs.R_rumpf)* DT.F_I + NP.y_SP_R*DT.F_R)/(DT.F_A+DT.F_I+DT.F_R);

NP.x_NP_ges = ((NP.x_NP_A) * DT.F_A + (DT.l_i_R + NP.versatz_HK - NP.x_NP_I)...
    *DT.F_I + (DT.l_i_R + NP.versatz_HK - NP.x_NP_R)*DT.F_R)/((DT.F_A+DT.F_I+DT.F_R));

% x_SP_ges = ((x_SP_A) * F_A + (l_i_R + versatz_HK - x_SP_I)*F_I + (l_i_R + versatz_HK - x_SP_R)*F_R)/((F_A+F_I+F_R));
% x_SP_ges = ((versatz_HK + l_i_R)-(tan(phi_VK_max)*s_I - x_SP_A)*F_A + ((versatz_HK + l_i_R)-x_SP_I)*F_I + ((versatz_HK + l_i_R)-x_SP_R*F_R)/((F_A+F_I+F_R));
NP.x_SP_ges = ((NP.x_SP_A) * DT.F_A + (NP.x_SP_I)*DT.F_I + (NP.x_SP_R)*DT.F_R)/((DT.F_A+DT.F_I+DT.F_R));

%% Tankvolumen

% Außentrapez A

Tank.Fq1_a_A = 0.13 * 0.9 * DT.l_a * (0.65 * DT.l_a - 0.15 * DT.l_a);
Tank.Fq2_i_A = 0.13 * 0.9 * DT.l_i_A * (0.65 * DT.l_i_A- 0.15 * DT.l_i_A);

Tank.a1_A = 0.13 * 0.9 * DT.l_a;
Tank.b1_A = DT.l_a - (DT.l_a * 0.5);

Tank.a2_A = 0.13 * 0.9 * DT.l_i_A;
Tank.b2_A = DT.l_i_A - (DT.l_i_A * 0.5);

Tank.c_aussen = DT.s_A;

Tank.V_OB_A = (Tank.c_aussen/3) * (Tank.Fq1_a_A + Tank.Fq2_i_A + (Tank.a1_A * Tank.b2_A + Tank.a2_A * Tank.b1_A)/(2));


% InnenTrapez I
Tank.Fq1_a_I = 0.13 * 0.9 * DT.l_a_I * (0.65 * DT.l_a_I - 0.15 * DT.l_a_I);
Tank.Fq2_i_I = 0.13 * 0.9 * DT.l_i_I * (0.65 * DT.l_i_I - 0.15 * DT.l_i_I);

Tank.a1_I = 0.13 * 0.9 * DT.l_a_I;
Tank.b1_I = DT.l_a_I - (DT.l_a_I * 0.5);

Tank.a2_I = 0.13 * 0.9 * DT.l_i_I;
Tank.b2_I = DT.l_i_I - (DT.l_i_I * 0.5);

Tank.c_innen = DT.s_I;

Tank.V_OB_I = (Tank.c_innen/3) * (Tank.Fq1_a_I + Tank.Fq2_i_I + (Tank.a1_I * Tank.b2_I + Tank.a2_I * Tank.b1_I)/(2));

% Rumpftrapez

Tank.Fq1_a_R = 0.13 * 0.9 * DT.l_a_R * (0.65 * DT.l_a_R - 0.15 * DT.l_a_R);
 % Fq2_i_R = 0.13 * 0.9 * l_a_R * (0.65 * l_a_R - 0.15 * l_a_R);

Tank.a1_R = 0.13 * 0.9 * DT.l_a_R;
Tank.b1_R = DT.l_a_R - (DT.l_a_R * 0.5);

Tank.a2_R = 0.13 * 0.9 * DT.l_i_R;
Tank.b2_R = DT.l_i_R - (DT.l_i_R * 0.5);

Tank.c_rumpf = specs.R_rumpf; 

Tank.V_OB_R = (Tank.c_rumpf/3) * (Tank.Fq1_a_R * 2 + (Tank.a1_R * Tank.b2_R + Tank.a2_R * Tank.b1_R)/(2));

% Gesamtvolumen Obelist

Tank.V_OB_ges = Tank.V_OB_A + Tank.V_OB_I + Tank.V_OB_R; % fuer shrink
Tank.V_OB_Basis = Tank.V_OB_A + Tank.V_OB_I;
% Tankvolumen
Tank.ka = 0.85;
Tank.kb = 0.95;
Tank.V_Tank = Tank.V_OB_ges * Tank.ka * Tank.kb * 1000; % in Liter
Tank.V_Tank_Basis = Tank.V_OB_Basis * Tank.ka * Tank.kb * 1000; % in Liter
% treibstoffmasse
Tank.rho_kerosin = 0.785; % kg/l

Tank.masse_Treibstoff_max = Tank.rho_kerosin * Tank.V_Tank*2;
Tank.masse_Treibstoff_Basis = Tank.rho_kerosin * Tank.V_Tank_Basis * 2;
Tank.Test_basis_Fuelmass = Ergebnis_basis_m.m_fuel < Tank.masse_Treibstoff_Basis;
Tank.Test_shrink_Fuelmass = Ergebnis_shrink_m.m_fuel < Tank.masse_Treibstoff_max;




%% Flügeltiefen abhaengig von eta(0 bis 1)

% Bedingungen  von eta(0) bis eta(R_rumpf) = l_i_R
% Bedingung von eta(R-Rumpf) bis eta(s_I) = l_i_R - tan(phi_VK)*s_I(eta)
% Bedingung von eta(l_i_A) bis eta(l_a_A) = l_i_A - tan(phi_VK_max)*((((n*10^(-2))*(b/2)) - (R_rumpf + s_I))) + tan(phi_HK_dt)*(((n*10^(-2))*(b/2)) - (R_rumpf + s_I));

for n = 1:1000
   if (n*10^(-3))*(Ergebnisse_Fluegel.b/2)<= specs.R_rumpf
       Ergebnisse_Fluegel.Fluegeltiefen_eta(1,1) = DT.l_i_R; 
       Ergebnisse_Fluegel.Fluegeltiefen_eta(1,n+1) = DT.l_i_R;
       Ergebnisse_Fluegel.zaehlvariabele_eta_Ru = n;
   elseif (n*10^(-3))*(Ergebnisse_Fluegel.b/2)< specs.R_rumpf + DT.s_I
       Ergebnisse_Fluegel.Fluegeltiefen_eta(1,n+1) = DT.l_i_R - tan(Ergebnisse_Fluegel.phi_VK_max)*(((n*10^(-3))*(Ergebnisse_Fluegel.b/2))-specs.R_rumpf);
   else
       Ergebnisse_Fluegel.Fluegeltiefen_eta(1,n+1) = DT.l_i_A - tan(Ergebnisse_Fluegel.phi_VK_max)*((((n*10^(-3))*(Ergebnisse_Fluegel.b/2))...
           - (specs.R_rumpf + DT.s_I))) + tan(DT.phi_HK_dt)*(((n*10^(-3))*(Ergebnisse_Fluegel.b/2)) - (specs.R_rumpf + DT.s_I));
   end
end
Ergebnisse_Fluegel.Fluegeltiefen_eta_oR = Ergebnisse_Fluegel.Fluegeltiefen_eta(1, Ergebnisse_Fluegel.zaehlvariabele_eta_Ru:length(Ergebnisse_Fluegel.Fluegeltiefen_eta));
<<<<<<< Updated upstream

Ergebnisse_Fluegel.F = Ergebnisse_stat_Flaechenbelastung.F;

Ergebnisse_Fluegel.l_mue = (Ergebnisse_Fluegel.b / Ergebnisse_Fluegel.F)* trapz(Ergebnisse_Fluegel.Fluegeltiefen_eta.^2) * 10^(-3); 
Ergebnisse_Fluegel.l_m = trapz(Ergebnisse_Fluegel.Fluegeltiefen_eta)*10^(-3);
=======
Ergebnisse_Fluegel.Fluegeltiefen_eta_Ru = Ergebnisse_Fluegel.Fluegeltiefen_eta(1, 1:Ergebnisse_Fluegel.zaehlvariabele_eta_Ru);
>>>>>>> Stashed changes

%% Safe

save Ergebnisse_Fluegel_Tank_NP.mat Ergebnisse_Fluegel ET DT NP Tank BO


