%% PS4 detaillierte Widerstandsabschaetzung nach Diederich
% Aufgeteilt in die Kapitel:
    % Fluegelwiderstand:    - Profilwiderstand
    %                       - Ind. Widerstand
    %                       - Transsonischer Widerstand
    % Rumpfwiderstand :     - Rumpfwiderstand
    % Triebwerkswiderstand
    % Leitwerkswiderstand:  - Profilwiderstand
    %                       - Trimmwiderstand
    %                       - Zusatzwiderstand
    % Interferenzwiderstand
clc
clear all
close all



%% Laden von Werten
load Projekt_specs.mat;
load Ergebnisse_Massen_FE2.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Auftrieb_Momente.mat;
load Ergebnisse_Leitwerke.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;


%% Allgemeine Variabelen

Flughoehe_CR = specs.flight_level * 10^2 ;     % in ft
hoehe_CR = round(unitsratio('m','ft')*Flughoehe_CR);



%% Getroffene Annahmen um Rechnungen vor berechnung der richtigen Werte durchfuehren zu können
        % es fehlen Werte als PS2 / PS3
    
    % Profilwiderstand des Fluegels  
phi_50 = tan((NP.versatz_HK + 0.5*DT.l_i_I - 0.5*DT.l_a)/(Ergebnisse_Fluegel.b/2 - specs.R_rumpf));          % Muss noch erfragen wie phi_50 brterchnet werden soll
v_air = specs.Ma_CR * ISA.a(hoehe_CR);
xu_l = 0.035;    % Wert zwischen 0.02 und 0.05
x_u = Ergebnisse_Fluegel.Fluegeltiefen_eta * xu_l;      % Annahme, dass l die Fluegeltiefen an der jeweiligen Fosition auf dem Fluegel sind
    % induzierter Widerstand des Flügels
c_A_F = Ergebnisse_stat_Flaechenbelastung.C_A_CR;       % aus PS4 Formel 11


%% ------------------------------------------------------------------------
%  -------------Fluegelwiderstand nach Diederich---------------------------
%  ------------------------------------------------------------------------


% Profilwiderstand des Fluegels

% PS4 S.2, Formel 4 
k = 0.27 * specs.d_l + 100 * (specs.d_l)^4;

%PS4 S.3, Formel 7
Re_u = x_u * v_air / ISA.kin_visk(hoehe_CR);
Re = Ergebnisse_Fluegel.Fluegeltiefen_eta * v_air / ISA.kin_visk(hoehe_CR);     % Annahme pleas confirm

% PS4 S.2, Formel 6
c_f_la_xu = 1.328./(sqrt(Re_u));
c_f_tu_xu = 0.455./(log(Re_u).^(2.58));
c_f_tu_l = 0.455./(log(Re).^(2.58));
% PS4 S.2, Formel 5
c_f = c_f_tu_l  - xu_l * (c_f_tu_xu - c_f_la_xu);

% PS4 S.2, Formel 3 
c_w_p_min_Re = 2 * c_f * (1 + k * cos(phi_50)^2);

% PS4 S.2, Formel 2 
c_w_p_eta = c_w_p_min_Re + 0.03 * (Ergebnisse_Auftriebsverteilung.c_a_eta).^6;




% PS4 S.3, Formel 10 %%% Sieht alles noch sehr inkorekt aus

test_integ = @(eta) c_w_p_eta .*  (Ergebnisse_Fluegel.Fluegeltiefen_eta) ./ (GRA.l_m);
c_w_p = integral(test_integ, 0, 1, 1001);
test_trapz = trapz(c_w_p_eta .*  (Ergebnisse_Fluegel.Fluegeltiefen_eta) ./ (GRA.l_m))
test_sum = sum(c_w_p)
% schnelltest plot
plot(c_w_p,0:0.001:1)




% Induzierter Widerstand des Fluegels

% PS4 S.4, Formel 12
 c0 = VWA.c_AF_anstieg^2 * (0.0088 * Ergebnisse_Fluegel.lambda - 0.0051 * Ergebnisse_Fluegel.lambda^2) * (1 - 0.0006 * Ergebnisse_Fluegel.streckung_phi25_max^2);

% PS4 S.4, Formel 13
% c1 = GRA.c_a_anstieg * (0.0134 * (Ergebnisse_Fluegel.lambda - 0.3) - 0.0037 * Ergebnisse_Fluegel.lambda^2);
c1 = VWA.c_AF_anstieg * (0.0134 * (Ergebnisse_Fluegel.lambda - 0.3) - 0.0037 * Ergebnisse_Fluegel.lambda^2); %% nicht sicer mit welchem Auftriebsanstieg gerechnet werden muss
% PS4 S.4, Formel 15
tau = 1 - Ergebnisse_Fluegel.streckung_phi25_max * (0.002 + 0.0084 * (Ergebnisse_Fluegel.lambda-0.2)^2);

% PS4 S.4, Formel 14
c2 = (1/tau) * (1 + (5*10^(-6)) * (rad2deg(Ergebnisse_Fluegel.phi_25_max))^3); 

% % delta eps keine ahnung wie das berechnet werden soll
% delta_eps = abs()
% 
% % PS4 S.4, Formel 11
% c_w_ind = c2 * (c_A_F^2)/(pi * Ergebnisse_Fluegel.streckung_phi25_max) +...
%     c1 * c_A_F * delta_eps + c0 * delta_eps^2;
% 


