5%% PS4 detaillierte Widerstandsabschaetzung nach Diederich
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
phi_50 = Ergebnisse_Fluegel.phi_50; % tan((NP.versatz_HK + 0.5*DT.l_i_I - 0.5*DT.l_a)/(Ergebnisse_Fluegel.b/2 - specs.R_rumpf));          % Muss noch erfragen wie phi_50 brterchnet werden soll
v_air = specs.Ma_CR * ISA.a(hoehe_CR);
xu_l = 0.035;    % Wert zwischen 0.02 und 0.05
x_u = Ergebnisse_Fluegel.Fluegeltiefen_eta_oR * xu_l;      % Annahme, dass l die Fluegeltiefen an der jeweiligen Fosition auf dem Fluegel sind
    % induzierter Widerstand des Flügels
%c_A_F = linspace(0,1); % Ergebnisse_stat_Flaechenbelastung.C_A_CR;       % aus PS4 Formel 11

    % Transsonischer Widerstand
Ma_unendlich = specs.Ma_CR;

    % Rumpfwiderstand
S_G_Ru = Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.Sg_12;
c_A_alpha_F = GRA.c_a_anstieg;      % Annahme bitte ueberpruefen !!!!!!!!!!!!!!!!!!!!

z_abstand = 3; % Abstand zwischen Profilsehnen angenommen vergleiche Torenbeek s480
zaehlvariabele_itt = 1;

    % Triebwerke
S_G_TW = Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Turbine_Area; % Umspuehlte TW oberflaeche
l_TW = Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Turbine_length;
d_TW = Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Turbine_Diameter;
TW_Einbauwinkel = 0; % Kann gegen einen beliebigen Einbauwinkel ausgetauscht werden ACHTUNG!!! in [GRAD] !!!!!


    % Profilwiderstand Leitwerk
d_l_HLW = specs.HLW_d2l;
d_l_SLW = specs.HLW_d2l;

xu_l_HLW = 0.035;
xu_l_SLW = 0.035;

x_u_HLW = HLW.Fluegeltiefen_eta_oR * xu_l_HLW;
x_u_SLW = SLW.Fluegeltiefen_eta_oR * xu_l_SLW;

        % Trimmwiderstand
l_mue = Ergebnisse_Fluegel.l_mue;  
dx_SP = 1.5; %%%%%% Ein random wert angenommen!!!!!!!!!!!
c_M_0_F = FM.c_M_NP_F0;

    % Zusatzwidertsand
c_A_F_laufvar = 0:0.01:1; 

    % Inerferenzwiderstand
l_int_F = DT.l_i_I;
n_int_F = 4;

l_int_HLW = HLW.Fluegeltiefen_eta_oR(1,1);
n_int_HLW = 4;

l_int_SLW = SLW.Fluegeltiefen_eta_oR(1,1);
n_int_SLW = 2;

l_int_NC = Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.l_nacell_pylon;
n_int_NC = 4;

l_int_PYL = Ergebnisse_Fluegel.Fluegeltiefen_eta(1,0.3*10^3)*0.65;
n_int_PYL = 4;



%% ------------------------------------------------------------------------
%  -------------Fluegelwiderstand nach Diederich---------------------------
%  ------------------------------------------------------------------------


for n_iteration = 1:100 

    c_A_F = n_iteration * 10^(-2);

%% Profilwiderstand des Fluegels

% PS4 S.2, Formel 4 
k = 0.27 * specs.d_l + 100 * (specs.d_l)^4;

%PS4 S.3, Formel 7
Re_u = x_u * v_air / ISA.kin_visk(hoehe_CR);
Re = Ergebnisse_Fluegel.Fluegeltiefen_eta_oR * v_air / ISA.kin_visk(hoehe_CR);     % Annahme pleas confirm

% PS4 S.2, Formel 6
c_f_la_xu = 1.328./(sqrt(Re_u));
c_f_tu_xu = 0.455./(log(Re_u).^(2.58));
c_f_tu_l = 0.455./(log(Re).^(2.58));

% PS4 S.2, Formel 5
c_f = c_f_tu_l  - xu_l * (c_f_tu_xu - c_f_la_xu);

% PS4 S.2, Formel 3 
c_w_p_min_Re = 2 * c_f * (1 + k * cos(phi_50)^2);

% PS4 S.2, Formel 2 
c_w_p_eta = c_w_p_min_Re + 0.03 .*...
    (Ergebnisse_Auftriebsverteilung.c_a_eta(1,Ergebnisse_Fluegel.zaehlvariabele_eta_Ru:length(Ergebnisse_Auftriebsverteilung.c_a_eta))).^6;




% PS4 S.3, Formel 10 %%% Sieht alles noch sehr inkorekt aus

% test_integ = @(eta) c_w_p_eta .*  (Ergebnisse_Fluegel.Fluegeltiefen_eta_oR) ./ (Ergebnisse_Fluegel.l_m);
% c_w_p = integral(test_integ, 0, 1, 1001);
c_w_p = trapz(c_w_p_eta .*  (Ergebnisse_Fluegel.Fluegeltiefen_eta_oR) ./ (Ergebnisse_Fluegel.l_m)) * 10^(-3);
%test_sum = sum(c_w_p)
% % schnelltest plot
% plot(c_w_p)




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

% delta eps keine ahnung wie das berechnet werden soll
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% achtung nachfragen!!!!!!!!!!!!!!!!
eta_Ru = length(VWA.epsilon_eta_Ru)*10^(-3);
delta_eps = abs(VWA.epsilon - VWA.epsilon .* eta_Ru);

% PS4 S.4, Formel 11
c_w_ind(1,n_iteration) = c2 .* (c_A_F.^2)./(pi * Ergebnisse_Fluegel.streckung_phi25_max) +...
    c1 .* c_A_F .* delta_eps + c0 .* delta_eps^2;



% Transsonischer Widerstand

% PS4 S.5 Formel 18
k_vector = [0.758, 0.1, -0.090, 0, -0.100];

for n_DD = 1:5
    M_DD_profil_phi25_vec(1,n_DD) = k_vector(1,n_DD) .* c_A_F.^(n_DD-1);
end    
M_DD_profil_phi25(1,n_iteration) = sum(M_DD_profil_phi25_vec);

% PS4 S.5 Formel 17
delta_Ma = Ma_unendlich - M_DD_profil_phi25/(sqrt(cos(Ergebnisse_Fluegel.phi_25_max)));

% PS4 S.4 Formel 16
delta_c_WM = 0.002 * exp(60 * delta_Ma);



%% ------------------------------------------------------------------------

% Rumpfwiderstand

Re_Ru = (specs.l_rumpf * (specs.Ma_CR * ISA.a(hoehe_CR)))/(ISA.kin_visk(hoehe_CR));
% PS4 S.5 Formel 19
c_f_tu_Ru = (0.455)/(log(Re_Ru)^(2.58));

% PS4 S.5 Formel 21
k_Rumpf = 2.2 * (specs.D_rumpf/specs.l_rumpf)^(3/2) + 3.8 * (specs.D_rumpf/specs.l_rumpf)^(3);


% PS4 S.5 Formel 21
c_w_Ru_min = c_f_tu_Ru * (1 + k_Rumpf) * S_G_Ru/Ergebnisse_stat_Flaechenbelastung.F;

% Kommentar: Ich habe keine Berechnungen fuer den Widerstand unter
% betrachtung des Anstellwinkes des Rumpfes gemacht, das bedeutet, dise
% Rechnung gilt nur fuer CR-Zustand Wenn das hinzugefuegt werden muss PS4
% S.5 und fogend

% PS4 S.6 Formel 26
c_A_alpha_H = (pi * HLW.streckung_phi25)/(1 + sqrt(1 + 0.25 * HLW.streckung_phi25^2 * (tan(HLW.phi_50)^2 + (1 - specs.Ma_CR^2))));

% PS4 S.6 Formel 25 Abwindfaktor = delta_alpha_w/delta_alpha_oH 
Abwindfaktor = 1.75 * ((c_A_alpha_F)/(pi * Ergebnisse_Fluegel.streckung_phi25_max *...
    (Ergebnisse_Fluegel.lambda * (HLW.r/(Ergebnisse_Fluegel.b/2))^0.25 *...
    (1+ (abs(z_abstand))/((Ergebnisse_Fluegel.b/2))) )));

% PS4 S.5 Formel 24 annahme c_A_alpha_F = c_A_alpha_oH da der Rumpf kein
% wirklichen Auftrieb erzeugt
c_A_alpha = c_A_alpha_F * (1+ ((c_A_alpha_H)/(c_A_alpha_F)) *...
    (HLW.F/Ergebnisse_stat_Flaechenbelastung.F) * (1 - Abwindfaktor) );

% PS4 S.5 Formel 23
if zaehlvariabele_itt <= 1
    c_A_ges(1, n_iteration) = 0.5;

elseif  zaehlvariabele_itt > 1
    c_A_ges(1, n_iteration) = c_A_ges(1, n_iteration - 1);
end    
alpha_Rumpf_grad(1, n_iteration) = ((c_A_ges(1, n_iteration) - Ergebnisse_stat_Flaechenbelastung.C_A_CR)/(c_A_alpha)) * ...
    (180 / pi);
% PS4 S.5 Formel 22
c_w_R_zu_c_w_Rmin(1, n_iteration) = 0.000208 * abs(alpha_Rumpf_grad(1, n_iteration)).^3 + 0.00125 * abs(alpha_Rumpf_grad(1, n_iteration)).^2 + 0.029 * abs(alpha_Rumpf_grad(1, n_iteration)) + 1;

c_w_R(1, n_iteration) = c_w_R_zu_c_w_Rmin(1, n_iteration) .* c_w_Ru_min;



%% -------------------------------------------------------------------------
% Triebwerke analog zu Rumpf
k_TW = 0.2;
Re_TW = (l_TW * (specs.Ma_CR * ISA.a(hoehe_CR))) / (ISA.kin_visk(hoehe_CR));

% PS4 S.5 Formel 19
c_f_tu_TW = (0.455)/(log(Re_TW)^(2.58));

% PS4 S.5 Formel 23
alpha_TW_grad = alpha_Rumpf_grad + TW_Einbauwinkel;

% PS4 S.5 Formel 20
c_w_TW_min = c_f_tu_TW * (1 + k_TW) * S_G_TW/Ergebnisse_stat_Flaechenbelastung.F;

% PS4 S.5 Formel 22
c_w_TW_zu_c_w_TWmin = 0.000208 * abs(alpha_TW_grad).^3 + 0.00125 * abs(alpha_TW_grad).^2 + 0.029 * abs(alpha_TW_grad) + 1;

c_w_TW = c_w_TW_zu_c_w_TWmin * c_w_TW_min;

%% --------------------------------------------------------------------------
% Leitwerkswiderstand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Profilwiderstand Leitwerk


% PS4 S.3 Formel 7 angewendet auf Leitwerke
Re_u_HLW = x_u_HLW * v_air / ISA.kin_visk(hoehe_CR);
Re_u_SLW = x_u_SLW * v_air / ISA.kin_visk(hoehe_CR);

Re_HLW = HLW.Fluegeltiefen_eta_oR * v_air / ISA.kin_visk(hoehe_CR);
Re_SLW = SLW.Fluegeltiefen_eta_oR * v_air / ISA.kin_visk(hoehe_CR);


% PS4 S.2 Formel 6 angewendet auf Leitwerke
c_f_la_xu_HLW = (1.328)./(sqrt(Re_u_HLW));
c_f_la_xu_SLW = (1.328)./(sqrt(Re_u_SLW));

c_f_tu_xu_HLW = 0.455./(log(Re_u_HLW).^(2.58));
c_f_tu_xu_SLW = 0.455./(log(Re_u_SLW).^(2.58));

c_f_tu_l_HLW = 0.455./(log(Re_HLW).^(2.58));
c_f_tu_l_SLW = 0.455./(log(Re_SLW).^(2.58));

% PS4 S.2 Formel 5 angewendet auf Leitwerke
c_f_HLW = c_f_tu_l_HLW - xu_l_HLW .* (c_f_tu_xu_HLW - c_f_la_xu_HLW);
c_f_SLW = c_f_tu_l_SLW - xu_l_SLW .* (c_f_tu_xu_SLW - c_f_la_xu_SLW);

% PS4 S.6 Formel 28 % Annahme 
k_HLW = 2.7 .* d_l_HLW + 100 .* d_l_HLW.^4;
k_SLW = 2.7 .* d_l_SLW + 100 .* d_l_SLW.^4;

% PS4 S.6 Formel 27
c_w_HLW_min = 2 .* c_f_HLW .* (1+ k_HLW .* cos(HLW.phi_50).^2) .* ((HLW.F)/(Ergebnisse_stat_Flaechenbelastung.F));
c_w_SLW_min = 2 .* c_f_SLW .* (1+ k_SLW .* cos(SLW.phi_50).^2) .* ((SLW.F)/(Ergebnisse_stat_Flaechenbelastung.F));


% ------------------Trimwiderstand Leitwerke--------------------


% PS4 S.7 Formel 33 
qH_q = 0.95;
c_A_H(1,n_iteration) = (c_M_0_F + c_A_F * ((dx_SP)/(l_mue))) / ...
    (qH_q * ((HLW.F)/(Ergebnisse_Fluegel.F)) * ((HLW.r - dx_SP)/(l_mue)));

% PS4 S.7 Formel 29
c_A_ges(1,n_iteration) = c_A_F + c_A_H(1,n_iteration) * qH_q * ((HLW.F)/(Ergebnisse_stat_Flaechenbelastung.F));
c_A_ges_vec(zaehlvariabele_itt,1) = c_A_ges(1,n_iteration);

zaehlvariabele_itt = zaehlvariabele_itt + 1;


% PS4 S.4 Formel 15
tau_H = 1 - HLW.streckung_phi25 * (0.002 + 0.0084 * (HLW.lambda - 0.2)^2);

% PS4 S.7 Formel 30
c_w_trim(1,n_iteration) = ((c_A_H(1,n_iteration).^2)/(pi * HLW.streckung_phi25)) * ...
    ((1 + (5*10^(-6)) * (rad2deg(HLW.phi_25))^3)/(tau_H)) *...
    ((HLW.F)/(Ergebnisse_Fluegel.F));


%%--------------------------------------------------------------------------
% Zusatzwiderstand


% PS4 S.8 Formel 38
% Achtung c_A_F_laufvar soll eien Laufvariabele sein 
d_alpha_oH = c_A_F ./ c_A_alpha_F;

% PS4 S.8 Formel 37
d_alpha_w = Abwindfaktor .* d_alpha_oH;

% PS4 S.8 Formel 36
delta_c_w_H(1,n_iteration) = c_A_H(1,n_iteration) .* sin(d_alpha_w) .* qH_q .* ((HLW.F)./(Ergebnisse_stat_Flaechenbelastung.F));


%%--------------------------------------------------------------------------
% Interferenzwiderstand
    
        % Fluegel
Re_F_wurzel = l_int_F * v_air / ISA.kin_visk(hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_F = ((0.1369)/(Re_F_wurzel)) * l_int_F^2 * n_int_F;
        

        % HLW
Re_HLW_wurzel = l_int_HLW * v_air / ISA.kin_visk(hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_HLW = ((0.1369)/(Re_HLW_wurzel)) * l_int_HLW^2 * n_int_HLW;


        % SLW
Re_SLW_wurzel = l_int_SLW * v_air / ISA.kin_visk(hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_SLW = ((0.1369)/(Re_SLW_wurzel)) * l_int_SLW^2 * n_int_SLW;


        % TW Nacell
Re_NC_wurzel = l_int_NC * v_air / ISA.kin_visk(hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_NC = ((0.1369)/(Re_NC_wurzel)) * l_int_NC^2 * n_int_NC;


        % Pyl Nacell
Re_PYL_wurzel = l_int_PYL * v_air / ISA.kin_visk(hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_PYL = ((0.1369)/(Re_PYL_wurzel)) * l_int_PYL^2 * n_int_PYL;

    % zusammenfassung

c_w_int_fs = (c_w_int_F + c_w_int_HLW + ...
    c_w_int_SLW +c_w_int_NC + c_w_int_PYL)/...
    Ergebnisse_stat_Flaechenbelastung.F;

n_iteration_vec(1,n_iteration) = n_iteration * 10^(-2);

end


% Plot design

figure(1)
hold on
grid on
xlim([0, 0.05])
ylim([0, 1])


% Plot Profilwiderstand
%plot(c_w_p,   '.r') %n_iteration_vec,

% plot Induzierter Widersatnd
plot(c_w_ind, n_iteration_vec, '-red')
% 
% 
% Plot Transsonischer Widersatnd
plot(delta_c_WM, n_iteration_vec, '-green')
 
% Plot Rumpfwiderstand
plot(c_w_R, n_iteration_vec, '-k')
 
% Plot Widerstand Triebwerk
plot(c_w_TW,n_iteration_vec, '-m');


% Plot Widerstand Seitenleitwerk
%p(1) =
%plot(c_w_SLW_min);
% 
% % Plot Widerstand Hoehenleitwerk 
% % p(2) = plot()
% 
% % Plot Trimwiderstand HLW
% %p(3) = 
 plot(c_w_trim, n_iteration_vec,'-blue');
% 
% % Plot Zusatzwidersatnd
% %p(4) = plot()





