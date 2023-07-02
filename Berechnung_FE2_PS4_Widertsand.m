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

Annahmen.Flughoehe_CR = specs.flight_level * 10^2 ;     % in ft
Annahmen.hoehe_CR = round(unitsratio('m','ft')*Annahmen.Flughoehe_CR);



%% Getroffene Annahmen um Rechnungen vor berechnung der richtigen Werte durchfuehren zu können
        % es fehlen Werte als PS2 / PS3
    
    % Profilwiderstand des Fluegels  
Annahmen.phi_50 = Ergebnisse_Fluegel.phi_50; % tan((NP.versatz_HK + 0.5*DT.l_i_I - 0.5*DT.l_a)/(Ergebnisse_Fluegel.b/2 - specs.R_rumpf));          % Muss noch erfragen wie Annahmen.phi_50 brterchnet werden soll
Annahmen.v_air = specs.Ma_CR * ISA.a(Annahmen.hoehe_CR);
Annahmen.xu_l = 0.035;    % Wert zwischen 0.02 und 0.05
Annahmen.x_u = Ergebnisse_Fluegel.Fluegeltiefen_eta_oR * Annahmen.xu_l;      % Annahme, dass l die Fluegeltiefen an der jeweiligen Fosition auf dem Fluegel sind
    % induzierter Widerstand des Flügels
%c_A_F = linspace(0,1); % Ergebnisse_stat_Flaechenbelastung.C_A_CR;       % aus PS4 Formel 11

    % Transsonischer Widerstand
Annahmen.Ma_unendlich = specs.Ma_CR;

    % Rumpfwiderstand
Annahmen.S_G_Ru = Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.Sg_12;
Annahmen.c_A_alpha_F = GRA.c_a_anstieg;      % Annahme bitte ueberpruefen !!!!!!!!!!!!!!!!!!!!

Annahmen.z_abstand = 3; % Abstand zwischen Profilsehnen angenommen vergleiche Torenbeek s480
Annahmen.zaehlvariabele_itt = 1;

    % Triebwerke
Annahmen.S_G_TW = Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Turbine_Area; % Umspuehlte TW oberflaeche
Annahmen.l_TW = Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Turbine_length;
Annahmen.d_TW = Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.Turbine_Diameter;
Annahmen.TW_Einbauwinkel = 0; % Kann gegen einen beliebigen Einbauwinkel ausgetauscht werden ACHTUNG!!! in [GRAD] !!!!!


    % Profilwiderstand Leitwerk
Annahmen.d_l_HLW = specs.HLW_d2l;
Annahmen.d_l_SLW = specs.HLW_d2l;

Annahmen.xu_l_HLW = 0.035;
Annahmen.xu_l_SLW = 0.035;

Annahmen.x_u_HLW = HLW.Fluegeltiefen_eta_oR * Annahmen.xu_l_HLW;
Annahmen.x_u_SLW = SLW.Fluegeltiefen_eta_oR * Annahmen.xu_l_SLW;

        % Trimmwiderstand
Annahmen.l_mue = Ergebnisse_Fluegel.l_mue;  
Annahmen.x_SP_MAC = 1.5; %%%%%% Ein random wert angenommen!!!!!!!!!!!
Annahmen.c_M_0_F = FM.c_M_NP_F0;%* 0.1;
Annahmen.qH_q = 0.85; % hat Kristof gesagt urspruenglich mit 0.95 angenommen

    % Zusatzwidertsand
Annahmen.c_A_F_laufvar = 0:0.01:1; 

    % Inerferenzwiderstand
Annahmen.l_int_F = DT.l_i_I;
Annahmen.n_int_F = 4;

Annahmen.l_int_HLW = HLW.Fluegeltiefen_eta_oR(1,1);
Annahmen.n_int_HLW = 4;

Annahmen.l_int_SLW = SLW.Fluegeltiefen_eta_oR(1,1);
Annahmen.n_int_SLW = 2;

Annahmen.l_int_NC = Anteile_einzel_Massen_FE2.Airplane_Structure.Nacelle_group.l_nacell_pylon;
Annahmen.n_int_NC = 4;

Annahmen.l_int_PYL = Ergebnisse_Fluegel.Fluegeltiefen_eta(1,0.3*10^3)*0.65;
Annahmen.n_int_PYL = 4;

% Off design
stuetzstellen = 30;
Annahmen.Ma_off_D = linspace(0.1, 1, stuetzstellen);
Annahmen.kappa = 1.4;



% ------------------------------------------------------------------------
%% Fluegelwiderstand nach Diederich---------------------------
%  ------------------------------------------------------------------------


%% Funktionen direktzugriff


%PS4 S.3, Formel 7
FUN.Re_CR_fun = @(l_Re,v_Re) (l_Re .* v_Re) ./ ISA.kin_visk(Annahmen.hoehe_CR);

% PS4 S.2, Formel 6
FUN.c_f_la_fun = @(Re) 1.328./(sqrt(Re));
FUN.c_f_tu_fun = @(Re) 0.455./(log(Re).^(2.58));

% PS4 S.4, Formel 15
FUN.tau_fun = @(Streckung, lambda) 1 - Streckung * (0.002 + 0.0084 * (lambda - 0.2)^2);


save Getroffene_Annahmen_und_FUN.mat Annahmen FUN


c_A_F_testvec = linspace(0, 1, stuetzstellen);

for n_iteration = 1:stuetzstellen

c_A_F = c_A_F_testvec(1,n_iteration);   % n_iteration * 10^(-2);
    
c_A_F_off_D = ((2)/(Annahmen.kappa * ISA.p(Annahmen.hoehe_CR) * Annahmen.Ma_off_D(1,n_iteration))) * Ergebnisse_stat_Flaechenbelastung.Fleachenbelastung; 
v_air_off_D = Annahmen.Ma_off_D(1,n_iteration) * ISA.a(Annahmen.hoehe_CR);



c_A_F_off_D_vec(1,n_iteration) = c_A_F_off_D;

%% Profilwiderstand des Fluegels

% % PS4 S.2, Formel 4 
% k = 0.27 * specs.d_l + 100 * (specs.d_l)^4;
% 
% %PS4 S.3, Formel 7
% Re_u = Annahmen.x_u * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% Re = Ergebnisse_Fluegel.Fluegeltiefen_eta_oR * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);     % Annahme pleas confirm
% 
% % PS4 S.2, Formel 6
% c_f_la_xu = 1.328./(sqrt(Re_u));
% c_f_tu_xu = 0.455./(log(Re_u).^(2.58));
% c_f_tu_l = 0.455./(log(Re).^(2.58));
% 
% % PS4 S.2, Formel 5
% c_f = c_f_tu_l  - Annahmen.xu_l * (c_f_tu_xu - c_f_la_xu);
% 
% % PS4 S.2, Formel 3 
% c_w_p_min_Re = 2 * c_f * (1 + k * cos(Annahmen.phi_50)^2);
% 
% % PS4 S.2, Formel 2 
% c_w_p_eta = c_w_p_min_Re + 0.03 .*...
%     (Ergebnisse_Auftriebsverteilung.c_a_eta(1,Ergebnisse_Fluegel.zaehlvariabele_eta_Ru:length(Ergebnisse_Auftriebsverteilung.c_a_eta))).^6;
% 
% 
% 
% 
% % PS4 S.3, Formel 10 %%% Sieht alles noch sehr inkorekt aus
% 
% % test_integ = @(eta) c_w_p_eta .*  (Ergebnisse_Fluegel.Fluegeltiefen_eta_oR) ./ (Ergebnisse_Fluegel.l_m);
% % c_w_p = integral(test_integ, 0, 1, 1001);
% c_w_p(1,n_iteration) = trapz(c_w_p_eta .*  (Ergebnisse_Fluegel.Fluegeltiefen_eta_oR) ./ (Ergebnisse_Fluegel.l_m)) * 10^(-3); % Achtung Potenz kann inkoreckt sein
% %test_sum = sum(c_w_p)
% % % schnelltest plot
% % plot(c_w_p)


[c_w_p(1,n_iteration), c_w_p_min_Re(n_iteration,:)] = Profilwiderstand(Annahmen.v_air);

[c_w_p_off_D(1,n_iteration), c_w_p_min_Re_off_D(n_iteration,:)] = Profilwiderstand(v_air_off_D);

%% Induzierter Widerstand des Fluegels

% % PS4 S.4, Formel 12
%  c0 = VWA.c_AF_anstieg^2 * (0.0088 * Ergebnisse_Fluegel.lambda - 0.0051 * Ergebnisse_Fluegel.lambda^2) * (1 - 0.0006 * Ergebnisse_Fluegel.streckung_phi25_max^2);
% 
% % PS4 S.4, Formel 13
% % c1 = GRA.c_a_anstieg * (0.0134 * (Ergebnisse_Fluegel.lambda - 0.3) - 0.0037 * Ergebnisse_Fluegel.lambda^2);
% c1 = VWA.c_AF_anstieg * (0.0134 * (Ergebnisse_Fluegel.lambda - 0.3) - 0.0037 * Ergebnisse_Fluegel.lambda^2); %% nicht sicer mit welchem Auftriebsanstieg gerechnet werden muss
% % PS4 S.4, Formel 15
% tau = 1 - Ergebnisse_Fluegel.streckung_phi25_max * (0.002 + 0.0084 * (Ergebnisse_Fluegel.lambda-0.2)^2);
% 
% % PS4 S.4, Formel 14
% c2 = (1/tau) * (1 + (5*10^(-6)) * (rad2deg(Ergebnisse_Fluegel.phi_25_max))^3); 
% 
% % delta eps keine ahnung wie das berechnet werden soll
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% achtung nachfragen!!!!!!!!!!!!!!!!
% eta_Ru = length(VWA.epsilon_eta_Ru)*10^(-3);
% delta_eps = abs(VWA.epsilon - VWA.epsilon .* eta_Ru);
% 
% % PS4 S.4, Formel 11
% c_w_ind(1,n_iteration) = c2 .* (c_A_F.^2)./(pi * Ergebnisse_Fluegel.streckung_phi25_max) +...
%     c1 .* c_A_F .* delta_eps + c0 .* delta_eps^2;

c_w_ind(1,n_iteration) = Induzierter_W(c_A_F);
c_w_ind_off_D(1,n_iteration) = Induzierter_W(c_A_F_off_D);

%% Transsonischer Widerstand

% % PS4 S.5 Formel 18
% k_vector = [0.758, 0.1, -0.090, 0, -0.100];
% 
% for n_DD = 1:5
%     M_DD_profil_phi25_vec(1,n_DD) = k_vector(1,n_DD) .* c_A_F.^(n_DD-1);
% end    
% M_DD_profil_phi25(1,n_iteration) = sum(M_DD_profil_phi25_vec);
% 
% % PS4 S.5 Formel 17
% delta_Ma = Annahmen.Ma_unendlich - M_DD_profil_phi25/(sqrt(cos(Ergebnisse_Fluegel.phi_25_max)));
% 
% % PS4 S.4 Formel 16
% delta_c_WM = 0.002 * exp(60 * delta_Ma);

[delta_c_WM(1,n_iteration), delta_Ma(1,n_iteration)] = Transsonischer_W(Annahmen.Ma_unendlich,c_A_F);

[delta_c_WM_off_D(1,n_iteration), delta_Ma_off_D(1,n_iteration)] = Transsonischer_W(Annahmen.Ma_off_D(1,n_iteration), c_A_F_off_D);


% ------------------------------------------------------------------------
%% Rumpfwiderstand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Annahmen.zaehlvariabele_itt <= 1
    c_A_ges(1, n_iteration) = 0.522;
    c_A_ges_off_D(1, n_iteration) = 0.522;

elseif  Annahmen.zaehlvariabele_itt > 1
    c_A_ges(1, n_iteration) = c_A_ges(1, n_iteration - 1);
    c_A_ges(1,1) = c_A_ges(1, n_iteration);
    c_A_ges_off_D(1, n_iteration) = c_A_ges_off_D(1, n_iteration - 1);
    c_A_ges_off_D(1,1) = c_A_ges_off_D(1, n_iteration);
end   


% PS4 S.6 Formel 25 Abwindfaktor = delta_alpha_w/delta_alpha_oH 
Abwindfaktor = 1.75 * ((Annahmen.c_A_alpha_F)/(pi * Ergebnisse_Fluegel.streckung_phi25_max *...
    (Ergebnisse_Fluegel.lambda * (HLW.r/(Ergebnisse_Fluegel.b/2))^0.25 *...
    (1+ (abs(Annahmen.z_abstand))/((Ergebnisse_Fluegel.b/2))) )));

% Re_Ru = (specs.l_rumpf * (specs.Ma_CR * ISA.a(Annahmen.hoehe_CR)))/(ISA.kin_visk(Annahmen.hoehe_CR));
% % PS4 S.5 Formel 19
% c_f_tu_Ru = (0.455)/(log(Re_Ru)^(2.58));
% 
% % PS4 S.5 Formel 21
% k_Rumpf = 2.2 * (specs.D_rumpf/specs.l_rumpf)^(3/2) + 3.8 * (specs.D_rumpf/specs.l_rumpf)^(3);
% 
% 
% % PS4 S.5 Formel 21
% c_w_Ru_min = c_f_tu_Ru * (1 + k_Rumpf) * Annahmen.S_G_Ru/Ergebnisse_Fluegel.F;
% 
% % Kommentar: Ich habe keine Berechnungen fuer den Widerstand unter
% % betrachtung des Anstellwinkes des Rumpfes gemacht, das bedeutet, dise
% % Rechnung gilt nur fuer CR-Zustand Wenn das hinzugefuegt werden muss PFluegeltiefen_eta_oRS4
% % S.5 und fogend
% 
% % PS4 S.6 Formel 26
% c_A_alpha_H = (pi * HLW.streckung_phi25)/(1 + sqrt(1 + 0.25 * HLW.streckung_phi25^2 * (tan(HLW.phi_50)^2 + (1 - specs.Ma_CR^2))));
% 
% 
% 
% % PS4 S.5 Formel 24 annahme Annahmen.c_A_alpha_F = c_A_alpha_oH da der Rumpf kein
% % wirklichen Auftrieb erzeugt
% c_A_alpha = Annahmen.c_A_alpha_F * (1+ ((c_A_alpha_H)/(Annahmen.c_A_alpha_F)) *...
%     (HLW.F/Ergebnisse_Fluegel.F) * (1 - Abwindfaktor) );
% 
% % PS4 S.5 Formel 23
% 
% alpha_Rumpf_grad(1, n_iteration) = ((c_A_ges(1, n_iteration) - Ergebnisse_stat_Flaechenbelastung.C_A_CR)/(c_A_alpha)) * ...
%     (180 / pi);
% % PS4 S.5 Formel 22
% c_w_R_zu_c_w_Rmin(1, n_iteration) = 0.000208 * abs(alpha_Rumpf_grad(1, n_iteration)).^3 + 0.00125 * abs(alpha_Rumpf_grad(1, n_iteration)).^2 + 0.029 * abs(alpha_Rumpf_grad(1, n_iteration)) + 1;
% 
% c_w_R(1, n_iteration) = c_w_R_zu_c_w_Rmin(1, n_iteration) .* c_w_Ru_min;

[c_w_R(1, n_iteration), alpha_Rumpf_grad(1, n_iteration)] =...
    Rumpfwiderstand(specs.Ma_CR,Abwindfaktor,c_A_ges(1,n_iteration));
[c_w_R_off_D(1, n_iteration), alpha_Rumpf_grad_off_D(1, n_iteration)] =...
    Rumpfwiderstand(Annahmen.Ma_off_D(1,n_iteration),Abwindfaktor,c_A_ges_off_D(1,n_iteration));

% -------------------------------------------------------------------------
%% Triebwerke analog zu Rumpf
% -------------------------------------------------------------------------


% k_TW = 0.2;
% Re_TW = (Annahmen.l_TW * (specs.Ma_CR * ISA.a(Annahmen.hoehe_CR))) / (ISA.kin_visk(Annahmen.hoehe_CR));
% 
% % PS4 S.5 Formel 19
% c_f_tu_TW = (0.455)/(log(Re_TW)^(2.58));
% 
% % PS4 S.5 Formel 23
% alpha_TW_grad = alpha_Rumpf_grad + Annahmen.TW_Einbauwinkel;
% 
% % PS4 S.5 Formel 20
% c_w_TW_min = c_f_tu_TW * (1 + k_TW) * (Annahmen.S_G_TW * 2)/Ergebnisse_Fluegel.F;
% 
% % PS4 S.5 Formel 22
% c_w_TW_zu_c_w_TWmin = 0.000208 * abs(alpha_TW_grad).^3 + 0.00125 * abs(alpha_TW_grad).^2 + 0.029 * abs(alpha_TW_grad) + 1;
% 
% c_w_TW = c_w_TW_zu_c_w_TWmin .* c_w_TW_min;


c_w_TW(1, n_iteration) = Triebwerkswiderstand(Annahmen.v_air, alpha_Rumpf_grad(1, n_iteration));
c_w_TW_off_D(1, n_iteration) = Triebwerkswiderstand(v_air_off_D, alpha_Rumpf_grad_off_D(1, n_iteration));



% --------------------------------------------------------------------------
%% Leitwerkswiderstand
%---------------------------------------------------------------------------

% % PS4 S.3 Formel 7 angewendet auf Leitwerke
% Re_u_HLW = Annahmen.x_u_HLW * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% Re_u_SLW = Annahmen.x_u_SLW * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% 
% Re_HLW = HLW.Fluegeltiefen_eta_oR * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% Re_SLW = SLW.Fluegeltiefen_eta_oR * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% 
% 
% % PS4 S.2 Formel 6 angewendet auf Leitwerke
% c_f_la_xu_HLW = (1.328)./(sqrt(Re_u_HLW));
% c_f_la_xu_SLW = (1.328)./(sqrt(Re_u_SLW));
% 
% c_f_tu_xu_HLW = 0.455./(log(Re_u_HLW).^(2.58));
% c_f_tu_xu_SLW = 0.455./(log(Re_u_SLW).^(2.58));
% 
% c_f_tu_l_HLW = 0.455./(log(Re_HLW).^(2.58));
% c_f_tu_l_SLW = 0.455./(log(Re_SLW).^(2.58));
% 
% % PS4 S.2 Formel 5 angewendet auf Leitwerke
% c_f_HLW = c_f_tu_l_HLW - Annahmen.xu_l_HLW .* (c_f_tu_xu_HLW - c_f_la_xu_HLW);
% c_f_SLW = c_f_tu_l_SLW - Annahmen.xu_l_SLW .* (c_f_tu_xu_SLW - c_f_la_xu_SLW);
% 
% % PS4 S.6 Formel 28 % Annahme 
% k_HLW = 2.7 .* Annahmen.d_l_HLW + 100 .* Annahmen.d_l_HLW.^4;
% k_SLW = 2.7 .* Annahmen.d_l_SLW + 100 .* Annahmen.d_l_SLW.^4;
% 
% % PS4 S.6 Formel 27
% c_w_HLW_min = 2 .* c_f_HLW .* (1+ k_HLW .* cos(HLW.phi_50).^2) .* ((HLW.F)/(Ergebnisse_Fluegel.F));
% c_w_SLW_min = 2 .* c_f_SLW .* (1+ k_SLW .* cos(SLW.phi_50).^2) .* ((SLW.F)/(Ergebnisse_Fluegel.F));

[c_w_HLW_min, c_w_SLW_min] = Leitwerke_W(Annahmen.v_air);
[c_w_HLW_min_off_D, c_w_SLW_min_off_D] = Leitwerke_W(v_air_off_D);

c_W_HLW(1,n_iteration) = trapz(c_w_HLW_min)*10^(-2);
c_W_SLW(1,n_iteration) = trapz(c_w_SLW_min)*10^(-2);

c_W_HLW_off_D(1,n_iteration) = trapz(c_w_HLW_min_off_D)*10^(-2);
c_W_SLW_off_D(1,n_iteration) = trapz(c_w_SLW_min_off_D)*10^(-2);


%--------------------------------------------------------------------------
%% Trimwiderstand Leitwerke
% -------------------------------------------------------------------------

% PS4 S.7 Formel 33 

% % c_A_H(1,n_iteration) = (Annahmen.c_M_0_F + c_A_F * ((Annahmen.x_SP_MAC)/(Annahmen.l_mue))) / ...
% %     (Annahmen.qH_q * ((HLW.F)/(Ergebnisse_Fluegel.F)) * ((HLW.r - Annahmen.x_SP_MAC)/(Annahmen.l_mue)));
% 
% c_A_H(1,n_iteration) = (Annahmen.c_M_0_F + c_A_F * (0.1)) ./ ...Annahmen.zaehlvariabele_itt
%                         (Annahmen.qH_q * ((HLW.F)/(Ergebnisse_Fluegel.F)) *...
%                         (((HLW.r/Annahmen.l_mue) - (0.1))));
% 
% 
% % PS4 S.7 Formel 29
% c_A_ges(1,n_iteration) = c_A_F + c_A_H(1,n_iteration) * Annahmen.qH_q * ((HLW.F)/(Ergebnisse_Fluegel.F));
% 
% Annahmen.zaehlvariabele_itt = Annahmen.zaehlvariabele_itt + 1;
% 
% 
% % PS4 S.4 Formel 15
% tau_H = 1 - HLW.streckung_phi25 * (0.002 + 0.0084 * (HLW.lambda - 0.2)^2);
% 
% % PS4 S.7 Formel 30
% c_w_trim(1,n_iteration) = ((c_A_H(1,n_iteration).^2) ./ (pi * HLW.streckung_phi25)) .* ...
%     ((1 + (5*10^(-6)) .* (abs(rad2deg(HLW.phi_25))).^3)./(tau_H)) .*...
%     ((HLW.F)./(Ergebnisse_Fluegel.F));

[c_A_H(1,n_iteration), c_A_ges(1,n_iteration), c_w_trim(1,n_iteration)] = Leitwerke_Trim_W(c_A_F);
[c_A_H_off_D(1,n_iteration), c_A_ges_off_D(1,n_iteration), c_w_trim_off_D(1,n_iteration)] = Leitwerke_Trim_W(c_A_F_off_D);

%--------------------------------------------------------------------------
%% Zusatzwiderstand
%--------------------------------------------------------------------------

% PS4 S.8 Formel 38
% Achtung Annahmen.c_A_F_laufvar soll eien Laufvariabele sein 
% d_alpha_oH(1,n_iteration) = c_A_F ./ Annahmen.c_A_alpha_F;
% 
% % PS4 S.8 Formel 37
% d_alpha_w(1,n_iteration) = Abwindfaktor .* d_alpha_oH(1,n_iteration);
% 
% % PS4 S.8 Formel 36
% delta_c_w_H(1,n_iteration) = c_A_H(1,n_iteration) .* sin(d_alpha_w(1,n_iteration)) .*...
%     Annahmen.qH_q .* ((HLW.F)./(Ergebnisse_Fluegel.F));

[delta_c_w_H(1,n_iteration)] = Zusatz_W(c_A_F, Abwindfaktor, c_A_H(1, n_iteration));
[delta_c_w_H_off_D(1,n_iteration)] = Zusatz_W(c_A_F_off_D, Abwindfaktor, c_A_H_off_D(1, n_iteration));

%--------------------------------------------------------------------------
%% Interferenzwiderstand
%--------------------------------------------------------------------------    
        % Fluegel
% Re_F_wurzel = Annahmen.l_int_F * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% % PS4 S.9 Formel 40
% c_w_int_F = ((0.1369)/(Re_F_wurzel)^0.4) * Annahmen.l_int_F^2 * Annahmen.n_int_F;
%         
% 
%         % HLW
% Re_HLW_wurzel = Annahmen.l_int_HLW * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% % PS4 S.9 Formel 40
% c_w_int_HLW = ((0.1369)/(Re_HLW_wurzel)^0.4) * Annahmen.l_int_HLW^2 * Annahmen.n_int_HLW;
% 
% 
%         % SLW
% Re_SLW_wurzel = Annahmen.l_int_SLW * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% % PS4 S.9 Formel 40
% c_w_int_SLW = ((0.1369)/(Re_SLW_wurzel)) * Annahmen.l_int_SLW^2 * Annahmen.n_int_SLW;
% 
% 
%         % TW Nacell
% Re_NC_wurzel = Annahmen.l_int_NC * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% % PS4 S.9 Formel 40
% c_w_int_NC = ((0.1369)/(Re_NC_wurzel)^0.4) * Annahmen.l_int_NC^2 * Annahmen.n_int_NC;
% 
% 
%         % Pyl Nacell
% Re_PYL_wurzel = Annahmen.l_int_PYL * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% % PS4 S.9 Formel 40
% c_w_int_PYL = ((0.1369)/(Re_PYL_wurzel)^0.4) * Annahmen.l_int_PYL^2 * Annahmen.n_int_PYL;
% 
%     % zusammenfassung
% 
% c_w_int_fs(1,n_iteration) = (c_w_int_F + c_w_int_HLW + ...
%     c_w_int_SLW +c_w_int_NC + c_w_int_PYL)/...
%     Ergebnisse_Fluegel.F;


c_w_int_fs(1,n_iteration) = Interferenz_W(Annahmen.v_air);
c_w_int_fs_off_D(1,n_iteration) = Interferenz_W(v_air_off_D);





end



%% Ergebnisse Speichern
n_iteration_vec = linspace(0, 1, stuetzstellen); %n_iteration * 10^(-2);

% Vector mit zu plottenden Werten

x_vector = [c_W_SLW; c_W_HLW; c_w_int_fs; c_w_R; c_w_TW; c_w_trim; delta_c_w_H; c_w_p; c_w_ind; delta_c_WM];

for n_vec = 1:10
    if n_vec == 1
        x_vector_sum(n_vec,:) = x_vector(n_vec,:);
    else
        x_vector_sum(n_vec,:) = x_vector_sum(n_vec-1,:) + x_vector(n_vec,:);
    end
end



% Offdesign Vector
off_D_vector = [c_W_SLW_off_D; c_W_HLW_off_D; c_w_int_fs_off_D; c_w_R_off_D; c_w_TW_off_D; c_w_trim_off_D; delta_c_w_H_off_D; c_w_p_off_D; c_w_ind_off_D; delta_c_WM_off_D];

for n_vec_off_D = 1:10
    if n_vec_off_D == 1
        x_vector_sum_off_D(n_vec_off_D,:) = off_D_vector(n_vec_off_D,:);
    else
        x_vector_sum_off_D(n_vec_off_D,:) = x_vector_sum_off_D(n_vec_off_D - 1,:) + off_D_vector(n_vec_off_D,:);
    end
end


Ergebnisse_Widerstand_FE2.c_A_F_off_D_vec = c_A_F_off_D_vec;
Ergebnisse_Widerstand_FE2.c_A_ges = c_A_ges;
Ergebnisse_Widerstand_FE2.c_A_ges_off_D = c_A_ges_off_D;
Ergebnisse_Widerstand_FE2.c_A_H = c_A_H;
Ergebnisse_Widerstand_FE2.c_A_H_off_D = c_A_H_off_D;
Ergebnisse_Widerstand_FE2.c_W_HLW = c_W_HLW;
Ergebnisse_Widerstand_FE2.c_W_HLW_off_D = c_W_HLW_off_D;
Ergebnisse_Widerstand_FE2.c_W_SLW = c_W_SLW;
Ergebnisse_Widerstand_FE2.c_W_SLW_off_D = c_W_SLW_off_D;
Ergebnisse_Widerstand_FE2.c_w_ind  = c_w_ind;
Ergebnisse_Widerstand_FE2.c_w_ind_off_D = c_w_ind_off_D;
Ergebnisse_Widerstand_FE2.c_w_int_fs =  c_w_int_fs;
Ergebnisse_Widerstand_FE2.c_w_int_fs_off_D =  c_w_int_fs_off_D;
Ergebnisse_Widerstand_FE2.c_w_p =  c_w_p;
Ergebnisse_Widerstand_FE2.c_w_p_off_D =  c_w_p_off_D;
Ergebnisse_Widerstand_FE2.c_w_R =  c_w_R;
Ergebnisse_Widerstand_FE2.c_w_R_off_D =  c_w_R_off_D;
Ergebnisse_Widerstand_FE2.c_w_trim =  c_w_trim;
Ergebnisse_Widerstand_FE2.c_w_trim_off_D = c_w_trim_off_D ;
Ergebnisse_Widerstand_FE2.c_w_TW =  c_w_TW;
Ergebnisse_Widerstand_FE2.c_w_TW_off_D =  c_w_TW_off_D;
Ergebnisse_Widerstand_FE2.delta_c_w_H =  delta_c_w_H;
Ergebnisse_Widerstand_FE2.delta_c_w_H_off_D =  delta_c_w_H_off_D;
Ergebnisse_Widerstand_FE2.delta_c_WM =  delta_c_WM;
Ergebnisse_Widerstand_FE2.delta_c_WM_off_D = delta_c_WM_off_D ;
Ergebnisse_Widerstand_FE2.n_iteration_vec =  n_iteration_vec;
Ergebnisse_Widerstand_FE2.alpha_Rumpf_grad =  alpha_Rumpf_grad;
Ergebnisse_Widerstand_FE2.alpha_Rumpf_grad_off_D =  alpha_Rumpf_grad_off_D;
Ergebnisse_Widerstand_FE2.Annahmen = Annahmen;
Ergebnisse_Widerstand_FE2.FUN = FUN;
Ergebnisse_Widerstand_FE2.x_vector = x_vector;
Ergebnisse_Widerstand_FE2.off_D_vector = off_D_vector;
Ergebnisse_Widerstand_FE2.x_vector_sum = x_vector_sum;
Ergebnisse_Widerstand_FE2.x_vector_sum_off_D = x_vector_sum_off_D;
Ergebnisse_Widerstand_FE2.c_W_ges = x_vector_sum(10,:);
Ergebnisse_Widerstand_FE2.c_W_ges_off_D = x_vector_sum_off_D(10,:);
Ergebnisse_Widerstand_FE2.c_W_ges_inkomp = x_vector_sum(9,:);
Ergebnisse_Widerstand_FE2.c_W_ges_off_D_inkomp = x_vector_sum_off_D(9,:);
Ergebnisse_Widerstand_FE2.delta_Ma = delta_Ma;
Ergebnisse_Widerstand_FE2.delta_Ma_off_D = delta_Ma_off_D;
Ergebnisse_Widerstand_FE2.c_w_p_min_Re = c_w_p_min_Re; % Matrix mit x=Anz Stützstellen y=Anz Re zahlen über den umspühlten Flügel ohne Rumpf  
Ergebnisse_Widerstand_FE2.c_w_p_min_Re_off_D = c_w_p_min_Re_off_D; % Matrix mit x=Anz Stützstellen y=Anz Re zahlen über den umspühlten Flügel ohne Rumpf  

save Ergebnisse_Widerstand_FE2.mat Ergebnisse_Widerstand_FE2;


%% Plot design

% figure(1)
% hold on
% grid on
% xlim([0, 0.05])
% ylim([0, 1])
% 
% 
% 
% % plot Leitwerke SLW und HLW
% 
% p(1) = plot(c_W_SLW, n_iteration_vec, '.-b');   % SLW
% 
% p(2) = plot(c_W_HLW, n_iteration_vec, '.-c');   % HLW
% 
% % Interferenzwiderstand
% 
% p(3) =  plot(c_w_int_fs, n_iteration_vec, '.green');
% 
% % Plot Rumpfwiderstand
% p(4) = plot(c_w_R, n_iteration_vec, '-k');
% 
% 
% % Plot Widerstand Triebwerk
% p(5) = plot(c_w_TW,n_iteration_vec, '-m');
% 
% 
% % Plot Trimwiderstand HLW
% 
% p(6) = plot(c_w_trim, n_iteration_vec,'-blue');
% 
% % Abwindwiderstand
% p(7) = plot(delta_c_w_H, n_iteration_vec, 'c');
% 
% % Plot Profilwiderstand 
% 
% p(8) = plot(c_w_p, n_iteration_vec, '-.k');
% 
% % plot Induzierter Widerstand
% p(9) = plot(c_w_ind, n_iteration_vec, '-red');
% 
% % Plot Transsonischer Widersatnd
% p(10) = plot(delta_c_WM, n_iteration_vec, '-green');
%  
% 
% % Testplot Widerstände aufaddiert
% widerstaende_aufaddiert = c_w_ind + delta_c_WM + c_w_R + c_w_TW + c_w_trim + c_w_int_fs + delta_c_w_H + c_W_HLW + c_W_SLW;
% plot(widerstaende_aufaddiert, n_iteration_vec, '*r')
% 
% 
% legend(p([1:10]),{'+ SLW', '+ HLW', '+ Interferenz', '+ Rumpf', '+ Triebwerk', '+ Trimmung', '+ Abwind', '+ Profil', '+ ind. Widerstand', '+ Wellenwiderstand'},'Location','southeast','FontSize',18);
% title('Kumulative Widerstandspolare')
% xlabel('c_{W}');
% ylabel('c_A');

% 
% % atomatisches Plotten
% figure(2)
% 
% hold on
% grid on
% xlim([0, 0.05])
% ylim([0, 1])
% 
% % Anzahl der Plots festlegen
% numPlots = 10 %6;     % muss veraendert werden um off Design noch zu plotten
% 
% % Vector mit zu plottenden Werten
% 
% 
% %x_vector = [c_W_SLW; c_W_HLW; c_w_int_fs; c_w_R; c_w_TW; c_w_trim; delta_c_w_H; c_w_p; c_w_ind; delta_c_WM; off_D];
% 
% % Farbverlauf definieren
% colorStart = [0, 1, 0];   % Startfarbe (RGB)
% colorEnd = [0, 0, 0];     % Endfarbe (RGB)
% 
% % Farbwerte für jeden Plot berechnen
% %colors = zeros(numPlots/2, 3);
% for n_color = 1:ceil(numPlots/2)
%     colors(n_color, :) = colorStart + (n_color-1) * (colorEnd - colorStart) / ((numPlots/2)-1);
% end
% 
% colors = vertcat(colors, colors);
% 
% 
% 
% % Linienarten definieren
% lineStyles = {'--', ':', '-.', '-', '--', ':', '-.', '-', '--', ':', '-.', '-' }; % Gestrichelt, Gepunktet
% 
% % Vector mit zu plottenden Werten
% 
% x_vector = [c_W_SLW; c_W_HLW; c_w_int_fs; c_w_R; c_w_TW; c_w_trim; delta_c_w_H; c_w_p; c_w_ind; delta_c_WM];
% 
% for n_vec = 1:numPlots
%     if n_vec == 1
%         x_vector_sum(n_vec,:) = x_vector(n_vec,:);
%     else
%         x_vector_sum(n_vec,:) = x_vector_sum(n_vec-1,:) + x_vector(n_vec,:)
%     end
% end
% 
% % Offdesign Vector
% off_D_vector = [c_W_SLW_off_D; c_W_HLW_off_D; c_w_int_fs_off_D; c_w_R_off_D; c_w_TW_off_D; c_w_trim_off_D; delta_c_w_H_off_D; c_w_p_off_D; c_w_ind_off_D; delta_c_WM_off_D];
% 
% for n_vec_off_D = 1:length(off_D_vector)
%     if n_vec_off_D == 1
%         x_vector_sum_off_D(n_vec_off_D,:) = off_D_vector(n_vec_off_D,:);
%     else
%         x_vector_sum_off_D(n_vec_off_D,:) = x_vector_sum_off_D(n_vec_off_D-1,:) + off_D_vector(n_vec_off_D,:);
%     end
% end
% 
% 
% % Plots erstellen
% 
% atoplots = cell(numPlots, 1);
% 
% for n_plot = 1:numPlots
%    autoplot(n_plot,1) = plot((x_vector_sum(n_plot,:)), n_iteration_vec, 'Color', colors(n_plot, :), 'LineStyle', lineStyles{1,n_plot});
% end
% plot(off_D_vector,c_A_F_off_D,'red')
% 
% 
% % schoen machen des Plots 
% 
% % legend(autoplot([1:10]),{'+ SLW', '+ HLW', '+ Interferenz', '+ Rumpf', '+ Triebwerk', '+ Trimmung', '+ Abwind', '+ Profil', '+ ind. Widerstand', '+ Wellenwiderstand', '+Off-Design'},...
% %     'Location','southeast','FontSize',18);
% title('Kumulative Widerstandspolare')
% xlabel('c_W');
% ylabel('c_A');
% 
% hold off;


%--------------------------------------------------------------------------
%% Funktionen zur Berechnung der einzelnen Widerstandsteile
%--------------------------------------------------------------------------

%% Funktion Profilwiderstand
function [c_w_p, c_w_p_min_Re] = Profilwiderstand(v_gegeben)

load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Auftrieb_Momente.mat;
load Getroffene_Annahmen_und_FUN.mat;

% PS4 S.2, Formel 4 
k = 0.27 * specs.d_l + 100 * (specs.d_l)^4;

%PS4 S.3, Formel 7
Re_u = FUN.Re_CR_fun(Annahmen.x_u, v_gegeben);                                           %Annahmen.x_u * v_gegeben / ISA.kin_visk(Annahmen.hoehe_CR);
Re_oR = FUN.Re_CR_fun(Ergebnisse_Fluegel.Fluegeltiefen_eta_oR, v_gegeben);       %Ergebnisse_Fluegel.Fluegeltiefen_eta_oR * v_gegeben / ISA.kin_visk(Annahmen.hoehe_CR);     % Annahme pleas confirm

% PS4 S.2, Formel 6
c_f_la_xu = FUN.c_f_la_fun(Re_u);                               %1.328./(sqrt(Re_u));
c_f_tu_xu = FUN.c_f_tu_fun(Re_u);                               %0.455./(log(Re_u).^(2.58));
c_f_tu_l = FUN.c_f_la_fun(Re_oR);                               % 0.455./(log(Re).^(2.58));

% PS4 S.2, Formel 5
c_f = c_f_tu_l  - Annahmen.xu_l * (c_f_tu_xu - c_f_la_xu);

% PS4 S.2, Formel 3 
c_w_p_min_Re = 2 * c_f * (1 + k * cos(Annahmen.phi_50)^2);

% PS4 S.2, Formel 2 
c_w_p_eta = c_w_p_min_Re + 0.03 .*...
    (Ergebnisse_Auftriebsverteilung.c_a_eta(1,Ergebnisse_Fluegel.zaehlvariabele_eta_Ru:length(Ergebnisse_Auftriebsverteilung.c_a_eta))).^6;

% PS4 S.3, Formel 10 %%% Sieht alles noch sehr inkorekt aus
c_w_p = trapz(c_w_p_eta .*  (Ergebnisse_Fluegel.Fluegeltiefen_eta_oR) ./ (Ergebnisse_Fluegel.l_m)) * 10^(-3); % Achtung Potenz kann inkoreckt sein

end

%% Funktion Induzierter Widerstand

function [c_w_ind] = Induzierter_W(c_A_F)

load Projekt_specs.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Auftrieb_Momente.mat;
load Getroffene_Annahmen_und_FUN.mat;


% PS4 S.4, Formel 12
c0 = VWA.c_AF_anstieg^2 * (0.0088 * Ergebnisse_Fluegel.lambda - 0.0051 * Ergebnisse_Fluegel.lambda^2) * (1 - 0.0006 * Ergebnisse_Fluegel.streckung_phi25_max^2);

% PS4 S.4, Formel 13
% c1 = GRA.c_a_anstieg * (0.0134 * (Ergebnisse_Fluegel.lambda - 0.3) - 0.0037 * Ergebnisse_Fluegel.lambda^2);
c1 = VWA.c_AF_anstieg * (0.0134 * (Ergebnisse_Fluegel.lambda - 0.3) - 0.0037 * Ergebnisse_Fluegel.lambda^2); %% nicht sicer mit welchem Auftriebsanstieg gerechnet werden muss
% PS4 S.4, Formel 15
%tau = 1 - Ergebnisse_Fluegel.streckung_phi25_max * (0.002 + 0.0084 * (Ergebnisse_Fluegel.lambda-0.2)^2);
tau = FUN.tau_fun(Ergebnisse_Fluegel.streckung_phi25_max, Ergebnisse_Fluegel.lambda);

% PS4 S.4, Formel 14
c2 = (1/tau) * (1 + (5*10^(-6)) * (rad2deg(Ergebnisse_Fluegel.phi_25_max))^3); 

% delta eps keine ahnung wie das berechnet werden soll
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% achtung nachfragen!!!!!!!!!!!!!!!!
eta_Ru = length(VWA.epsilon_eta_Ru)*10^(-3);
delta_eps = abs(VWA.epsilon - VWA.epsilon .* eta_Ru);

% PS4 S.4, Formel 11
c_w_ind = c2 .* (c_A_F.^2)./(pi * Ergebnisse_Fluegel.streckung_phi25_max) +...
    c1 .* c_A_F .* delta_eps + c0 .* delta_eps^2;

end

%% Funktion Wellenwiderstand

function [delta_c_WM, delta_Ma] = Transsonischer_W(Ma_unendlich, c_A_F)

load Ergebnisse_Fluegel_Tank_NP.mat;

% PS4 S.5 Formel 18
k_vector = [0.758, 0.1, -0.090, 0, -0.100];

for n_DD = 1:5
    M_DD_profil_phi25_vec(1,n_DD) = k_vector(1,n_DD) .* c_A_F.^(n_DD-1);
end    
M_DD_profil_phi25 = sum(M_DD_profil_phi25_vec);

% PS4 S.5 Formel 17
delta_Ma = Ma_unendlich - M_DD_profil_phi25./(sqrt(cos(Ergebnisse_Fluegel.phi_25_max)));

% PS4 S.4 Formel 16
delta_c_WM = 0.002 * exp(60 * delta_Ma);

end

%% Funktion Rumpfwiderstand

function [c_w_R, alpha_Rumpf_grad] = Rumpfwiderstand(Machzahl, Abwindfaktor, c_A_ges)


load Projekt_specs.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Auftrieb_Momente.mat;
load Ergebnisse_Leitwerke.mat;
load Getroffene_Annahmen_und_FUN.mat;



%Re_Ru = (specs.l_rumpf * (specs.Ma_CR * ISA.a(Annahmen.hoehe_CR)))/(ISA.kin_visk(Annahmen.hoehe_CR));
Re_Ru = FUN.Re_CR_fun(specs.l_rumpf, Annahmen.v_air);

% PS4 S.5 Formel 19
c_f_tu_Ru = FUN.c_f_tu_fun(Re_Ru); %(0.455)/(log(Re_Ru)^(2.58));

% PS4 S.5 Formel 21
k_Rumpf = 2.2 * (specs.D_rumpf / specs.l_rumpf)^(3/2) + 3.8 * (specs.D_rumpf/specs.l_rumpf)^(3);


% PS4 S.5 Formel 21
c_w_Ru_min = c_f_tu_Ru * (1 + k_Rumpf) * Annahmen.S_G_Ru/Ergebnisse_Fluegel.F;

% Kommentar: Ich habe keine Berechnungen fuer den Widerstand unter
% betrachtung des Anstellwinkes des Rumpfes gemacht, das bedeutet, dise
% Rechnung gilt nur fuer CR-Zustand Wenn das hinzugefuegt werden muss PFluegeltiefen_eta_oRS4
% S.5 und fogend

% PS4 S.6 Formel 26
c_A_alpha_H = (pi * HLW.streckung_phi25)/...
    (1 + sqrt(1 + 0.25 * HLW.streckung_phi25^2 * (tan(HLW.phi_50)^2 + (1 - Machzahl^2))));

% PS4 S.6 Formel 25 Abwindfaktor = delta_alpha_w/delta_alpha_oH 
% Abwindfaktor = 1.75 * ((Annahmen.c_A_alpha_F)/(pi * Ergebnisse_Fluegel.streckung_phi25_max *...
%     (Ergebnisse_Fluegel.lambda * (HLW.r/(Ergebnisse_Fluegel.b/2))^0.25 *...
%     (1+ (abs(Annahmen.z_abstand))/((Ergebnisse_Fluegel.b/2))) )));

% PS4 S.5 Formel 24 annahme Annahmen.c_A_alpha_F = c_A_alpha_oH da der Rumpf kein
% wirklichen Auftrieb erzeugt
c_A_alpha = Annahmen.c_A_alpha_F * (1+ ((c_A_alpha_H)/(Annahmen.c_A_alpha_F)) *...
    (HLW.F/Ergebnisse_Fluegel.F) * (1 - Abwindfaktor) );

% PS4 S.5 Formel 23
% if Annahmen.zaehlvariabele_itt <= 1
%     c_A_ges(1, n_iteration) = 0.522;
% 
% elseif  Annahmen.zaehlvariabele_itt > 1
%     c_A_ges(1, n_iteration) = c_A_ges(1, n_iteration - 1);
%     c_A_ges(1,1) = c_A_ges(1, n_iteration)
% end    
alpha_Rumpf_grad = ((c_A_ges - Ergebnisse_stat_Flaechenbelastung.C_A_CR)/(c_A_alpha)) * ...
    (180 / pi);
% PS4 S.5 Formel 22
c_w_R_zu_c_w_Rmin = 0.000208 * abs(alpha_Rumpf_grad).^3 + 0.00125 * abs(alpha_Rumpf_grad).^2 + 0.029 * abs(alpha_Rumpf_grad) + 1;

c_w_R = c_w_R_zu_c_w_Rmin .* c_w_Ru_min;

end

%% Funktion Triebwerkswiderstand

function [c_w_TW] = Triebwerkswiderstand(v_air, alpha_Rumpf_grad)

load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Getroffene_Annahmen_und_FUN.mat;

k_TW = 0.2;
Re_TW = FUN.Re_CR_fun(Annahmen.l_TW, v_air); %(Annahmen.l_TW * (specs.Ma_CR * ISA.a(Annahmen.hoehe_CR))) / (ISA.kin_visk(Annahmen.hoehe_CR));

% PS4 S.5 Formel 19
c_f_tu_TW = FUN.c_f_tu_fun(Re_TW); % (0.455)/(log(Re_TW)^(2.58));

% PS4 S.5 Formel 23
alpha_TW_grad = alpha_Rumpf_grad + Annahmen.TW_Einbauwinkel;

% PS4 S.5 Formel 20
c_w_TW_min = c_f_tu_TW * (1 + k_TW) * (Annahmen.S_G_TW * 2)/Ergebnisse_Fluegel.F;

% PS4 S.5 Formel 22
c_w_TW_zu_c_w_TWmin = 0.000208 * abs(alpha_TW_grad).^3 + 0.00125 * abs(alpha_TW_grad).^2 + 0.029 * abs(alpha_TW_grad) + 1;

c_w_TW = c_w_TW_zu_c_w_TWmin .* c_w_TW_min;

end

%% Funktion Leitwerkswiderstand

function  [c_w_HLW_min, c_w_SLW_min] = Leitwerke_W(v_air)

load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Leitwerke.mat;
load Getroffene_Annahmen_und_FUN.mat;

% Profilwiderstand Leitwerk


% PS4 S.3 Formel 7 angewendet auf Leitwerke
Re_u_HLW = FUN.Re_CR_fun(Annahmen.x_u_HLW, v_air);  %Annahmen.x_u_HLW * v_air / ISA.kin_visk(Annahmen.hoehe_CR);
Re_u_SLW = FUN.Re_CR_fun(Annahmen.x_u_SLW, v_air);  %Annahmen.x_u_SLW * v_air / ISA.kin_visk(Annahmen.hoehe_CR);

Re_HLW = FUN.Re_CR_fun(HLW.Fluegeltiefen_eta_oR, v_air);    % HLW.Fluegeltiefen_eta_oR .* v_air ./ ISA.kin_visk(Annahmen.hoehe_CR);
Re_SLW = FUN.Re_CR_fun(SLW.Fluegeltiefen_eta_oR, v_air);    % SLW.Fluegeltiefen_eta_oR .* v_air ./ ISA.kin_visk(Annahmen.hoehe_CR);


% PS4 S.2 Formel 6 angewendet auf Leitwerke
c_f_la_xu_HLW = FUN.c_f_la_fun(Re_u_HLW); % (1.328)./(sqrt(Re_u_HLW));
c_f_la_xu_SLW = FUN.c_f_la_fun(Re_u_SLW); % (1.328)./(sqrt(Re_u_SLW));

c_f_tu_xu_HLW = FUN.c_f_tu_fun(Re_u_HLW); % 0.455./(log(Re_u_HLW).^(2.58));
c_f_tu_xu_SLW = FUN.c_f_tu_fun(Re_u_SLW); % 0.455./(log(Re_u_SLW).^(2.58));

c_f_tu_l_HLW = FUN.c_f_tu_fun(Re_HLW); % 0.455./(log(Re_HLW).^(2.58));
c_f_tu_l_SLW = FUN.c_f_tu_fun(Re_SLW); % 0.455./(log(Re_SLW).^(2.58));

% PS4 S.2 Formel 5 angewendet auf Leitwerke
c_f_HLW = c_f_tu_l_HLW - Annahmen.xu_l_HLW .* (c_f_tu_xu_HLW - c_f_la_xu_HLW);
c_f_SLW = c_f_tu_l_SLW - Annahmen.xu_l_SLW .* (c_f_tu_xu_SLW - c_f_la_xu_SLW);

% PS4 S.6 Formel 28 % Annahme 
k_HLW = 2.7 .* Annahmen.d_l_HLW + 100 .* Annahmen.d_l_HLW.^4;
k_SLW = 2.7 .* Annahmen.d_l_SLW + 100 .* Annahmen.d_l_SLW.^4;

% PS4 S.6 Formel 27
c_w_HLW_min = 2 .* c_f_HLW .* (1+ k_HLW .* cos(HLW.phi_50).^2) .* ((HLW.F)/(Ergebnisse_Fluegel.F));
c_w_SLW_min = 2 .* c_f_SLW .* (1+ k_SLW .* cos(SLW.phi_50).^2) .* ((SLW.F)/(Ergebnisse_Fluegel.F));

end

%% Funktion Trimwiderstand Leitwerke

function [c_A_H, c_A_ges, c_w_trim] = Leitwerke_Trim_W(c_A_F)


load Projekt_specs.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Leitwerke.mat;
load Getroffene_Annahmen_und_FUN.mat;


% PS4 S.7 Formel 33 

% c_A_H(1,n_iteration) = (Annahmen.c_M_0_F + c_A_F * ((Annahmen.x_SP_MAC)/(Annahmen.l_mue))) / ...
%     (Annahmen.qH_q * ((HLW.F)/(Ergebnisse_Fluegel.F)) * ((HLW.r - Annahmen.x_SP_MAC)/(Annahmen.l_mue)));

c_A_H = (Annahmen.c_M_0_F + c_A_F * (0.1)) ./ ...
                        (Annahmen.qH_q * ((HLW.F)/(Ergebnisse_Fluegel.F)) *...
                        (((HLW.r/Annahmen.l_mue) - (0.1))));


% PS4 S.7 Formel 29
c_A_ges = c_A_F + c_A_H * Annahmen.qH_q * ((HLW.F)/(Ergebnisse_Fluegel.F));
% c_A_ges_vec(Annahmen.zaehlvariabele_itt,1) = c_A_ges(1,n_iteration);

Annahmen.zaehlvariabele_itt = Annahmen.zaehlvariabele_itt + 1;


% PS4 S.4 Formel 15
tau_H = FUN.tau_fun(HLW.streckung_phi25, HLW.lambda); %1 - HLW.streckung_phi25 * (0.002 + 0.0084 * (HLW.lambda - 0.2)^2);

% PS4 S.7 Formel 30
c_w_trim = ((c_A_H.^2) ./ (pi * HLW.streckung_phi25)) .* ...
    ((1 + (5*10^(-6)) .* (abs(rad2deg(HLW.phi_25))).^3)./(tau_H)) .*...
    ((HLW.F)./(Ergebnisse_Fluegel.F));


end

%% Funktion Zusatzwiderstand

function [delta_c_w_H] = Zusatz_W(c_A_F, Abwindfaktor, c_A_H)


load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Leitwerke.mat;
load Getroffene_Annahmen_und_FUN.mat;

% PS4 S.8 Formel 38
% Achtung Annahmen.c_A_F_laufvar soll eien Laufvariabele sein 
d_alpha_oH = c_A_F ./ Annahmen.c_A_alpha_F;

% PS4 S.8 Formel 37
d_alpha_w = Abwindfaktor .* d_alpha_oH;

% PS4 S.8 Formel 36
delta_c_w_H = c_A_H .* sin(d_alpha_w) .*...
    Annahmen.qH_q .* ((HLW.F)./(Ergebnisse_Fluegel.F));

end

%% Funktion Interferenzwiderstand

function [c_w_int_fs] = Interferenz_W(v_air)

load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Getroffene_Annahmen_und_FUN.mat;

Re_F_wurzel = FUN.Re_CR_fun(Annahmen.l_int_F, v_air); % Annahmen.l_int_F * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_F = ((0.1369)/(Re_F_wurzel)^0.4) * Annahmen.l_int_F^2 * Annahmen.n_int_F;
        

        % HLW
Re_HLW_wurzel = FUN.Re_CR_fun(Annahmen.l_int_HLW, v_air); %  Annahmen.l_int_HLW * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_HLW = ((0.1369)/(Re_HLW_wurzel)^0.4) * Annahmen.l_int_HLW^2 * Annahmen.n_int_HLW;


        % SLW
Re_SLW_wurzel = FUN.Re_CR_fun(Annahmen.l_int_SLW, v_air); %  Annahmen.l_int_SLW * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_SLW = ((0.1369)/(Re_SLW_wurzel)) * Annahmen.l_int_SLW^2 * Annahmen.n_int_SLW;


        % TW Nacell
Re_NC_wurzel = FUN.Re_CR_fun(Annahmen.l_int_NC, v_air); %  Annahmen.l_int_NC * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_NC = ((0.1369)/(Re_NC_wurzel)^0.4) * Annahmen.l_int_NC^2 * Annahmen.n_int_NC;


        % Pyl Nacell
Re_PYL_wurzel = FUN.Re_CR_fun(Annahmen.l_int_PYL, v_air); %  Annahmen.l_int_PYL * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_PYL = ((0.1369)/(Re_PYL_wurzel)^0.4) * Annahmen.l_int_PYL^2 * Annahmen.n_int_PYL;

    % zusammenfassung

c_w_int_fs = (c_w_int_F + c_w_int_HLW + ...
    c_w_int_SLW +c_w_int_NC + c_w_int_PYL)/...
    Ergebnisse_Fluegel.F;

end
