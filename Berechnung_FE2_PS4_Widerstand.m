function Berechnung_FE2_PS4_Widerstand

clc
clear all
close all


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





%% Laden von Werten
load Projekt_specs.mat;
load Ergebnisse_Massen_FE2.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Auftrieb_Momente.mat;
load Ergebnisse_Leitwerke.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Schwerpunkt.mat;

% Bereitstellen der Funktionen zur Widerstandsberechnung aus Unterordner 
addpath('Unterfunktionen Widerstand');


%% Allgemeine Variabelen

Annahmen.Flughoehe_CR = specs.flight_level * 10^2 ;     % in ft
Annahmen.hoehe_CR = round(unitsratio('m','ft')*Annahmen.Flughoehe_CR);

stuetzstellen = 500;
c_A_F = linspace(0, 1, stuetzstellen);
Ma_off_D = linspace(0, 2, stuetzstellen).';


%% Getroffene Annahmen um Rechnungen vor berechnung der richtigen Werte durchfuehren zu können
        % es fehlen Werte als PS2 / PS3
    
    % Profilwiderstand des Fluegels  

    
    %% Wichtig MUSS NOCH VERAENDERT WERDEN
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Annahmen.x_SP_MAC = 1.5; %%%%%% Ein random wert angenommen!!!!!!!!!!!
%     Annahmen.x_NP_MAC_oH = 0.5;
    Faktor = 0.5; %% Achtung wert angepasst, muss nochmal mit leon drüber schauen (stand 28.07. 0Uhr) 0.5
    Annahmen.dx_SP_lmue =  Delta_CG_MAC_durch_lmue * Faktor ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Annahmen.z_abstand = StatStab.z_abstand; % Abstand zwischen Profilsehnen angenommen vergleiche Torenbeek s480
 
    Annahmen.c_M_0_F = FM.c_M_NP_F0; %-0.1; % %* 0.1; % Wert aus FE1 ist zu klein 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    

    
    
    Annahmen.phi_50 = Ergebnisse_Fluegel.phi_50; % tan((NP.versatz_HK + 0.5*DT.l_i_I - 0.5*DT.l_a)/(Ergebnisse_Fluegel.b/2 - specs.R_rumpf));          % Muss noch erfragen wie Annahmen.phi_50 brterchnet werden soll

Annahmen.xu_l = 0.035;    % Wert zwischen 0.02 und 0.05
Annahmen.x_u = Ergebnisse_Fluegel.Fluegeltiefen_eta_oR .* Annahmen.xu_l;      % Annahme, dass l die Fluegeltiefen an der jeweiligen Position auf dem Fluegel sind
    % induzierter Widerstand des Flügels
%c_A_F = linspace(0,1); % Ergebnisse_stat_Flaechenbelastung.C_A_CR;       % aus PS4 Formel 11

    % Transsonischer Widerstand
Annahmen.Ma_unendlich = specs.Ma_CR;

    % Rumpfwiderstand
Annahmen.S_G_Ru = Anteile_einzel_Massen_FE2.Airplane_Structure.Fuselage_group.Sg_12;
Annahmen.c_A_alpha_F = VWA.c_AF_anstieg;      % Annahme bitte ueberpruefen !!!!!!!!!!!!!!!!!!!!


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
Annahmen.xu_l_Py = 0.035;

Annahmen.x_u_HLW = HLW.Fluegeltiefen_eta_oR * Annahmen.xu_l_HLW;
Annahmen.x_u_SLW = SLW.Fluegeltiefen_eta_oR * Annahmen.xu_l_SLW;
Annahmen.x_u_Py = specs.lh_TW * Annahmen.xu_l_Py;

        % Trimmwiderstand
Annahmen.l_mue = Ergebnisse_Fluegel.l_mue;  

Annahmen.qH_q = 0.85; % hat Kristof gesagt urspruenglich mit 0.95 angenommen 

    % Zusatzwidertsand


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

Annahmen.stuetzstellen = stuetzstellen;
Annahmen.kappa = 1.4;



% ------------------------------------------------------------------------
%% Fluegelwiderstand nach Diederich---------------------------
%  ------------------------------------------------------------------------


%% Funktionen direktzugriff


%PS4 S.3, Formel 7
FUN.Re_H_fun = @(l_Re,v_Re,hoehe) (l_Re .* v_Re) ./ (ISA.kin_visk(hoehe));

% PS4 S.2, Formel 6
FUN.c_f_la_fun = @(Re) 1.328./(sqrt(Re));
FUN.c_f_tu_fun = @(Re) 0.455./(log10(Re).^(2.58));

% PS4 S.4, Formel 15
FUN.tau_fun = @(Streckung, lambda) 1 - Streckung * (0.002 + 0.0084 * (lambda - 0.2)^2);

% PS 9 Formel 25 FE1 Berechnung von Auftriebsbeiwertverteilung
FUN.c_a_eta_fun = @(c_A_F) (GRA.gamma_a_eta .* c_A_F .* GRA.l_m) ./ (Ergebnisse_Fluegel.Fluegeltiefen_eta);


save Getroffene_Annahmen_und_FUN.mat Annahmen FUN




v_air = ones(stuetzstellen,1) .* specs.Ma_CR .* ISA.a(Annahmen.hoehe_CR);

%-------------------- Off_D

v_air_off_D = Ma_off_D .* ISA.a(Annahmen.hoehe_CR);

c_A_F_off_D = (((2)./(Annahmen.kappa .* ISA.p(Annahmen.hoehe_CR) .* Ma_off_D.^2)) .* Ergebnisse_stat_Flaechenbelastung.Fleachenbelastung).'; 


%% Berechnungen

% Leitwerke
[c_w_HLW_min, c_w_SLW_min] = Leitwerke_W(v_air, Annahmen.hoehe_CR);

c_W_HLW = trapz(c_w_HLW_min.').*10.^(-3);
c_W_SLW = trapz(c_w_SLW_min.').*10.^(-3);

[c_w_HLW_min_off_D, c_w_SLW_min_off_D] = Leitwerke_W(v_air_off_D, Annahmen.hoehe_CR);

c_W_HLW_off_D = trapz(c_w_HLW_min_off_D.').*10.^(-3);
c_W_SLW_off_D = trapz(c_w_SLW_min_off_D.').*10.^(-3);

% Interferenz

c_w_int_fs = Interferenz_W(v_air, Annahmen.hoehe_CR).';
c_w_int_fs_off_D = Interferenz_W(v_air_off_D, Annahmen.hoehe_CR).';

% Leitwerk Trim

[c_A_H, c_A_ges, c_w_trim] = Leitwerke_Trim_W(c_A_F);
[c_A_H_off_D, c_A_ges_off_D, c_w_trim_off_D] = Leitwerke_Trim_W(c_A_F_off_D);

% Rumpf
% PS4 S.6 Formel 25 Abwindfaktor = delta_alpha_w/delta_alpha_oH 
% Abwindfaktor = 1.75 * ((Annahmen.c_A_alpha_F)/(pi * Ergebnisse_Fluegel.streckung_phi25_max *...
%     (Ergebnisse_Fluegel.lambda * (HLW.r/(Ergebnisse_Fluegel.b/2))^(0.25) *...
%     (1+ (abs(Annahmen.z_abstand))/((Ergebnisse_Fluegel.b/2))) )));
Abwindfaktor = 1.75 * (Annahmen.c_A_alpha_F/(pi * Ergebnisse_Fluegel.streckung_phi25_max *...
    (Ergebnisse_Fluegel.lambda * (HLW.r/(Ergebnisse_Fluegel.b/2)))^0.25 *...
    (1+ (abs(Annahmen.z_abstand/(Ergebnisse_Fluegel.b/2))))));


[c_w_R_interm, alpha_Rumpf_grad_interm, c_A_alpha] = Rumpfwiderstand(specs.Ma_CR, Abwindfaktor, c_A_ges, v_air, Annahmen.hoehe_CR);
c_w_R = diag(c_w_R_interm).';
alpha_Rumpf_grad = alpha_Rumpf_grad_interm;



[c_w_R_off_D_interm, alpha_Rumpf_grad_off_D_interm, c_A_alpha_off_D] = Rumpfwiderstand(Ma_off_D, Abwindfaktor, c_A_ges_off_D, v_air_off_D, Annahmen.hoehe_CR);
c_w_R_off_D = diag(c_w_R_off_D_interm).';
alpha_Rumpf_grad_off_D = diag(alpha_Rumpf_grad_off_D_interm).';


% Triebwerke
c_w_TW = Triebwerkswiderstand(v_air, alpha_Rumpf_grad, Annahmen.hoehe_CR);
c_w_TW_off_D = Triebwerkswiderstand(v_air_off_D, alpha_Rumpf_grad_off_D, Annahmen.hoehe_CR);

% Zusatzwiderstand

[delta_c_w_H] = Zusatz_W(c_A_F, Abwindfaktor, c_A_H);
[delta_c_w_H_off_D] = Zusatz_W(c_A_F_off_D, Abwindfaktor, c_A_H_off_D);

% Profilwiderstand

[c_w_p, c_w_p_min_Re, c_w_p_test] = Profilwiderstand(v_air,c_A_F, Annahmen.hoehe_CR);
% c_w_p = c_w_p_interm; %trapz(c_w_p_interm.');
[c_w_p_off_D, c_w_p_min_Re_off_D, c_w_p_test_off_D] = Profilwiderstand(v_air_off_D, c_A_F_off_D, Annahmen.hoehe_CR);
% c_w_p_off_D = c_w_p_off_D_interm; %trapz(c_w_p_off_D_interm.');

% Induzierter Widerstand

c_w_ind = Induzierter_W(c_A_F);
c_w_ind_off_D = Induzierter_W(c_A_F_off_D);


% Wellenwiderstand / transsonischer Widerstand

[delta_c_WM, delta_Ma] = Transsonischer_W(specs.Ma_CR, c_A_F);

[delta_c_WM_off_D_interm, delta_Ma_off_D_interm] = Transsonischer_W(Ma_off_D, c_A_F_off_D);
delta_c_WM_off_D = diag(delta_c_WM_off_D_interm).';
delta_Ma_off_D = diag(delta_Ma_off_D_interm).';




%% Ergebnisse berechnen



% Anzahl der Plots festlegen
x_vector = [c_W_SLW; c_W_HLW; c_w_int_fs; c_w_R; c_w_TW; c_w_trim; delta_c_w_H; c_w_p; c_w_ind; delta_c_WM];
sz = size(x_vector);
numPlots = sz(1,1); %6;     % muss veraendert werden um off Design noch zu plotten


for n_vec = 1:numPlots
    if n_vec == 1
        x_vector_sum(n_vec,:) = x_vector(n_vec,:);
    else
        x_vector_sum(n_vec,:) = x_vector_sum(n_vec-1,:) + x_vector(n_vec,:);
    end
end

% Offdesign Vector
off_D_vector = [c_W_SLW_off_D; c_W_HLW_off_D; c_w_int_fs_off_D; c_w_R_off_D; c_w_TW_off_D; c_w_trim_off_D; delta_c_w_H_off_D; c_w_p_off_D; c_w_ind_off_D; delta_c_WM_off_D];

for n_vec_off_D = 1:numPlots
    if n_vec_off_D == 1
        x_vector_sum_off_D(n_vec_off_D,:) = off_D_vector(n_vec_off_D,:);
    else
        x_vector_sum_off_D(n_vec_off_D,:) = x_vector_sum_off_D(n_vec_off_D-1,:) + off_D_vector(n_vec_off_D,:);
    end
end

% Reziproke Gleitzahlen

Gleitverhaeltnis_Des = c_A_F ./ x_vector_sum(numPlots,:);

Gleitverhaeltnis_off_D = c_A_F_off_D ./ x_vector_sum_off_D(numPlots,:);



%% Ergebnisse speichern
schnittpunkt_off_D_c_A_CR = InterX([x_vector_sum_off_D(numPlots,:); c_A_F_off_D], [[0, 1]; [Ergebnisse_stat_Flaechenbelastung.C_A_CR, Ergebnisse_stat_Flaechenbelastung.C_A_CR]]);
gelitzahl = schnittpunkt_off_D_c_A_CR(1,1)./schnittpunkt_off_D_c_A_CR(2,1);


Ergebnisse_Widerstand_FE2.c_A_F = c_A_F;
Ergebnisse_Widerstand_FE2.c_A_F_off_D = c_A_F_off_D;
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
Ergebnisse_Widerstand_FE2.v_air = v_air;
Ergebnisse_Widerstand_FE2.v_air_off_D = v_air_off_D;
Ergebnisse_Widerstand_FE2.Ma_off_D = Ma_off_D;
Ergebnisse_Widerstand_FE2.cW_cA_off_D = gelitzahl;
Ergebnisse_Widerstand_FE2.Abwindfaktor = Abwindfaktor;
Ergebnisse_Widerstand_FE2.c_A_alpha = c_A_alpha;
Ergebnisse_Widerstand_FE2.stuetzstellen = stuetzstellen;
Ergebnisse_Widerstand_FE2.Gleitverhaeltnis_Des = Gleitverhaeltnis_Des;
Ergebnisse_Widerstand_FE2.Gleitverhaeltnis_off_D = Gleitverhaeltnis_off_D;



save Ergebnisse_Widerstand_FE2.mat Ergebnisse_Widerstand_FE2;





