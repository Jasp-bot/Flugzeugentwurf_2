%% Funktion Rumpfwiderstand

%function [c_w_R, alpha_Rumpf_grad, c_A_alpha] = Rumpfwiderstand(Machzahl, Abwindfaktor, c_A_ges, v_air)


% clc
% clear all



load Projekt_specs.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Auftrieb_Momente.mat;
load Ergebnisse_Leitwerke.mat;
load Getroffene_Annahmen_und_FUN.mat;
load Ergebnisse_Widerstand_FE2.mat;
load Ergebnisse_ISA_DATA.mat;



Machzahl = specs.Ma_CR;
Abwindfaktor = Ergebnisse_Widerstand_FE2.Abwindfaktor;
c_A_ges = Ergebnisse_Widerstand_FE2.c_A_ges;
v_air = Ergebnisse_Widerstand_FE2.v_air;

%Re_Ru = (specs.l_rumpf * (specs.Ma_CR * ISA.a(Annahmen.hoehe_CR)))/(ISA.kin_visk(Annahmen.hoehe_CR));
Re_Ru = FUN.Re_CR_fun(specs.l_rumpf, v_air);
% Re_Ru = FUN.Re_CR_fun(73.8, v_air);
% PS4 S.5 Formel 19
c_f_tu_Ru = FUN.c_f_tu_fun(Re_Ru); %(0.455)/(log(Re_Ru)^(2.58));

% PS4 S.5 Formel 21
k_Rumpf = 2.2 * (specs.D_rumpf / specs.l_rumpf)^(3/2) + 3.8 * (specs.D_rumpf/specs.l_rumpf)^(3);


% PS4 S.5 Formel 21
c_w_Ru_min = c_f_tu_Ru * (1 + k_Rumpf) * (Annahmen.S_G_Ru)/Ergebnisse_Fluegel.F;

% Kommentar: Ich habe keine Berechnungen fuer den Widerstand unter
% betrachtung des Anstellwinkes des Rumpfes gemacht, das bedeutet, dise
% Rechnung gilt nur fuer CR-Zustand Wenn das hinzugefuegt werden muss PFluegeltiefen_eta_oRS4
% S.5 und fogend

% PS4 S.6 Formel 26
c_A_alpha_H = (pi * HLW.streckung_phi25)./...
    (1 + sqrt(1 + 0.25 .* HLW.streckung_phi25.^2 * (tan(HLW.phi_50).^2 + (1 - Machzahl.^2)))); % phi_50

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
alpha_Rumpf_grad = (((c_A_ges - Ergebnisse_stat_Flaechenbelastung.C_A_CR)./(c_A_alpha)) .* (180 ./ pi));
% PS4 S.5 Formel 22
c_w_R_zu_c_w_Rmin = 0.000208 .* (abs(alpha_Rumpf_grad)).^3 + 0.00125 * (abs(alpha_Rumpf_grad)).^2 + 0.029 .* (abs(alpha_Rumpf_grad)) + 1;

c_w_R = (c_w_R_zu_c_w_Rmin .* c_w_Ru_min); %.*5;


%end
