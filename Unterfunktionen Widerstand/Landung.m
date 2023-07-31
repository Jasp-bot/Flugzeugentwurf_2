function [x_vector_sum, x_vector] = Landung(v_eingang, hoehe_LDG, c_A_F)
% Achtung v_eingang bitte nur als ein Skalar
% hoehe_LDG in Meter gerundet auf ganze Zahlen
% c_A_F muss ein zeilenvektor sein


load Projekt_specs.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Getroffene_Annahmen_und_FUN.mat;
load Ergebnisse_Leitwerke.mat;
load Ergebnisse_ISA_DATA.mat;
addpath('Unterfunktionen Widerstand');

v_air =  ones(length(c_A_F),1) .* v_eingang; % ones(1,1000).* 

% hoehe_LDG = round(unitsratio('m','ft') * 1500);

% c_A_F = linspace(0,3,1000);

% Leitwerke
[c_w_HLW_min, c_w_SLW_min] = Leitwerke_W(v_air, hoehe_LDG);


c_W_HLW = trapz(c_w_HLW_min.').*10.^(-3);
c_W_SLW = trapz(c_w_SLW_min.').*10.^(-3);

% Interferenz

c_w_int_fs = Interferenz_W(v_air, hoehe_LDG).';

% Leitwerk Trim

[c_A_H, c_A_ges, c_w_trim] = Leitwerke_Trim_W(c_A_F);

% Rumpf
Abwindfaktor = 1.75 * (Annahmen.c_A_alpha_F/(pi * Ergebnisse_Fluegel.streckung_phi25_max *...
    (Ergebnisse_Fluegel.lambda * (HLW.r/(Ergebnisse_Fluegel.b/2)))^0.25 *...
    (1+ (abs(Annahmen.z_abstand/(Ergebnisse_Fluegel.b/2))))));

Machzahl = v_air ./ ISA.a(hoehe_LDG);

[c_w_R_interm, alpha_Rumpf_grad_interm, c_A_alpha] = Rumpfwiderstand(Machzahl, Abwindfaktor, c_A_ges, v_air, hoehe_LDG);
c_w_R = diag(c_w_R_interm).';
alpha_Rumpf_grad = alpha_Rumpf_grad_interm;

% Triebwerke
c_w_TW = Triebwerkswiderstand(v_air, alpha_Rumpf_grad, hoehe_LDG);

% Zusatzwiderstand

[delta_c_w_H] = Zusatz_W(c_A_F, Abwindfaktor, c_A_H);

% Profilwiderstand

[c_w_p, c_w_p_min_Re, c_w_p_test] = Profilwiderstand(v_air,c_A_F, hoehe_LDG);

% Induzierter Widerstand

c_w_ind = Induzierter_W(c_A_F);

% Wellenwiderstand / transsonischer Widerstand

[delta_c_WM_mat, delta_Ma] = Transsonischer_W(Machzahl, c_A_F);
delta_c_WM = diag(delta_c_WM_mat).';
% Anzahl der Plots festlegen                                                   %c_w_p;
x_vector = [c_W_SLW; c_W_HLW; c_w_int_fs; c_w_R; c_w_TW; c_w_trim; delta_c_w_H;  c_w_ind; delta_c_WM];

sz = size(x_vector);
numPlots = sz(1,1); %6;     % muss veraendert werden um off Design noch zu plotten


for n_vec = 1:numPlots
    if n_vec == 1
        x_vector_sum(n_vec,:) = x_vector(n_vec,:);
    else
        x_vector_sum(n_vec,:) = x_vector_sum(n_vec-1,:) + x_vector(n_vec,:);
    end
end

Gleitverhaeltnis_Des = c_A_F ./ x_vector_sum(numPlots,:);


end