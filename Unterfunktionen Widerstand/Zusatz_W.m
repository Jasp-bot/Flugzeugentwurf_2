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
