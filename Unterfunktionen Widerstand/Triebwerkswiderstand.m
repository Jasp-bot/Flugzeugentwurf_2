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
