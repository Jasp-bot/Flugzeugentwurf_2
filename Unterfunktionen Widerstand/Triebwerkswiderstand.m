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
c_w_TW_min = c_f_tu_TW * (1 + k_TW) * (Annahmen.S_G_TW)/Ergebnisse_Fluegel.F;

% PS4 S.5 Formel 22
c_w_TW_zu_c_w_TWmin = 0.000208 * abs(alpha_TW_grad).^3 + 0.00125 * abs(alpha_TW_grad).^2 + 0.029 * abs(alpha_TW_grad) + 1;


% PYLON
% PS4 S.3 Formel 7 angewendet auf Leitwerke
Re_u_Py = FUN.Re_CR_fun(Annahmen.x_u_Py, v_air);  %Annahmen.x_u_HLW * v_air / ISA.kin_visk(Annahmen.hoehe_CR);

Re_Py = FUN.Re_CR_fun(specs.lh_TW, v_air);    % HLW.Fluegeltiefen_eta_oR .* v_air ./ ISA.kin_visk(Annahmen.hoehe_CR);

% PS4 S.2 Formel 6 angewendet auf Leitwerke
c_f_la_xu_Py = FUN.c_f_la_fun(Re_u_Py); % (1.328)./(sqrt(Re_u_HLW));
c_f_tu_xu_Py = FUN.c_f_tu_fun(Re_u_Py); % 0.455./(log(Re_u_HLW).^(2.58));
c_f_tu_l_Py = FUN.c_f_tu_fun(Re_Py); % 0.455./(log(Re_HLW).^(2.58));

% PS4 S.2 Formel 5 angewendet auf Pylon
c_f_Py = c_f_tu_l_Py - Annahmen.xu_l_Py .* (c_f_tu_xu_Py - c_f_la_xu_Py);

% PS4 S.6 Formel 28 % Annahme 
% Würde ich jetzt mal als das selbe Profil annehmen
k_Py = 2.7 .* Annahmen.d_l_HLW + 100 .* Annahmen.d_l_HLW.^4;

% PS4 S.6 Formel 27
c_w_Py_min = 2 .* c_f_Py .* (1+ k_Py .* cosd(specs.pylonangle_deg).^2) .* ((specs.F_pylon)/(Ergebnisse_Fluegel.F));


c_w_TW = diag(c_w_TW_zu_c_w_TWmin .* c_w_TW_min .* 2 ).' + (2 * c_w_Py_min).'; % 2x Pylons
end
