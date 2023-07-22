% Funktion Trimwiderstand Leitwerke
% [c_A_H, c_A_ges, c_w_trim] = Leitwerke_Trim_W(c_A_F)

function [c_A_H, c_A_ges, c_w_trim] = Leitwerke_Trim_W(c_A_F)


load Projekt_specs.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Leitwerke.mat;
load Getroffene_Annahmen_und_FUN.mat;


% PS4 S.7 Formel 33 

% c_A_H(1,n_iteration) = (Annahmen.c_M_0_F + c_A_F * ((Annahmen.x_SP_MAC)/(Annahmen.l_mue))) / ...
%     (Annahmen.qH_q * ((HLW.F)/(Ergebnisse_Fluegel.F)) * ((HLW.r - Annahmen.x_SP_MAC)/(Annahmen.l_mue)));

c_A_H = (Annahmen.c_M_0_F + c_A_F .* Annahmen.dx_SP_lmue) ./ ...
                        (Annahmen.qH_q .* ((HLW.F)./(Ergebnisse_Fluegel.F)) .*...
                        (((HLW.r./Annahmen.l_mue) - Annahmen.dx_SP_lmue)));


% PS4 S.7 Formel 29
c_A_ges = c_A_F + c_A_H * Annahmen.qH_q * ((HLW.F)/(Ergebnisse_Fluegel.F));
% c_A_ges_vec(Annahmen.zaehlvariabele_itt,1) = c_A_ges(1,n_iteration);

Annahmen.zaehlvariabele_itt = Annahmen.zaehlvariabele_itt + 1;


% PS4 S.4 Formel 15
tau_H = FUN.tau_fun(HLW.streckung_phi25, HLW.lambda); %1 - HLW.streckung_phi25 * (0.002 + 0.0084 * (HLW.lambda - 0.2)^2);

% PS4 S.7 Formel 30
% c_w_trim = ((c_A_H.^2) ./ (pi * HLW.streckung_phi25)) .* ...
%     ((1 + (5*10^(-6)) .* (abs(rad2deg(HLW.phi_25))).^3)./(tau_H)) .*...
%     ((HLW.F)./(Ergebnisse_Fluegel.F));

c_w_trim = ((c_A_H.^2)/(pi * HLW.streckung_phi25)) .* ((1 + 5*10^(-6) * (abs(rad2deg(HLW.phi_25))^3))./(tau_H)) .* (HLW.F/Ergebnisse_Fluegel.F); 

end
