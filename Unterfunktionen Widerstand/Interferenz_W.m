%% Funktion Interferenzwiderstand
% [c_w_int_fs] = Interferenz_W(v_air)


function [c_w_int_fs] = Interferenz_W(v_air, hoehe)

load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Getroffene_Annahmen_und_FUN.mat;

Re_F_wurzel = FUN.Re_H_fun(Annahmen.l_int_F, v_air, hoehe); % Annahmen.l_int_F * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_F = ((0.1369)./(Re_F_wurzel).^0.4) .* Annahmen.l_int_F^2 .* Annahmen.n_int_F;
        

        % HLW
Re_HLW_wurzel = FUN.Re_H_fun(Annahmen.l_int_HLW, v_air, hoehe); %  Annahmen.l_int_HLW * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_HLW = ((0.1369)./(Re_HLW_wurzel).^0.4) .* Annahmen.l_int_HLW.^2 .* Annahmen.n_int_HLW;


        % SLW
Re_SLW_wurzel = FUN.Re_H_fun(Annahmen.l_int_SLW, v_air, hoehe); %  Annahmen.l_int_SLW * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_SLW = ((0.1369)./(Re_SLW_wurzel).^0.4) .* Annahmen.l_int_SLW^2 * Annahmen.n_int_SLW;


        % TW Nacell
Re_NC_wurzel = FUN.Re_H_fun(Annahmen.l_int_NC, v_air, hoehe); %  Annahmen.l_int_NC * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_NC = ((0.1369)./(Re_NC_wurzel).^0.4) .* Annahmen.l_int_NC^2 .* Annahmen.n_int_NC;


        % Pyl Nacell
Re_PYL_wurzel = FUN.Re_H_fun(Annahmen.l_int_PYL, v_air, hoehe); %  Annahmen.l_int_PYL * Annahmen.v_air / ISA.kin_visk(Annahmen.hoehe_CR);
% PS4 S.9 Formel 40
c_w_int_PYL = ((0.1369)./(Re_PYL_wurzel).^0.4) .* Annahmen.l_int_PYL^2 .* Annahmen.n_int_PYL;

    % zusammenfassung

c_w_int_fs = (c_w_int_F + c_w_int_HLW + ...
    c_w_int_SLW +c_w_int_NC + c_w_int_PYL)./...
    Ergebnisse_Fluegel.F;

end
