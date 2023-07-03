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