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
