%% Funktion Induzierter Widerstand

function [c_w_ind] = Induzierter_W(c_A_F)

load Projekt_specs.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Auftrieb_Momente.mat;
load Getroffene_Annahmen_und_FUN.mat;


% PS4 S.4, Formel 12
c0 = VWA.c_AF_anstieg^2 * (0.0088 * Ergebnisse_Fluegel.lambda - 0.0051 * Ergebnisse_Fluegel.lambda^2) * (1 - 0.0006 * Ergebnisse_Fluegel.streckung_phi25_max^2);

% PS4 S.4, Formel 13
% c1 = GRA.c_a_anstieg * (0.0134 * (Ergebnisse_Fluegel.lambda - 0.3) - 0.0037 * Ergebnisse_Fluegel.lambda^2);
c1 = VWA.c_AF_anstieg * (0.0134 * (Ergebnisse_Fluegel.lambda - 0.3) - 0.0037 * Ergebnisse_Fluegel.lambda^2); %% nicht sicer mit welchem Auftriebsanstieg gerechnet werden muss
% PS4 S.4, Formel 15
%tau = 1 - Ergebnisse_Fluegel.streckung_phi25_max * (0.002 + 0.0084 * (Ergebnisse_Fluegel.lambda-0.2)^2);
tau = FUN.tau_fun(Ergebnisse_Fluegel.streckung_phi25_max, Ergebnisse_Fluegel.lambda);

% PS4 S.4, Formel 14
c2 = (1/tau) * (1 + (5*10^(-6)) * (rad2deg(Ergebnisse_Fluegel.phi_25_max))^3); 

% delta eps keine ahnung wie das berechnet werden soll
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% achtung nachfragen!!!!!!!!!!!!!!!!
eta_Ru = length(VWA.epsilon_eta_Ru)*10^(-3);
delta_eps = abs(VWA.epsilon - VWA.epsilon .* eta_Ru);

% PS4 S.4, Formel 11
c_w_ind = c2 .* (c_A_F.^2)./(pi * Ergebnisse_Fluegel.streckung_phi25_max) +...
    c1 .* c_A_F .* delta_eps + c0 .* delta_eps^2;

end
