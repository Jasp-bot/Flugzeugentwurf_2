% PS9 Ergebnisse Auftriebsverteilung, Momentenbeiwert, Leitwerke

% function Berechnung_PS9_Auftrieb_Momente

clc
clear all
close all

%% Laden der Werte

load Projekt_specs.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load CaVerteilungFE12.mat;
load c_m_NP_verteilung_engauge.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat

%% Grundrissabhaengiger Anteil

% Formel 5 PS 9
% gamma_a(eta) = c1 *((l(eta)/l_m)+ c2*(3*sqrt(1-eta^2))/(pi)) + c3 * f(eta)



eta_SP = [NP.y_SP_A; NP.y_SP_I; NP.y_SP_R; NP.y_SP_ges];

% Berechung für c1 c2 c3

GRA.Ma_P_CR = specs.Ma_CR * sqrt(cos(Ergebnisse_Fluegel.phi_25_max)); % profilmachzahl zum Ablesen der c_a_anstiegs werte

GRA.c_a_anstieg = (CaVerteilungFE12(10,2)-CaVerteilungFE12(2,2))/...
    (deg2rad(CaVerteilungFE12(10,1))-deg2rad(CaVerteilungFE12(2,1))); % aus Profilkatalog kein plan ob richtig

FF = Ergebnisse_Fluegel.streckung_phi25_max / ((GRA.c_a_anstieg/(2*pi)) * cos(Ergebnisse_Fluegel.phi_25_max));

c1 = 0.003647 + FF * (0.05614 - 0.0009842 * FF);
c2 = 1 + FF *(0.002721 * FF - 0.1098);
c3 = -0.008265 + FF * (0.05493 - 0.001816 * FF);


% Pfeilungskorrektur f(eta)
winkel = ([-45 -30 -15 0 15 30 45 60]);  %% kein plan ob richtig

koeff_pfeilung = [2.1511 1.7988 1.53045 1.2621 1.08715 0.9122 0.6533 0.559;
    0.009 0.4009 0.2337 0.0665 0.1467 0.2209 1.3704 0.8908;
    -7.3509 -7.887 -3.75395 0.3791 4.39155 8.404 6.9543 11.602;
    7.3094 17.856 5.9553 -5.9454 -16.4117 -26.878 -25.579 -37.973;
    -2.3048 -23.375 -6.447 10.481 20.644 30.807 29.78 44.477;
    -0.4104 10.867 2.31195 -6.2431 -9.88955 -13.536 -13.249 -19.627];


phi_eff = atan((tan(Ergebnisse_Fluegel.phi_25_max))/(sqrt(1-specs.Ma_CR^2)));
% rad2deg(phi_eff)
eta = (0:0.001:1);
startspalte =5 ;
endspalte =6 ;
% Berechnung von linearer Interpoltion
% for startwert = 1:6
%      lin_interp(1,startwert) = koeff_pfeilung(startwert,startspalte) + (koeff_pfeilung(startwert, endspalte)-koeff_pfeilung(startwert, startspalte)/(winkel(1,endspalte)- winkel(1,startspalte))) * (rad2deg(phi_eff) - winkel(1,startspalte));
%  end
% 
%  
%  % Berechnung der Summe für Pfeilungskorrektur f(eta)
% 
% for start = 1:6
%     f_eta_mat(:,start) = (lin_interp(1,start) .*(eta.').^(start-1));
%     %f_eta_coeff(:,start = )
% end


% Berechnung von linearer Interpoltion
for startwert = 1:6
    inter_coef(1,startwert) =  interp1(winkel, koeff_pfeilung(startwert,:),rad2deg(phi_eff));
end

f_eta = 0;
for n = 1:6
    f_eta1 =  inter_coef(1,n) .* ((eta.')).^(n-1);
    f_eta = f_eta + f_eta1;
end

% plot(eta,f_eta);
% plot(CaVerteilungFE12(:,1),CaVerteilungFE12(:,2));

GRA.l_m = Ergebnisse_stat_Flaechenbelastung.F/Ergebnisse_Fluegel.b;


% !!!!! eventuell ist einfach nur eta einsetzten nicht richtig
GRA.gamma_a_eta = c1 .* (((Ergebnisse_Fluegel.Fluegeltiefen_eta))./(GRA.l_m)) + c2 .* (4.*sqrt(1 - eta.^2))./pi + c3 .* f_eta.'; 

%% Verwindungsabhaengiger teil gamma_b(eta)

% gamma_b_eta = k1 * c_AF_anstieg * gamma_a_eta * (eps(eta) - integal(eps(eta)*gamma_a_eta(eta) d_eta)); 


k0 = 0.2692 * FF - 0.03873 * FF^2 + 0.002951 * FF^3 - 0.0001106 *FF^4 + 1.559 * 10^(-6) *FF^5;
k1 = 0.3064 + FF * (0.05185 - 0.0014 * FF);

VWA.c_AF_anstieg = GRA.c_a_anstieg * k0 * cos(Ergebnisse_Fluegel.phi_25_max);

VWA.epsilon = deg2rad(-5);% geometrische Verwindung des Fuegels
VWA.epsilon_eta = (VWA.epsilon) .* eta;
VWA.epsilon_eta_oR = VWA.epsilon_eta(1,Ergebnisse_Fluegel.zaehlvariabele_eta_Ru : length(VWA.epsilon_eta));
VWA.epsilon_eta_Ru = VWA.epsilon_eta(1,1:Ergebnisse_Fluegel.zaehlvariabele_eta_Ru);

VWA.gamma_b_fun = @(eta) ((VWA.epsilon) .* eta) * (c1 .* ((eta)./(GRA.l_m)) + c2 .* (4.*sqrt(1 - eta.^2))./pi + c3 .* f_eta.');
VWA.test_integral = integral(VWA.gamma_b_fun, 0, 1, 1001); % trapezintegral
%test = trapz(eta,gamma_b_fun)

VWA.gamma_b_eta = k1 .* VWA.c_AF_anstieg .* GRA.gamma_a_eta .*(VWA.epsilon_eta - VWA.test_integral);

%% berechnung von gamma(eta)
c_AF = Ergebnisse_stat_Flaechenbelastung.C_A_CR;
gamma_eta_ges = GRA.gamma_a_eta * c_AF + VWA.gamma_b_eta;



c_a_eta = (GRA.gamma_a_eta .* Ergebnisse_stat_Flaechenbelastung.C_A_CR .* GRA.l_m) ./ (Ergebnisse_Fluegel.Fluegeltiefen_eta);

%% Part 2 Fluegelmomente

% Grundrissabhaengiger Anteil (Fluegelmomente)

% c_M_NP_B = (2)/(F * l_eta) * integral([0 bis b/2] * c_M_NP(y) * (l(y)^2) * d_y);

% fo = fit(CMNPBDiagramm25(:,1),CMNPBDiagramm25(:,2),'cubicinterp');
% 
% plot(CMNPBDiagramm25(:,1),CMNPBDiagramm25(:,2))

FM.c_M_NP_eta = mean(CMNPBDiagramm25(:,1));

% %fun_cMNP_B = @(eta)  c_M_NP_eta  * (Fluegeltiefen_eta./ l_m ).^2;
% % test2_int = integral(fun_cMNP_B, 0, 1, 101);

for laufvariable_C_M_NP_B = 1:1:length(Ergebnisse_Fluegel.Fluegeltiefen_eta);
    FM.teilergebnis(1,laufvariable_C_M_NP_B) = FM.c_M_NP_eta  .* (Ergebnisse_Fluegel.Fluegeltiefen_eta(1,laufvariable_C_M_NP_B)./ GRA.l_m ).^2; %
end  

% summe = sum(teilergebnis);
FM.summe = trapz(FM.teilergebnis)*10^(-3);
%FM.summe = trapezoidal_area(eta,FM.teilergebnis);
FM.c_M_NP_B = (GRA.l_m/ NP.l_mue_ges) .* FM.summe; %integral(fun_cMNP_B, 0, 1, 101);


% Verwindungsanteil

% fun_delta_c_M_NP_eps = @(x) gamma_b_eta .* x;
for laufvar_delta_C_M_NP_ets = 1:1:length(Ergebnisse_Fluegel.Fluegeltiefen_eta)
    FM.teilergebnis2(1,laufvar_delta_C_M_NP_ets) = VWA.gamma_b_eta(1,laufvar_delta_C_M_NP_ets) * (laufvar_delta_C_M_NP_ets-1)*10^(-3) ; 
end
% summe2 = sum(teilergebnis2);
FM.summe2 =  trapz(FM.teilergebnis2)*10^(-4); %%%%%%%%%%%%%%%%%%%%%%%%%%%%% eigentlich 10^(-3)
%FM.summe2 = trapezoidal_area(FM.teilergebnis2,eta); % trapz(FM.teilergebnis2);%*10^(-3);
FM.delta_c_M_NP_eps = - ((Ergebnisse_Fluegel.streckung_phi25_max .* tan(Ergebnisse_Fluegel.phi_25_max) .* GRA.l_m)./(2 .* NP.l_mue_ges)) .* FM.summe2; %integral(fun_delta_c_M_NP_eps, 0, 1, 101); 


% Nullmomentenbeiwert

FM.c_M_NP_F0 = FM.c_M_NP_B + FM.delta_c_M_NP_eps;







%% Ergebnisse der Auftriebsverteilung speichern

Ergebnisse_Auftriebsverteilung.gamma_eta_ges = gamma_eta_ges;
Ergebnisse_Auftriebsverteilung.gamma_b_eta = VWA.gamma_b_eta;
Ergebnisse_Auftriebsverteilung.gamma_a_eta = GRA.gamma_a_eta;
Ergebnisse_Auftriebsverteilung.VWA = VWA;
Ergebnisse_Auftriebsverteilung.GRA = GRA;
Ergebnisse_Auftriebsverteilung.eta = eta;
Ergebnisse_Auftriebsverteilung.c_a_eta = c_a_eta;
save Ergebnisse_Auftrieb_Momente.mat Ergebnisse_Auftriebsverteilung VWA GRA FM












