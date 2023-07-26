%% Fluleistung 2 PS 8 Nutzlast Reichweiten Diagramm

clc
clear all
close all


load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Basis_stat_m.mat;
load Ergebnisse_Widerstand_FE2.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Start_Landeanforderungen.mat;
load Ergebnisse_Massen_FE2.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_FLugleistung_1.mat;
addpath('Unterfunktionen Widerstand');



%% Annahmen

% geopotentialhoehe


hoehe_m = Ergebnisse_Flugleistung_1.hoehe_m;
hoehe_CR = round(unitsratio('m','ft')*(specs.flight_level*10^2));
v_CR = specs.Ma_CR * ISA.a(hoehe_CR) 
TO_Masse =  Ergebnisse_Massen_FE2.M_TO;
G_TO = TO_Masse * specs.g;
k_CR = 0.98;
S0 = k_CR * G_TO * (Ergebnisse_Widerstand_FE2.cW_cA_off_D) / (schub_CR.S_S0_CR * schub_CR.S_S0_E);



%--------------------------------------------------------------------------
%% Rechnungen
%--------------------------------------------------------------------------

% Berechnung, wo H_CR in unserem vector hoehe_m liegt

for zv1 = 1 : length(hoehe_m)
    if hoehe_m(1,zv1) <= hoehe_CR
        punkt_H_CR = zv1 + 1;
    else
    end
end


%% Kann FEHLER erzeugen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_CR_index = find(round(Ergebnisse_Flugleistung_1.v_TAS_j_CL(punkt_H_CR,:))==round(v_CR+0.2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_TAS_j_CL_index = 1;

for zv2 = Ergebnisse_Flugleistung_1.Hochpunkte.SEP_CL_x : length(Ergebnisse_Flugleistung_1.v_TAS_j_CL)
    if Ergebnisse_Flugleistung_1.v_TAS_j_CL(punkt_H_CR,zv2) <= v_CR
       v_TAS_j_CL_index = zv2+1;
    elseif v_TAS_j_CL_index > length(Ergebnisse_Flugleistung_1.v_TAS_j_CL)
       disp('Achtung zv2 ist groesser als die laenge des Vektors v_TAS_j_CL')
    end
end



%% Energiehoehe HE fuer Steigzeit


% PS8 Formel 2


HE_SEP_max_CL = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,3) + (Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,1).^2)/(2* specs.g);
HE_CL_max_vector_SEP = HE_SEP_max_CL;
% HE_CL_max_vector_SEP(length(HE_CL_max_vector_SEP)+1) = HE_SEP_max_CL(length(HE_SEP_max_CL),1);

%plot(HE_SEP_max_CL, 1./Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(:,2),'*k')

SEP_CL_max_vector = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,2); %./10;
% SEP_CL_max_vector(length(SEP_CL_max_vector)+1) =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(length(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec),2);%./10;


%t = trapz(SEP_CL_max_vector, HE_CL_max_vector_SEP)
t = trapz(HE_CL_max_vector_SEP, 1./SEP_CL_max_vector) % Ergebnis macht mit unseren inkorrekten werten sinn denke ich



%% Steigstrecke

% PS8 Formel 2

% Da wir nur sie Kostenoptimale geschwindigkeit betrachen

% HE_SET_max_CL = Ergebnisse_Flugleistung_1.TAS_SET_H_CL_vec(1:punkt_H_CR,3) + (Ergebnisse_Flugleistung_1.TAS_SET_H_CL_vec(1:punkt_H_CR,1).^2)/(2* specs.g);
% HE_CL_max_vector_SET = HE_SET_max_CL;
% HE_CL_max_vector_SET(length(HE_CL_max_vector_SET)+1) = HE_SEP_max_CL(length(HE_SEP_max_CL),1);

%plot(HE_SEP_max_CL, 1./Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(:,2),'*k')

% da wir am kostenguenstigsten steigen wollen wir aus SEP SET errechnet
SEP_CL_max_vector_R = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,2) ./ Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,1); 
% SET_CL_max_vector(length(SET_CL_max_vector)+1) =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(length(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec),2);


%t = trapz(SET_CL_max_vector, HE_CL_max_vector_SET)
R_CL = trapz(HE_CL_max_vector_SEP, 1./SEP_CL_max_vector_R) % Ergebnis macht mit unseren inkorrekten werten sinn denke ich

%% Steigkraftstoffverbrauch

[rho_rho0_H, T_H, a_H, p_p0_H, rho] = Atmos_H(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,3));


Ma_h_CL =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,1) ./ a_H;
 
[~,~,b_s_CL_1_s] = SFC_vec(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,3), Ma_h_CL, specs.bypass);
b_s_CL = b_s_CL_1_s .* (1/specs.g);
%b_s_CL_ = b_s_CL(:,1);
% b_s_CL_kg_Ns_vec(9,1) = b_s_CL_kg_Ns(length(b_s_CL_kg_Ns),1);

[S_S0_CL] = S_S0_KF_j(0.9, rho_rho0_H, Ma_h_CL, p_p0_H, specs.bypass);
S_S0_CL_vec = S_S0_CL(:,1);
% S_S0_CL_vec = S_S0_CL(length(S_S0_CL),1);
zu_integrierende_werte_Steigkraftstoff = (b_s_CL .* S_S0_CL_vec .* (S0)) ./ (SEP_CL_max_vector);

m_F_CL = trapz(HE_CL_max_vector_SEP, zu_integrierende_werte_Steigkraftstoff)
% trapz( zu_integrierende_werte_Steigkraftstoff,HE_CL_max_vector_SEP)




%% NRD
m_OE = Ergebnisse_Massen_FE2.M_OE;
m_TO = Ergebnisse_Massen_FE2.M_TO;
m_fuel = Ergebnisse_Massen_FE2.M_F;
m_ZF = Ergebnisse_Massen_FE2.M_ZF;
% Bestimmung der neuen Massen mf3, mf4, mf9

% m_3 = Ergebnisse_Massen_FE2.M_TO * FF.mf2;
m_alt_u_HLD = (1 - FF.Mff_6_10) * Ergebnisse_Massen_FE2.M_TO;
% mf3_neu = 1 - m_F_CL/m_3;
m_P = specs.m_cargo + specs.m_pax;
% Punkt A--------------------------------
m_TO_A = Ergebnisse_Massen_FE2.M_OE + m_P + m_alt_u_HLD;
m_F_HLD_A = (1 - FF.mf9) * m_TO_A;
R_A = 0;
m_RF_A = m_alt_u_HLD;
m_F_A = m_RF_A;
m_TF_A = 0;
A = [0; m_P; m_F_A; m_TF_A; m_RF_A];

% m_TO_A = m_P + m_OE + Fuel.m_F_ALT + Fuel.m_F_HLD;
% m_F_HLD_A = (1-Fuel.Fractions.m_f9) * m_TO_A;
% R_A = 0;
% m_RF_A = Fuel.m_F_ALT + m_F_HLD_A;
% m_F_A = m_RF_A;
% m_TF_A = 0;
% A = [0;m_P_alleco; m_F_A; m_TF_A; m_RF_A];

% B--------------------------------------


m_P_B = m_P;
m_TO_B = Ergebnisse_Massen_FE2.M_TO;
m_F_B = m_TO_B -Ergebnisse_Massen_FE2.M_OE - m_P_B;

m_ZF_B = m_P_B + m_OE;
mf4_B = (1.05 - m_F_B/m_TO)/(FF.mf2 * FF.mf3 * FF.mf5 * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9 + 0.05);
R_CR_B = v_CR / (FF.sfc_CR_1PERs * Ergebnisse_Widerstand_FE2.cW_cA_off_D * specs.g) * log(1/mf4_B);
R_B = (R_CL + R_CR_B)/1000;

m_F_C_B = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_B * FF.mf5)) * m_TO_B;
m_RF_OC_B = m_ZF_B * (1/((1-m_F_C_B/m_ZF_B) * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9) -1);

m_RF_B = m_RF_OC_B + m_F_C_B;
m_TF_B = m_F_B - m_RF_B;

B = [R_B; m_P_B; m_F_B; m_TF_B; m_RF_B];

% m_P_B = m_P_alleco;

% m_TO_B = m_TO;
% m_F_B = m_TO_B - m_OE - m_P_alleco;
% m_ZF_B = m_P_B + m_OE;
% 
% Fuel.Fractions.m_f4_B = (1.05 - m_F_B/m_TO)/(Fuel.Fractions.m_f2 * Fuel.Fractions.m_f3_new * Fuel.Fractions.m_f5 * Fuel.Fractions.m_f6 * Fuel.Fractions.m_f7 * Fuel.Fractions.m_f8 * Fuel.Fractions.m_f9 + 0.05);
% R_CR_B = v_CR/(sfc_CR * Widerstand.epsilon_CR*g0) * log(1/Fuel.Fractions.m_f4_B);
% R_B = (R_CL + R_CR_B)/1000;
% 
% m_F_C_B = 0.05 *(1 - (Fuel.Fractions.m_f2 * Fuel.Fractions.m_f3_new * Fuel.Fractions.m_f4_B * Fuel.Fractions.m_f5)) * m_TO_B;
% m_RF_OC_B = m_ZF_B * (1/((1-m_F_C_B/m_ZF_B) * Fuel.Fractions.m_f6 * Fuel.Fractions.m_f7 * Fuel.Fractions.m_f8 * Fuel.Fractions.m_f9) -1);
% m_RF_B = m_RF_OC_B + m_F_C_B;
% m_TF_B = m_F_B - m_RF_B;
% 
% B = [R_B; m_P_B; m_F_B; m_TF_B; m_RF_B];

% C --------------------------------------

m_F_C = m_fuel;
m_TO_C = m_TO;
m_P_C = m_TO_C - (m_OE + m_F_C);

m_ZF_C = m_OE + m_P_C;
mf4_C = (1.05 - m_F_C/m_TO_C)/(FF.mf2 * FF.mf3 * FF.mf5 * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9 + 0.05);
R_CR_C = v_CR / (FF.sfc_CR_1PERs * Ergebnisse_Widerstand_FE2.cW_cA_off_D *specs.g) * log(1/mf4_C);
R_C = (R_CL + R_CR_C)/1000;


% m_F_C_C = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_C * FF.mf5)) * m_TO_C;
m_F_C_C = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_B * FF.mf5)) * m_TO_C;
m_RF_OC_C = m_ZF_C * (1/((1-m_F_C_C/m_ZF_C) * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9) -1);
m_RF_C = m_RF_OC_C + m_F_C_C;
m_TF_C = m_F_C - m_RF_C;

C = [R_C; m_P_C; m_F_C; m_TF_C; m_RF_C];

% m_F_C = Fluegel.Tank.M_fuel;
% m_TO_C = m_TO;
% m_P_C = m_TO_C - (m_OE + m_F_C);
% 
% m_ZF_C = m_OE + m_P_C;
% Fuel.Fractions.m_f4_C = (1.05 - m_F_C/m_TO_C)/(Fuel.Fractions.m_f2 * Fuel.Fractions.m_f3_new * Fuel.Fractions.m_f5 * Fuel.Fractions.m_f6 * Fuel.Fractions.m_f7 * Fuel.Fractions.m_f8 * Fuel.Fractions.m_f9 + 0.05);
% R_CR_C = v_CR/(sfc_CR * Widerstand.epsilon_CR*g0) * log(1/Fuel.Fractions.m_f4_C);
% R_C = (R_CL + R_CR_C)/1000;
% 
% 
% m_F_C_C = 0.05 *(1 - (Fuel.Fractions.m_f2 * Fuel.Fractions.m_f3_new * Fuel.Fractions.m_f4_B * Fuel.Fractions.m_f5)) * m_TO_C;
% m_RF_OC_C = m_ZF_C * (1/((1-m_F_C_C/m_ZF_C) * Fuel.Fractions.m_f6 * Fuel.Fractions.m_f7 * Fuel.Fractions.m_f8 * Fuel.Fractions.m_f9) -1);
% m_RF_C = m_RF_OC_C + m_F_C_C;
% m_TF_C = m_F_C - m_RF_C;
% 
% C = [R_C; m_P_C; m_F_C; m_TF_C; m_RF_C];

% D --------------------------------------

m_F_D = m_fuel;
m_P_D = 0;
m_TO_D = m_OE + m_F_D ;
m_ZF_D = m_OE + m_P_D;
mf4_D = (1.05 - m_F_D / m_TO_D)/(FF.mf2 * FF.mf3* FF.mf5 * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9 + 0.05);
R_CR_D = v_CR / (FF.sfc_CR_1PERs * Ergebnisse_Widerstand_FE2.cW_cA_off_D *specs.g) * log(1/mf4_D);
R_D = (R_CL + R_CR_D)/1000;

% m_F_C_D = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_D * FF.mf5)) * m_TO_D;
m_F_C_D = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_B * FF.mf5)) * m_TO_D;
m_RF_OC_D = m_ZF_D * (1/((1-m_F_C_D/m_ZF_D) * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9) -1);
m_RF_D = m_RF_OC_D + m_F_C_D;
m_TF_D = m_F_D - m_RF_D;
 
D = [R_D; m_P_D; m_F_D; m_TF_D; m_RF_D];

% m_F_D = Fluegel.Tank.M_fuel;
% m_P_D = 0;
% m_TO_D = m_OE + m_F_D ;
% m_ZF_D = m_OE + m_P_D;
% Fuel.Fractions.m_f4_D = (1.05 - m_F_D/m_TO_D)/(Fuel.Fractions.m_f2 * Fuel.Fractions.m_f3_new * Fuel.Fractions.m_f5 * Fuel.Fractions.m_f6 * Fuel.Fractions.m_f7 * Fuel.Fractions.m_f8 * Fuel.Fractions.m_f9 + 0.05);
% R_CR_D = v_CR/(sfc_CR * Widerstand.epsilon_CR*g0) * log(1/Fuel.Fractions.m_f4_D);
% R_D = (R_CL + R_CR_D)/1000;
% 
% 
% m_F_C_D = 0.05 *(1 - (Fuel.Fractions.m_f2 * Fuel.Fractions.m_f3_new * Fuel.Fractions.m_f4_B * Fuel.Fractions.m_f5)) * m_TO_D;
% m_RF_OC_D = m_ZF_D * (1/((1-m_F_C_D/m_ZF_D) * Fuel.Fractions.m_f6 * Fuel.Fractions.m_f7 * Fuel.Fractions.m_f8 * Fuel.Fractions.m_f9) -1);
% m_RF_D = m_RF_OC_D + m_F_C_D;
% m_TF_D = m_F_D - m_RF_D;
% 
% D = [R_D; m_P_D; m_F_D; m_TF_D; m_RF_D];



figure(1) 
hold on 
grid on 
ylim([100000 m_TO+10000])
xlim([0 20000])

p1(1) = plot([A(1) B(1) C(1) D(1)].*10,[A(2)+m_OE, B(2)+m_OE, C(2)+m_OE, D(2)+m_OE]);
p1(2) = plot([A(1) B(1) C(1) D(1)].*10,[A(3)+m_ZF, B(3)+m_ZF, C(3)+m_ZF, D(3)+m_ZF]);
p1(3) = plot([A(1) B(1) C(1) D(1)].*10,[A(4)+m_ZF, B(4)+m_ZF, C(4)+m_ZF, D(4)+m_ZF]);
p1(4) = plot([A(1) B(1) C(1) D(1)].*10,[A(5)+m_ZF, B(5)+m_ZF, C(5)+m_ZF, D(5)+m_ZF]);
p1(5) = plot([0, 20000],[m_TO, m_TO], Color=[0.5 0.5 0.5], LineStyle="--");
p1(6) = plot([0, 20000],[m_OE, m_OE], Color=[0.5 0.5 0.5], LineStyle="-.");
% % plot(R_DP, m_P,'rx')

title('Nutzlast-Reichweiten-Diagramm', 'FontSize',25)
legend(p1([1:6]),{'Nutzlast', 'Treibstoffmasse', 'Reisekraftstoffmasse', 'Reserve', 'M_{TO}', 'M_{OE}'},...
     'Location','southwest','FontSize',25);
xlabel('Reichweite in km','FontSize',20)
ylabel('Masse in kg','FontSize',20)


%% Funktionen


function [rho_rho0_H, T_H, a_H, p_p0_H, rho] = Atmos_H(hoehe_m)
    load Ergebnisse_ISA_DATA.mat
    for zv = 1:size(hoehe_m)
    rho_rho0_H(zv,1) = ISA.rho(hoehe_m(zv,1)) / ISA.rho_0;
    T_H(zv,1) = ISA.T(hoehe_m(zv,1));
    a_H(zv,1)= ISA.a(hoehe_m(zv,1));
    p_p0_H(zv,1) = ISA.p(hoehe_m(zv,1)) / ISA.p0;
    rho(zv,1) = ISA.rho(hoehe_m(zv,1));
    end
end



%% Funktion zur berechnung von Gleichung 12 PS7 s.3 
function [S_S0] = S_S0_KF_j(D, rho_rho0_H, Ma_j, p_p0, bypass)
    load Ergebnisse_FLugleistung_1.mat
    S_S0_KF_j = D .* rho_rho0_H .* exp(-0.35 .* Ma_j .* p_p0 .* sqrt(bypass));
    S_S0_E = 1 - (1.3 + 0.25 * bypass) * 0.02;
    S_S0 = S_S0_KF_j .* S_S0_E;
end


