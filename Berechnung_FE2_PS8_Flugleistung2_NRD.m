%% Fluleistung 2 PS 8 Nutzlast Reichweiten Diagramm

%% Kommentar
% !!! Annahme, dass M_OE von ECO und 3 KLassen gleich ist



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
load Schwerpunkt.mat;
addpath('Unterfunktionen Widerstand');



%% Annahmen

% geopotentialhoehe


hoehe_m = Ergebnisse_Flugleistung_1.hoehe_m;
hoehe_CR = round(unitsratio('m','ft')*(specs.flight_level*10^2));
hoehe_ALT = round(unitsratio('m','ft')*(specs.flight_level_ALT*10^2));
hoehe_HLD = round(unitsratio('m', 'ft')*1500);
v_CR = specs.Ma_CR * ISA.a(hoehe_CR);
v_ALT = specs.Ma_CR * ISA.a(hoehe_ALT);


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

for zv3 = 1 : length(hoehe_m)
    if hoehe_m(1,zv3) <= hoehe_ALT
        punkt_H_ALT = zv3 + 1;
    else
    end
end


%% Kann FEHLER erzeugen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v_CR_index = find(round(Ergebnisse_Flugleistung_1.v_TAS_j_CL(punkt_H_CR,:))==round(v_CR+0.2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_TAS_j_CL_index = 1;
HP_x_SEP_CL_H_CR = Ergebnisse_Flugleistung_1.Hochpunkte.SEP_CL_x(punkt_H_CR,1); 
for zv2 = HP_x_SEP_CL_H_CR : length(Ergebnisse_Flugleistung_1.v_TAS_j_CL)
    if Ergebnisse_Flugleistung_1.v_TAS_j_CL(punkt_H_CR,zv2) <= v_CR
       v_TAS_j_CL_index = zv2+1;
    elseif v_TAS_j_CL_index > length(Ergebnisse_Flugleistung_1.v_TAS_j_CL)
       disp('Achtung zv2 ist groesser als die laenge des Vektors v_TAS_j_CL')
    end
end

v_TAS_j_CL_ALT_index = 1;
HP_x_SEP_CL_H_ALT = Ergebnisse_Flugleistung_1.Hochpunkte.SEP_CL_ALT_x(punkt_H_ALT,1); 
for zv4 = HP_x_SEP_CL_H_ALT : length(Ergebnisse_Flugleistung_1.v_TAS_j_CL_ALT)
    if Ergebnisse_Flugleistung_1.v_TAS_j_CL_ALT(punkt_H_CR,zv2) <= v_ALT
       v_TAS_j_CL_ALT_index = zv4+1;
    elseif v_TAS_j_CL_ALT_index > length(Ergebnisse_Flugleistung_1.v_TAS_j_CL_ALT)
       disp('Achtung zv4 ist groesser als die laenge des Vektors v_TAS_j_CL_ALT')
    end
end



%% Energiehoehe HE fuer Steigzeit bis auf CR


% PS8 Formel 2
HE_SEP_max_CL = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,3) + (Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,1).^2)/(2* specs.g);
HE_CL_max_vector_SEP = HE_SEP_max_CL;

% fuer Beschleunigungsanteil
v_beschl_CL = Ergebnisse_Flugleistung_1.v_TAS_j_CL(punkt_H_CR, HP_x_SEP_CL_H_CR:v_TAS_j_CL_index).';
HE_SEP_besch = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(punkt_H_CR,3) + (v_beschl_CL .^2)./(2* specs.g);


SEP_CL_max_vector = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,2); %./10;
% SEP_CL_max_vector(length(SEP_CL_max_vector)+1) =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(length(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec),2);%./10;
SEP_CL_max_beschl_vec = Ergebnisse_Flugleistung_1.SEP_CL(punkt_H_CR, HP_x_SEP_CL_H_CR : v_TAS_j_CL_index).';

%t = trapz(SEP_CL_max_vector, HE_CL_max_vector_SEP)
Steigflug.t_CL = trapz(HE_CL_max_vector_SEP, 1./SEP_CL_max_vector); % Ergebnis macht mit unseren inkorrekten werten sinn denke ich
Steigflug.t_beschl = trapz(HE_SEP_besch, 1./SEP_CL_max_beschl_vec);

Steigflug.t = Steigflug.t_CL + Steigflug.t_beschl;

%% Energiehoehe HE fuer Steigzeit bis auf ALT


% PS8 Formel 2
HE_SEP_max_CL_ALT = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec(1:punkt_H_ALT,3) + (Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec(1:punkt_H_ALT,1).^2)/(2* specs.g);
HE_CL_max_vector_SEP_ALT = HE_SEP_max_CL_ALT;

% fuer Beschleunigungsanteil
v_beschl_CL_ALT = Ergebnisse_Flugleistung_1.v_TAS_j_CL_ALT(punkt_H_ALT, HP_x_SEP_CL_H_ALT:v_TAS_j_CL_ALT_index).';
HE_SEP_besch_ALT = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec(punkt_H_ALT,3) + (v_beschl_CL_ALT .^2)./(2* specs.g);


SEP_CL_ALT_max_vector = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec(1:punkt_H_ALT,2); %./10;
% SEP_CL_max_vector(length(SEP_CL_max_vector)+1) =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(length(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec),2);%./10;
SEP_CL_ALT_max_beschl_vec = Ergebnisse_Flugleistung_1.SEP_CL_ALT(punkt_H_ALT, HP_x_SEP_CL_H_ALT : v_TAS_j_CL_ALT_index).';

%t = trapz(SEP_CL_max_vector, HE_CL_max_vector_SEP)
Steigflug.t_CL_ALT = trapz(HE_CL_max_vector_SEP_ALT, 1./SEP_CL_ALT_max_vector); % Ergebnis macht mit unseren inkorrekten werten sinn denke ich
Steigflug.t_beschl_ALT = trapz(HE_SEP_besch_ALT, 1./SEP_CL_ALT_max_beschl_vec);

Steigflug.t_ALT = Steigflug.t_CL_ALT + Steigflug.t_beschl_ALT;


%% Steigstrecke CL auf CR

% % PS8 Formel 2
% 
% % Da wir nur sie Kostenoptimale geschwindigkeit betrachen
% 
% % HE_SET_max_CL = Ergebnisse_Flugleistung_1.TAS_SET_H_CL_vec(1:punkt_H_CR,3) + (Ergebnisse_Flugleistung_1.TAS_SET_H_CL_vec(1:punkt_H_CR,1).^2)/(2* specs.g);
% % HE_CL_max_vector_SET = HE_SET_max_CL;
% % HE_CL_max_vector_SET(length(HE_CL_max_vector_SET)+1) = HE_SEP_max_CL(length(HE_SEP_max_CL),1);
% 
% %plot(HE_SEP_max_CL, 1./Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(:,2),'*k')

% da wir am kostenguenstigsten steigen wollen wir aus SET SEP errechnet
SEP_CL_max_vector_R = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,2) ./ Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,1); 
% SET_CL_max_vector(length(SET_CL_max_vector)+1) =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(length(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec),2);
SEP_CL_max_beschl_vector_R = SEP_CL_max_beschl_vec ./ v_beschl_CL;

%t = trapz(SET_CL_max_vector, HE_CL_max_vector_SET)
Steigflug.R_CL1 = trapz(HE_CL_max_vector_SEP, 1./SEP_CL_max_vector_R); % Ergebnis macht mit unseren inkorrekten werten sinn denke ich
Steigflug.R_beschl = trapz(HE_SEP_besch, 1./SEP_CL_max_beschl_vector_R);

Steigflug.R_CL = Steigflug.R_CL1 + Steigflug.R_beschl;

%% Steigstrecke CL auf ALT

% da wir am kostenguenstigsten steigen wollen wir aus SEP SET errechnet
SEP_CL_ALT_max_vector_R = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec(1:punkt_H_ALT,2) ./ Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec(1:punkt_H_ALT,1); 
% % SET_CL_max_vector(length(SET_CL_max_vector)+1) =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(length(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec),2);
SEP_CL_ALT_max_beschl_vector_R = SEP_CL_ALT_max_beschl_vec ./ v_beschl_CL_ALT;
% 
% %t = trapz(SET_CL_max_vector, HE_CL_max_vector_SET)
Steigflug.R_CL1_ALT = trapz(HE_CL_max_vector_SEP_ALT, 1./SEP_CL_ALT_max_vector_R); % Ergebnis macht mit unseren inkorrekten werten sinn denke ich
Steigflug.R_beschl_ALT = trapz(HE_SEP_besch_ALT, 1./SEP_CL_ALT_max_beschl_vector_R);

Steigflug.R_CL_ALT = Steigflug.R_CL1_ALT + Steigflug.R_beschl_ALT;


%% Steigkraftstoffverbrauch CL auf CR

[rho_rho0_H, T_H, a_H, p_p0_H, rho] = Atmos_H(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,3));


Ma_h_CL =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,1) ./ a_H;
Ma_h_CL_beschl = v_beschl_CL ./ a_H(punkt_H_CR,1);

[~,~,b_s_CL_1_s] = SFC_vec(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:punkt_H_CR,3), Ma_h_CL, specs.bypass);
[~,~,b_s_CL_1_s_beschl] = SFC_vec(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(punkt_H_CR,3), Ma_h_CL_beschl, specs.bypass);
b_s_CL = b_s_CL_1_s .* (1/specs.g);
b_s_CL_beschl = b_s_CL_1_s_beschl .* (1/specs.g);
%b_s_CL_ = b_s_CL(:,1);
% b_s_CL_kg_Ns_vec(9,1) = b_s_CL_kg_Ns(length(b_s_CL_kg_Ns),1);

[S_S0_CL] = S_S0_KF_j(0.9, rho_rho0_H, Ma_h_CL, p_p0_H, specs.bypass);
[S_S0_CL_beschl] = S_S0_KF_j(0.9, rho_rho0_H(punkt_H_CR), Ma_h_CL_beschl, p_p0_H(punkt_H_CR), specs.bypass);
S_S0_CL_vec = S_S0_CL(:,1);
S_S0_CL_vec_beschl = S_S0_CL_beschl(:,1);
% S_S0_CL_vec = S_S0_CL(length(S_S0_CL),1);
zu_integrierende_werte_Steigkraftstoff = (b_s_CL .* S_S0_CL_vec .* (S0)) ./ (SEP_CL_max_vector);

zu_integrierende_werte_Steigkraftstoff_beschl = (b_s_CL_beschl .* S_S0_CL_vec_beschl .* (S0)) ./ (SEP_CL_max_beschl_vec);

Steigflug.m_F_CL1 = trapz(HE_CL_max_vector_SEP, zu_integrierende_werte_Steigkraftstoff);
% trapz( zu_integrierende_werte_Steigkraftstoff,HE_CL_max_vector_SEP)

Steigflug.m_F_beschl = trapz(HE_SEP_besch, zu_integrierende_werte_Steigkraftstoff_beschl);

Steigflug.m_F_CL = Steigflug.m_F_CL1 + Steigflug.m_F_beschl;  


%% Steigkraftstoffverbrauch CL auf CR
[rho_rho0_H_ALT, T_H_ALT, a_H_ALT, p_p0_H_ALT, rho_ALT] = Atmos_H(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec(1:punkt_H_ALT,3));

% machzahl fuer SFC
Ma_h_CL_ALT =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec(1:punkt_H_ALT,1) ./ a_H_ALT;
Ma_h_CL_ALT_beschl = v_beschl_CL_ALT ./ a_H_ALT(punkt_H_ALT,1);

[~,~,b_s_CL_ALT_1_s] = SFC_vec(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec(1:punkt_H_ALT,3), Ma_h_CL_ALT, specs.bypass);
[~,~,b_s_CL_ALT_1_s_beschl] = SFC_vec(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_ALT_vec(punkt_H_ALT,3), Ma_h_CL_ALT_beschl, specs.bypass);
b_s_CL_ALT = b_s_CL_ALT_1_s .* (1/specs.g);
b_s_CL_ALT_beschl = b_s_CL_ALT_1_s_beschl .* (1/specs.g);


[S_S0_CL_ALT] = S_S0_KF_j(0.9, rho_rho0_H_ALT, Ma_h_CL_ALT, p_p0_H_ALT, specs.bypass);
[S_S0_CL_ALT_beschl] = S_S0_KF_j(0.9, rho_rho0_H_ALT(punkt_H_ALT), Ma_h_CL_ALT_beschl, p_p0_H_ALT(punkt_H_ALT), specs.bypass);
S_S0_CL_ALT_vec = S_S0_CL_ALT(:,1);
S_S0_CL_ALT_vec_beschl = S_S0_CL_ALT_beschl(:,1);


% Fuelmassenberechnung
zu_integrierende_werte_Steigkraftstoff_ALT = (b_s_CL_ALT .* S_S0_CL_ALT_vec .* (S0)) ./ (SEP_CL_ALT_max_vector);

zu_integrierende_werte_Steigkraftstoff_beschl_ALT = (b_s_CL_ALT_beschl .* S_S0_CL_ALT_vec_beschl .* (S0)) ./ (SEP_CL_ALT_max_beschl_vec);


Steigflug.m_F_CL1_ALT = trapz(HE_CL_max_vector_SEP_ALT, zu_integrierende_werte_Steigkraftstoff_ALT);


Steigflug.m_F_beschl_ALT = trapz(HE_SEP_besch_ALT, zu_integrierende_werte_Steigkraftstoff_beschl_ALT);

Steigflug.m_F_CL_ALT = Steigflug.m_F_CL1_ALT + Steigflug.m_F_beschl_ALT;  




%% Fuel Fraction neu


% Climb auf CR 3 Klassen
FFneu.m3_neu =Ergebnisse_Massen_FE2.M_TO * FF.mf2;
FFneu.mf3 = 1 - (Steigflug.m_F_CL)/(FFneu.m3_neu);


% Reichweite CR 3 Klassen
FFneu.mf4 = (1.05 - (Ergebnisse_Massen_FE2.M_F)/(Ergebnisse_Massen_FE2.M_TO))/(FF.mf2 * FF.mf3 * FF.mf5 *(FF.mf6 * FF.mf7 * FF.mf8* FF.mf9 * FF.mf10 + 0.05));

% sfc_1_s *1/g0 3 Klassen
R_CR = (v_CR / (Ergebnisse_Widerstand_FE2.cW_cA_off_D  * FF.sfc_CR_1PERs )) * log(1/ FFneu.mf4);  

% CLimb alt 3 Klassen
FFneu.m6_neu = Ergebnisse_Massen_FE2.M_TO * FF.mf2 * FF.mf3 * FF.mf5 * FFneu.mf4;
FFneu.mf6 = 1 - (Steigflug.m_F_CL_ALT)/(FFneu.m6_neu);

% Holding 3 Klassen
[sfc_HLD_daNh, sfc_HLD_1PERh, sfc_HLD_1PERs] = SFC(hoehe_HLD, FF.Ma_HLD, specs.bypass);

FFneu.mf9 = 1/(exp(specs.t_HLD * Ergebnisse_Widerstand_FE2.cW_cA_off_D * sfc_HLD_1PERs )); % eventuell muss das g0 noch rausgekürtzt werden







%% NRD


m_OE = Ergebnisse_Massen_FE2.M_OE;
m_TO = Ergebnisse_Massen_FE2.M_TO;
m_fuel = Ergebnisse_Massen_FE2.M_F;
m_ZF = Ergebnisse_Massen_FE2.M_ZF;
% Bestimmung der neuen Massen mf3, mf4, mf9

% m_3 = Ergebnisse_Massen_FE2.M_TO * FF.mf2;
m_alt_u_HLD = (1 - FF.Mff_6_10) * Ergebnisse_Massen_FE2.M_TO;
% mf3_neu = 1 - m_F_CL/m_3;
m_P_A = specs.m_cargo + specs.m_pax;


% Punkt A--------------------------------
m_TO_A = Ergebnisse_Massen_FE2.M_OE + m_P_A + m_alt_u_HLD;
% m_F_HLD_A = (1 - FF.mf9) * m_TO_A;
R_A = 0;
m_RF_A = m_alt_u_HLD;
m_F_A = m_RF_A;
m_TF_A = 0;
A = [R_A; m_P_A; m_F_A; m_TF_A; m_RF_A];

    % All ECO


R_A_ECO = 0;
m_RF_A_ECO = m_alt_u_HLD;
m_F_A_ECO = m_RF_A_ECO;
m_TF_A_ECO = 0;
m_P_A_ECO = specs.m_cargo + specs.m_pax_all_eco;

delta_MP_MPeco = m_P_A_ECO - m_P_A; 
m_fuel_ECO = m_fuel - delta_MP_MPeco;

A_ECO = [R_A_ECO; m_P_A_ECO; m_F_A_ECO; m_TF_A_ECO; m_RF_A_ECO];



% m_TO_A = m_P + m_OE + Fuel.m_F_ALT + Fuel.m_F_HLD;
% m_F_HLD_A = (1-Fuel.Fractions.m_f9) * m_TO_A;
% R_A = 0;
% m_RF_A = Fuel.m_F_ALT + m_F_HLD_A;
% m_F_A = m_RF_A;
% m_TF_A = 0;
% A = [0;m_P_alleco; m_F_A; m_TF_A; m_RF_A];

% B--------------------------------------



[~,~,sfc_1PERs] = SFC_vec(hoehe_CR, specs.Ma_CR, specs.bypass);

m_P_B = m_P_A;
m_TO_B = Ergebnisse_Massen_FE2.M_TO;
m_F_B = m_TO_B -Ergebnisse_Massen_FE2.M_OE - m_P_B;

m_ZF_B = m_P_B + m_OE;
mf4_B = mf4_fun(m_F_B, 0);
R_CR_B = R_NRD_fun(v_CR, sfc_1PERs, mf4_B); %v_CR / (FF.sfc_CR_1PERs * Ergebnisse_Widerstand_FE2.cW_cA_off_D * specs.g) * log(1/mf4_B);
R_B = (Steigflug.R_CL + R_CR_B)/1000;

m_F_C_B = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_B * FF.mf5)) * m_TO_B;
m_RF_OC_B = m_ZF_B * (1/((1-m_F_C_B/m_ZF_B) * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9) -1);

m_RF_B = m_RF_OC_B + m_F_C_B;
m_TF_B = m_F_B - m_RF_B;

B = [R_B; m_P_B; m_F_B; m_TF_B; m_RF_B];

    % B ECO

m_P_B_ECO = m_P_A_ECO;
m_TO_B_ECO = Ergebnisse_Massen_FE2.M_TO;
m_F_B_ECO = m_TO_B_ECO -Ergebnisse_Massen_FE2.M_OE - m_P_B_ECO;


m_ZF_B_ECO = m_P_B_ECO + m_OE; % !!! Annahme, dass M_OE von ECO und 3 KLassen gleich ist
mf4_B_ECO = mf4_fun(m_F_B_ECO, 0);
R_CR_B_ECO = R_NRD_fun(v_CR, sfc_1PERs, mf4_B_ECO);
R_B_ECO = (Steigflug.R_CL + R_CR_B_ECO)/1000;

m_F_C_B_ECO = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_B_ECO * FF.mf5)) * m_TO_B_ECO;
m_RF_OC_B_ECO = m_ZF_B_ECO * (1/((1-m_F_C_B_ECO/m_ZF_B_ECO) * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9) -1);

m_RF_B_ECO = m_RF_OC_B_ECO + m_F_C_B_ECO;
m_TF_B_ECO = m_F_B_ECO - m_RF_B_ECO;


B_ECO = [R_B_ECO; m_P_B_ECO; m_F_B_ECO; m_TF_B_ECO; m_RF_B_ECO];



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

m_F_max = Betankung.Masse_Fuel_Aussen_Theoretisch + Betankung.Masse_Fuel_Innen;
delta_mFmax_MF = m_F_max - m_fuel;


m_F_C = m_F_max;
m_TO_C = m_TO;
m_P_C = m_TO_C - (m_OE + m_F_C);

m_ZF_C = m_OE + m_P_C;
mf4_C = mf4_fun(m_F_C, 0);
R_CR_C = R_NRD_fun(v_CR, sfc_1PERs, mf4_C);%v_CR / (FF.sfc_CR_1PERs * Ergebnisse_Widerstand_FE2.cW_cA_off_D *specs.g) * log(1/mf4_C);
R_C = (Steigflug.R_CL + R_CR_C)/1000;


% m_F_C_C = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_C * FF.mf5)) * m_TO_C;
m_F_C_C = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_B * FF.mf5)) * m_TO_C;
m_RF_OC_C = m_ZF_C * (1/((1-m_F_C_C/m_ZF_C) * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9) -1);
m_RF_C = m_RF_OC_C + m_F_C_C;
m_TF_C = m_F_C - m_RF_C;

C = [R_C; m_P_C; m_F_C; m_TF_C; m_RF_C];


    % C ECO
m_F_max_ECO = Betankung.Masse_Fuel_Aussen_Theoretisch + Betankung.Masse_Fuel_Innen;
delta_mFmax_MF_ECO = m_F_max_ECO - m_fuel_ECO;

m_F_C_ECO = m_F_max_ECO;
m_TO_C_ECO = m_TO;
m_P_C_ECO = m_TO_C_ECO - (m_OE + m_F_C_ECO);

m_ZF_C_ECO = m_OE + m_P_C_ECO;
mf4_C_ECO = mf4_fun(m_F_C_ECO, 0);

R_CR_C_ECO = R_NRD_fun(v_CR, sfc_1PERs, mf4_C_ECO);
R_C_ECO = (Steigflug.R_CL + R_CR_C_ECO)/1000;

m_F_C_C_ECO = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_B_ECO * FF.mf5)) * m_TO_C_ECO;
m_RF_OC_C_ECO = m_ZF_C_ECO * (1/((1-m_F_C_C_ECO/m_ZF_C_ECO) * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9) -1);
m_RF_C_ECO = m_RF_OC_C_ECO + m_F_C_C_ECO;
m_TF_C_ECO = m_F_C_ECO- m_RF_C_ECO;


C_ECO = [R_C_ECO; m_P_C_ECO; m_F_C_ECO; m_TF_C_ECO; m_RF_C_ECO];


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

m_F_D = m_F_max;
m_P_D = 0;
m_TO_D = m_OE + m_F_D;
m_ZF_D = m_OE + m_P_D;
reduktion_M_TO_D = m_TO - m_TO_D; 
mf4_D = mf4_fun(m_F_D, reduktion_M_TO_D);
R_CR_D = R_NRD_fun(v_CR, sfc_1PERs, mf4_D); %v_CR / (FF.sfc_CR_1PERs * Ergebnisse_Widerstand_FE2.cW_cA_off_D *specs.g) * log(1/mf4_D);
R_D = (Steigflug.R_CL + R_CR_D)/1000;

% m_F_C_D = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_D * FF.mf5)) * m_TO_D;
m_F_C_D = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_B * FF.mf5)) * m_TO_D;
m_RF_OC_D = m_ZF_D * (1/((1-m_F_C_D/m_ZF_D) * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9) -1);
m_RF_D = m_RF_OC_D + m_F_C_D;
m_TF_D = m_F_D - m_RF_D;
 
D = [R_D; m_P_D; m_F_D; m_TF_D; m_RF_D];

    % D ECO

m_F_D_ECO = m_F_max_ECO;
m_P_D_ECO = 0;
m_TO_D_ECO = m_OE + m_F_D_ECO + m_P_D_ECO;
m_ZF_D_ECO = m_OE + m_P_D_ECO;
reduktion_M_TO_D_ECO = m_TO - m_TO_D_ECO; 
mf4_D_ECO = mf4_fun(m_F_D_ECO, reduktion_M_TO_D_ECO);
R_CR_D_ECO = R_NRD_fun(v_CR, sfc_1PERs, mf4_D_ECO);
R_D_ECO = (Steigflug.R_CL + R_CR_D_ECO)/1000;

m_F_C_D_ECO = 0.05 *(1 - (FF.mf2 * FF.mf3 * mf4_B_ECO * FF.mf5)) * m_TO_D_ECO;
m_RF_OC_D_ECO = m_ZF_D_ECO * (1/((1-m_F_C_D_ECO/m_ZF_D_ECO) * FF.mf6 * FF.mf7 * FF.mf8 * FF.mf9) -1);
m_RF_D_ECO = m_RF_OC_D_ECO + m_F_C_D_ECO;
m_TF_D_ECO = m_F_D_ECO - m_RF_D_ECO;
 
D_ECO = [R_D_ECO; m_P_D_ECO; m_F_D_ECO; m_TF_D_ECO; m_RF_D_ECO];


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
% ylim([100000 m_TO+10000])
xlim([0 16000])
% X = [Reichweite; Payload; Fuelmasse; Tripfuel; Reservefuel];
p1(1) = plot([A_ECO(1) B_ECO(1) C_ECO(1) D_ECO(1)],[A_ECO(2), B_ECO(2), C_ECO(2), D_ECO(2)]); % 
p1(2) = plot([A_ECO(1) B_ECO(1) C_ECO(1) D_ECO(1)],[A_ECO(3), B_ECO(3), C_ECO(3), D_ECO(3)]); % +m_ZF-delta_mFmax_MF
p1(3) = plot([A_ECO(1) B_ECO(1) C_ECO(1) D_ECO(1)],[A_ECO(4), B_ECO(4), C_ECO(4), D_ECO(4)]); % +m_ZF-delta_mFmax_MF
p1(4) = plot([A_ECO(1) B_ECO(1) C_ECO(1) D_ECO(1)],[A_ECO(5), B_ECO(5), C_ECO(5), D_ECO(5)]); % +m_ZF-delta_mFmax_MF
% p1(5) = plot([0, 20000],[m_TO, m_TO], Color=[0.5 0.5 0.5], LineStyle="--");
% p1(6) = plot([0, 20000],[m_OE, m_OE], Color=[0.5 0.5 0.5], LineStyle="-.");
plot(specs.max_range_basis_km, m_P_A,'rx')

title('Nutzlast-Reichweiten-Diagramm', 'FontSize',25)
legend(p1(1:4),{'Nutzlast', 'Treibstoffmasse', 'Reisekraftstoffmasse', 'Reserve'},... % , 'M_{TO}', 'M_{OE}'
     'Location','eastoutside','FontSize',25);
xlabel('Reichweite in km','FontSize',20)
ylabel('Masse in kg','FontSize',20)




%% Speichern von Ergebnissen
NRD.A = A;
NRD.B = B;
NRD.C = C;
NRD.D = D;
NRD.A_ECO = A_ECO;
NRD.B_ECO = B_ECO;
NRD.C_ECO = C_ECO;
NRD.D_ECO = D_ECO;

NRD.R_A = R_A;
NRD.R_B = R_B;
NRD.R_C = R_C;
NRD.R_D = R_D;
NRD.R_A_ECO = R_A_ECO;
NRD.R_B_ECO = R_B_ECO;
NRD.R_C_ECO = R_C_ECO;
NRD.R_D_ECO = R_D_ECO;


NRD.m_P_A = m_P_A;
NRD.m_P_B = m_P_B;
NRD.m_P_C = m_P_C;
NRD.m_P_D = m_P_D;
NRD.m_P_A_ECO = m_P_A_ECO;
NRD.m_P_B_ECO = m_P_B_ECO;
NRD.m_P_C_ECO = m_P_C_ECO;
NRD.m_P_D_ECO = m_P_D_ECO;


NRD.m_F_A = m_F_A;
NRD.m_F_B = m_F_B;
NRD.m_F_C = m_F_C;
NRD.m_F_D = m_F_D;
NRD.m_F_A_ECO = m_F_A_ECO;
NRD.m_F_B_ECO = m_F_B_ECO;
NRD.m_F_C_ECO = m_F_C_ECO;
NRD.m_F_D_ECO = m_F_D_ECO;

NRD.m_TF_A = m_TF_A;
NRD.m_TF_B = m_TF_B;
NRD.m_TF_C = m_TF_C;
NRD.m_TF_D = m_TF_D;
NRD.m_TF_A_ECO = m_TF_A_ECO;
NRD.m_TF_B_ECO = m_TF_B_ECO;
NRD.m_TF_C_ECO = m_TF_C_ECO;
NRD.m_TF_D_ECO = m_TF_D_ECO;


NRD.m_RF_A = m_RF_A;
NRD.m_RF_B = m_RF_B;
NRD.m_RF_C = m_RF_C;
NRD.m_RF_D = m_RF_D;
NRD.m_RF_A_ECO = m_RF_A_ECO;
NRD.m_RF_B_ECO = m_RF_B_ECO;
NRD.m_RF_C_ECO = m_RF_C_ECO;
NRD.m_RF_D_ECO = m_RF_D_ECO;




Ergebnisse_Flugleistung_2.Steigflug = Steigflug;
Ergebnisse_Flugleistung_2.FFneu = FFneu;
Ergebnisse_Flugleistung_2.NRD = NRD;

save Ergebnisse_Flugleistung_2.mat Ergebnisse_Flugleistung_2 FFneu NRD

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


%% Funktion mf4
function [mf4] = mf4_fun(m_F, M_TO_reduktion)
    load Ergebnisse_Massen_FE2.mat;
    mf4 = (1.05 - (m_F)/(Ergebnisse_Massen_FE2.M_TO - M_TO_reduktion))/(FF.mf2 * FF.mf3 * FF.mf5 *(FF.mf6 * FF.mf7 * FF.mf8* FF.mf9 * FF.mf10 + 0.05));
end

%% Funktion Reichweite NRD

function R = R_NRD_fun(v_CR, sfc_1PERs, mf4)
    load Ergebnisse_Widerstand_FE2.mat;
    load Projekt_specs.mat;

    R = (v_CR / (Ergebnisse_Widerstand_FE2.cW_cA_off_D * sfc_1PERs)) * log(1/ mf4);  % Ergebnisse_Widerstand_FE2.cW_cA_off_D

end



