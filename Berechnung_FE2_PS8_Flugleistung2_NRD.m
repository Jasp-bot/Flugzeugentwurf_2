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

TO_Masse =  Ergebnisse_Massen_FE2.M_TO;
G_TO = TO_Masse * specs.g;
k_CR = 0.98;
S0 = k_CR * G_TO * (1/schub_CR.Eta) / (schub_CR.S_S0_CR * schub_CR.S_S0_E);

%--------------------------------------------------------------------------
%% Rechnungen
%--------------------------------------------------------------------------

%% Energiehoehe HE fuer Steigzeit

HE_SEP_max_CL = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(:,3) + (Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(:,1).^2)/(2* specs.g);
HE_CL_max_vector_SEP = HE_SEP_max_CL(1:8);
HE_CL_max_vector_SEP(length(HE_CL_max_vector_SEP)+1) = HE_SEP_max_CL(length(HE_SEP_max_CL),1);

%plot(HE_SEP_max_CL, 1./Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(:,2),'*k')

SEP_CL_max_vector = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:8,2);
SEP_CL_max_vector(length(SEP_CL_max_vector)+1) =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(length(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec),2);


%t = trapz(SEP_CL_max_vector, HE_CL_max_vector_SEP)
t = trapz(HE_CL_max_vector_SEP, 1./SEP_CL_max_vector) % Ergebnis macht mit unseren inkorrekten werten sinn denke ich



%% Steigstrecke

HE_SET_max_CL = Ergebnisse_Flugleistung_1.TAS_SET_H_CL_vec(:,3) + (Ergebnisse_Flugleistung_1.TAS_SET_H_CL_vec(:,1).^2)/(2* specs.g);
HE_CL_max_vector_SET = HE_SET_max_CL(1:8);
HE_CL_max_vector_SET(length(HE_CL_max_vector_SET)+1) = HE_SEP_max_CL(length(HE_SEP_max_CL),1);

%plot(HE_SEP_max_CL, 1./Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(:,2),'*k')

SET_CL_max_vector = Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(1:8,2);
SET_CL_max_vector(length(SET_CL_max_vector)+1) =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(length(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec),2);


%t = trapz(SET_CL_max_vector, HE_CL_max_vector_SET)
R = trapz(HE_CL_max_vector_SET, 1./SET_CL_max_vector) % Ergebnis macht mit unseren inkorrekten werten sinn denke ich

%% Steigkraftstoffverbrauch

[rho_rho0_H, T_H, a_H, p_p0_H, rho] = Atmos_H(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(:,3));


Ma_h_CL =  Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(:,1) ./ a_H;
 
[b_s_CL_kg_daNh,~,~] = SFC_vec(Ergebnisse_Flugleistung_1.TAS_SEP_H_CL_vec(:,3), Ma_h_CL, specs.bypass);
b_s_CL_kg_Ns = b_s_CL_kg_daNh .* (1/10) .* (1/3600);
b_s_CL_kg_Ns_vec = b_s_CL_kg_Ns(1:8,1);
b_s_CL_kg_Ns_vec(9,1) = b_s_CL_kg_Ns(length(b_s_CL_kg_Ns),1);

[S_S0_CL] = S_S0_KF_j(0.9, rho_rho0_H, Ma_h_CL, p_p0_H, specs.bypass);
S_S0_CL_vec = S_S0_CL(1:8,1);
S_S0_CL_vec = S_S0_CL(length(S_S0_CL),1);
zu_integrierende_werte_Steigkraftstoff = (b_s_CL_kg_Ns_vec .* S_S0_CL_vec .* (S0/G_TO) .* G_TO) ./ (SEP_CL_max_vector);

m_F_CL = trapz(HE_CL_max_vector_SEP, zu_integrierende_werte_Steigkraftstoff)
% trapz( zu_integrierende_werte_Steigkraftstoff,HE_CL_max_vector_SEP)
%% NRD



















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


