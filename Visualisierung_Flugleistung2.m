%% Visualisierung PS8 Flugleistung 2 NRD 
clc
clear all
close all

Berechnung_FE2_PS8_Flugleistung2_NRD;

load Projekt_specs.mat;
load Ergebnisse_Flugleistung_2.mat;


figure(1) 
hold on 
grid on 
% ylim([100000 m_TO+10000])
xlim([0 16000])
% X = [Reichweite; Payload; Fuelmasse; Tripfuel; Reservefuel];
p1(1) = plot([NRD.A_ECO(1) NRD.B_ECO(1) NRD.C_ECO(1) NRD.D_ECO(1)],...
    [NRD.A_ECO(2), NRD.B_ECO(2), NRD.C_ECO(2), NRD.D_ECO(2)]); % 
p1(2) = plot([NRD.A_ECO(1) NRD.B_ECO(1) NRD.C_ECO(1) NRD.D_ECO(1)],...
    [NRD.A_ECO(3), NRD.B_ECO(3), NRD.C_ECO(3), NRD.D_ECO(3)]); % +m_ZF-delta_mFmax_MF
p1(3) = plot([NRD.A_ECO(1) NRD.B_ECO(1) NRD.C_ECO(1) NRD.D_ECO(1)],...
    [NRD.A_ECO(4), NRD.B_ECO(4), NRD.C_ECO(4), NRD.D_ECO(4)]); % +m_ZF-delta_mFmax_MF
p1(4) = plot([NRD.A_ECO(1) NRD.B_ECO(1) NRD.C_ECO(1) NRD.D_ECO(1)],...
    [NRD.A_ECO(5), NRD.B_ECO(5), NRD.C_ECO(5), NRD.D_ECO(5)]); % +m_ZF-delta_mFmax_MF
% p1(5) = plot([0, 20000],[m_TO, m_TO], Color=[0.5 0.5 0.5], LineStyle="--");
% p1(6) = plot([0, 20000],[m_OE, m_OE], Color=[0.5 0.5 0.5], LineStyle="-.");
p1(5) = plot(specs.max_range_basis_km, NRD.m_P_A,'rx');

title('Nutzlast-Reichweiten-Diagramm für All-Eco', 'FontSize',25)
legend(p1(1:5),{'Nutzlast', 'Treibstoffmasse', 'Reisekraftstoffmasse', 'Reservekraftstoff', 'DP'},... % , 'M_{TO}', 'M_{OE}'
     'Location','eastoutside','FontSize',25);
xlabel('Reichweite in km','FontSize',20)
ylabel('Masse in kg','FontSize',20)