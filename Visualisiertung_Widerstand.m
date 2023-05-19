% Visualisierung Widerstand

clc
clear all
close all

%% Ausfuerung des Programms
load Ergebnisse_Endwerte_Iteration_V1.mat
Endwerte_Iteration = Berechnungen_PS10_Widerstand;

%% Laden der Dateien

load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Basis_stat_m.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Start_Landeanforderungen.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Widerstand.mat;


%% Visualisierung


figure(1)
hold on
grid on
title('Lilienthalpolare')
plot(Widerstand.C_w_clean_all,Widerstand.y_CR,'k')
plot(Widerstand.C_w_LR_all,Widerstand.y_CR,'b')
plot(Widerstand.C_w_HS_all,Widerstand.y_CR,'--k')
legend('Cruise Clean Configuration','Cruise Long Range','Cruise High Speed','Location','southeast')

xlabel('Widerstandsbeiwert C_w')
ylabel('Auftriebsbeiwert C_A')
hold off



figure(2)
hold on
grid on
title('Lilienthalpolare')
plot(Widerstand.C_w_clean_all,Widerstand.y_CR,'k')
plot(Widerstand.C_w_TO_clean_all,Widerstand.y_to)
plot(Widerstand.C_w_TO_all,Widerstand.y_to)

% Landung

plot(Widerstand.C_w_LDG_clean_all,Widerstand.y)
plot(Widerstand.C_w_LDG_all,Widerstand.y)

legend('Cruise Clean Configuration','TO Clean Configuration','TO Gear Down','LD Clean Configuration','LD Gear Down', 'Location','southeast')

xlabel('Widerstandsbeiwert C_w')
ylabel('Auftriebsbeiwert C_A')
ylim([0 2.7])
hold off

figure(3)
hold on
grid on
title('Reziproke Gleitzahl')
plot(Widerstand.y_CR,GZ.CA_CW_Clean, 'b');
plot(Widerstand.y_CR,GZ.CA_CW_LR,'b--');
plot(Widerstand.y_CR,GZ.CA_CW_HS,'b-.');
plot(Widerstand.y_to,GZ.CA_CW_TO_Clean,'r');
plot(Widerstand.y_to,GZ.CA_CW_TO,'r--');
plot(Widerstand.y,GZ.CA_CW_LDG_Clean,'m');
plot(Widerstand.y,GZ.CA_CW_LDG,'m--')
xlabel("Auftriebsbeiwert C_A")
ylabel("Reziproke Gleitzahl E")

xlim([0 2.7])


plot(Ergebnisse_stat_Flaechenbelastung.C_A_CR, schub_CR.Eta, 'oblack','MarkerSize', 8)
plot(startschub.c_A_max_thrust_match, startschub.Eta_To_inv(3,1),'xblack','MarkerSize', 8)
plot(landeanvorderung.c_A_max_LDG,(1/landeanvorderung.Eta_LDG),'*black','MarkerSize', 8)

plot(Widerstand.y_CR(1,(round(Ergebnisse_stat_Flaechenbelastung.C_A_CR * 10^3))),...
    GZ.CA_CW_LR(1,(round(Ergebnisse_stat_Flaechenbelastung.C_A_CR * 10^3))),'or','MarkerSize', 8) % 523
plot(Widerstand.y_to(1,(round(startschub.c_A_max_thrust_match * 10^3))),...
    GZ.CA_CW_TO(1,(round(startschub.c_A_max_thrust_match * 10^3))),'xr','MarkerSize', 8) % 1801
plot(Widerstand.y(1,(round(landeanvorderung.c_A_max_LDG * 10^3))),...
    GZ.CA_CW_LDG(1,(round(landeanvorderung.c_A_max_LDG * 10^3))),'*red','MarkerSize', 8) % 2401

legend("Cruise Clean","Cruise LR","Cruise HS","TO Clean","TO Gear Down","LDG clean","LDG Gear Down",...
    'Reale Gleitzahl CR','Reale Gleitzahl TO','Reale Gleitzahl LDG','Ideale Gleitzahl CR','Ideale Gleitzahl TO','Ideale Gleitzahl LDG')

hold off