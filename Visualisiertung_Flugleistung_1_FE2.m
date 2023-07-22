% FE2 PS7 Visualisierung Flugleistung 1 


clc
clear all
close all


load Projekt_specs.mat
load Ergebnisse_FLugleistung_1.mat
load Ergebnisse_Widerstand_FE2.mat

% Anzahl der Plots
numPlots = length(Ergebnisse_Flugleistung_1.hoehe_m);

% Farbverlauf definieren
colorStart = [1, 0.2, 0];   % Startfarbe (RGB)
colorEnd = [0, 0, 1];     % Endfarbe (RGB)

% Farbwerte für jeden Plot berechnen
colors_1 = zeros(numPlots, 3);
colors_2 = zeros(numPlots, 3);
for i = 1:numPlots
    colors_1(i, :) = colorStart + (i-1) * (colorEnd - colorStart) / (numPlots-1);
    colors_2(i, :) = colorStart + (i-1) * (colorEnd - colorStart) / (numPlots-1);
end

% Linienarten definieren
lineStyles = {'--', '-', '-.', '-', '--', '-', '-.', '-', '--', '-', '-.', '-' }; % Gestrichelt, Gepunktet


%% Horizontalflugdiagramm

figure(1) % Horiontalflugdiagramm
hold on
grid on
title('Horizontalflugdiagramm G_{ICA}', 'FontSize',25);
ylabel('S/G', 'FontSize',20);
xlabel('v_{EAS}', 'FontSize',20);
xlim([0 300]);
ylim([0 1]);

% p1_1(n_datensatz) = plot(v_EAS_j, S_G_vorh_j(n_datensatz, :), 'Color',colors_1(n_datensatz,:), 'LineStyle','--');
% p1_2(n_datensatz) = plot(v_EAS_j, S_G_erf_j(n_datensatz, :),  'Color',colors_2(n_datensatz,:), 'LineStyle','-');
% Plots erstellen
for n_plot = 1:numPlots
    autoplot_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_EAS_j, Ergebnisse_Flugleistung_1.S_G_vorh_j(n_plot,:), 'LineStyle', lineStyles{1,n_plot}); % , 'Color', colors(n_plot, :), 'LineStyle', lineStyles{1,n_plot});
    autoplot_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_EAS_j, Ergebnisse_Flugleistung_1.S_G_erf_j(n_plot,:), 'LineStyle', lineStyles{1,n_plot}); % , 'Color', colors(n_plot, :), 'LineStyle', lineStyles{1,n_plot});
    Legend{2*n_plot-1} =  '(S/G)_{erf} Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot)) + 'm' ;
    Legend{2*n_plot} =  '(S/G)_{vorh} Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot)) + 'm';
    legend(Legend(1:n_plot*2),Location='eastoutside', FontSize=10);

end


% plot(v_EAS_j(n_datensatz, :), eps_inkomp_j(n_datensatz, :), '--b');
hold off


%% SET

figure(2) % SET
hold on 
grid on
title('Spezifischer Schubüberschuss G_{ICA}', 'FontSize',25)
ylabel('SET','FontSize',20);
xlabel('v_{TAS}','FontSize',20);
xlim([0 350]);
ylim([0 0.7]);

% p2_1(n_datensatz) = plot(v_TAS_j(n_datensatz, :), SET(n_datensatz, :), 'Color',colors_2(n_datensatz,:), 'LineStyle','-');
% p2_2(n_datensatz) = plot(v_TAS_j(n_datensatz, SET_x(n_datensatz,1)), SET_y(n_datensatz,1), '*k');
% Legend2{n_datensatz*2-1} = 'SET in [rad] Hoehe ' + sprintf("%d",hoehe_m(n_datensatz)) + 'm';
% Legend2{n_datensatz*2} = '';
% legend(Legend2(1:n_datensatz*2), Location='eastoutside', FontSize=16)


for n_plot = 1:numPlots
    autoplot2_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot,:), Ergebnisse_Flugleistung_1.SET(n_plot,:), 'LineStyle', lineStyles{1,n_plot}); % , 'Color', colors(n_plot, :), 'LineStyle', lineStyles{1,n_plot});
    autoplot2_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot, Ergebnisse_Flugleistung_1.Hochpunkte.SET_x(n_plot,1)), Ergebnisse_Flugleistung_1.Hochpunkte.SET_y(n_plot,:), '*b'); % , 'Color', colors(n_plot, :), 'LineStyle', lineStyles{1,n_plot});
    
    Legend2{2*n_plot-1} =  'SET in rad, Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot)) + 'm' ;
    Legend2{2*n_plot} =  '' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot)) + 'm';
    legend(Legend2(1:n_plot*2),Location = 'eastoutside', FontSize=10);

end

hold off

%% SEP

figure(3)
hold on 
grid on
title('Spezifischer Überschussleistung G_{ICA}', 'FontSize',25)
ylabel('SEP','FontSize',20);
xlabel('v_{TAS}','FontSize',20);
xlim([0 350]);
ylim([0 100]);

% p3_1(n_datensatz) = plot(v_TAS_j(n_datensatz, :), SEP(n_datensatz, :), 'Color',colors_2(n_datensatz,:), 'LineStyle','-');
% p3_2(n_datensatz) = plot(v_TAS_j(n_datensatz, SEP_x(n_datensatz,1)), SEP_y(n_datensatz,1), '*b');
for n_plot = 1:numPlots
    autoplot3_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot, :), Ergebnisse_Flugleistung_1.SEP(n_plot, :), 'Color',colors_2(n_plot,:), 'LineStyle','-');
    autoplot3_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot, Ergebnisse_Flugleistung_1.Hochpunkte.SEP_x(n_plot,1)), Ergebnisse_Flugleistung_1.Hochpunkte.SEP_y(n_plot,1), '*b');
    Legend3{n_plot*2-1} = 'SEP in m/s, Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot)) + 'm' ;
    Legend3{n_plot*2} ='';
    legend(Legend3(1:n_plot*2), Location = 'eastoutside', FontSize=16)
end

hold off

%% SR

figure(4) % SR

hold on 
grid on
title('Spezifische Reichweite G_{ICA}','FontSize',25)
ylabel('SR in m/kg','FontSize',20);
xlabel('v_{TAS} in m/s','FontSize',20);
xlim([0 300]);
ylim([0 20]);

for n_plot = 1:numPlots
    autoplot4_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot, :), Ergebnisse_Flugleistung_1.SR(n_plot, :), 'Color',colors_2(n_plot,:), 'LineStyle','-');
    autoplot4_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot, Ergebnisse_Flugleistung_1.Hochpunkte.SR_x(n_plot,1)), Ergebnisse_Flugleistung_1.Hochpunkte.SR_y(n_plot,1), '*k');
    Legend4{n_plot*2-1} = 'SEP in m/kg Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot)) + 'm' ;
    Legend4{n_plot*2} = '';
    legend(Legend4(1:n_plot*2), Location='eastoutside', FontSize=16)

end

hold off

%% Spezifische Flugdauer
figure(5) 

hold on 
grid on 
title('Spezifische Flugdauer G_{ICA}','FontSize',25)
ylabel('SE in s/kg','FontSize',20);
xlabel('v_{TAS} in m/s','FontSize',20);
xlim([0 350]);

for n_plot = 1:numPlots
autoplot5_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot, :) , Ergebnisse_Flugleistung_1.SE(n_plot,:), 'Color',colors_2(n_plot,:), 'LineStyle','-' );
autoplot5_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot, Ergebnisse_Flugleistung_1.Hochpunkte.SE_x(n_plot,1)), Ergebnisse_Flugleistung_1.Hochpunkte.SE_y(n_plot,1), '*k');
Legend5{n_plot*2-1} = ' in s/kg, Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot)) + 'm' ;
Legend5{n_plot*2} = '';
legend(Legend5(1:n_plot*2), Location='eastoutside', FontSize=16)

end
hold off

%% Optimalgeschwindigkeiten 

figure(6)
hold on 
grid on 

title('Optimalgeschwindigkeiten bei G_{LDG}','FontSize',25)

ylabel('H in m','FontSize',20);
xlabel('v_{TAS} in m/s','FontSize',20);
% xlim([0 350]);


plot(Ergebnisse_Flugleistung_1.TAS_SR_H_vec(1:10,1), Ergebnisse_Flugleistung_1.TAS_SR_H_vec(1:10,3),'k');
plot(Ergebnisse_Flugleistung_1.TAS_SEP_H_vec(1:10,1), Ergebnisse_Flugleistung_1.TAS_SEP_H_vec(1:10,3),'b');
plot(Ergebnisse_Flugleistung_1.TAS_SET_H_vec(1:10,1), Ergebnisse_Flugleistung_1.TAS_SET_H_vec(1:10,3),'r');

legend('SR', 'SEP', 'SET', 'FontSize',15)

hold off