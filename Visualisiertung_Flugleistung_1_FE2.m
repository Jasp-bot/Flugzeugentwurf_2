% FE2 PS7 Visualisierung Flugleistung 1 


%% Kommentar
% Es muessen noch die Axen angepasst werden (Schriftgroesse) !!!!!
% v_S_1g und v_min schneiden sich noch nocht im Flugbereichsdiagramm


Berechnung_FE2_PS7_Flugleistung1

clc
clear all
close all


%% Plots verwalten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wert zwischen 1 und 0 wählen 
% 1 = wird displayed || 0 = wird nicht displayed

Aufgabe_1 = 1;  % Horizontalflugdiagramme fuer CL, CR, DEC
Aufgabe_2 = 0;  % SET, SEP, SR, fuer CR
Aufgabe_3 = 0;  % Optimalgeschwindigheiten CR
Aufgabe_4 = 0;  % Spez Flugdauer SE fuer DEC
Aufgabe_5 = 1;  % Flugbereichsdiagramme fuer CR, DEC

Disp_vec = [Aufgabe_1, Aufgabe_2, Aufgabe_3, Aufgabe_4, Aufgabe_5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vorarbeiten fuer Plots


load Projekt_specs.mat
load Ergebnisse_ISA_DATA.mat
load Ergebnisse_FLugleistung_1.mat
load Ergebnisse_Widerstand_FE2.mat
load Ergebnisse_Massen_FE2.mat



% Anzahl der Plots
numPlots = length(Ergebnisse_Flugleistung_1.hoehe_m);

w = 10;     % Variabele zur reduktion der geplotteten graphen es empiehlt sich abhaengig von der anzahn der stuetzstellen eine zahl zu waehlen

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

%% Aufgabe 1
if Disp_vec(1,1) == 1
    %% Horizontalflugdiagramm CR
    
    figure(1) % Horiontalflugdiagramm
    hold on
    grid on
    title('Horizontalflugdiagramm G_{ICA}', 'FontSize',25);
    ylabel('S/G', 'FontSize',20);
    xlabel('v_{EAS}', 'FontSize',20);
    xlim([0 300]);
    ylim([0 1]);
    
    
    for n_plot = 1:numPlots/w
    
        autoplot_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_EAS_j, Ergebnisse_Flugleistung_1.S_G_vorh_j(n_plot*w,:), 'LineStyle', lineStyles{1,1}, 'Color', colors_1(n_plot*w, :)); %'LineStyle', lineStyles{1,n_plot});
        autoplot_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_EAS_j, Ergebnisse_Flugleistung_1.S_G_erf_j(n_plot*w,:), 'LineStyle', lineStyles{1,2}, 'Color', colors_1(n_plot*w, :)); %, 'LineStyle', lineStyles{1,n_plot});
        Legend{2*n_plot-1} =  '(S/G)_{erf} Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm' ;
        Legend{2*n_plot} =  '(S/G)_{vorh} Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm';
        legend(Legend(1:n_plot*2),Location='eastoutside', FontSize=10);
    
    end
    
    
    % plot(v_EAS_j(n_datensatz, :), eps_inkomp_j(n_datensatz, :), '--b');
    hold off
    
    %% Horizontalflugdiagramm CL
    figure(8)
    
    hold on
    grid on
    title('Horizontalflugdiagramm G_{TO}', 'FontSize',25);
    ylabel('S/G', 'FontSize',20);
    xlabel('v_{EAS}', 'FontSize',20);
    xlim([0 300]);
    ylim([0 1]);
    
    
    for n_plot = 1:numPlots/w
    
        autoplot_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_EAS_j_CL, Ergebnisse_Flugleistung_1.S_G_vorh_j_CL(n_plot*w,:), 'LineStyle', lineStyles{1,1}, 'Color', colors_1(n_plot*w, :)); %'LineStyle', lineStyles{1,n_plot});
        autoplot_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_EAS_j_CL, Ergebnisse_Flugleistung_1.S_G_erf_j_CL(n_plot*w,:), 'LineStyle', lineStyles{1,2}, 'Color', colors_1(n_plot*w, :)); %, 'LineStyle', lineStyles{1,n_plot});
        Legend{2*n_plot-1} =  '(S/G)_{erf} Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm' ;
        Legend{2*n_plot} =  '(S/G)_{vorh} Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm';
        legend(Legend(1:n_plot*2),Location='eastoutside', FontSize=10);
    
    end
    
    hold off
    
    %% Horizontalflugdiagramm DEC
    figure(9)
    
    hold on
    grid on
    title('Horizontalflugdiagramm G_{LDG}', 'FontSize',25);
    ylabel('S/G', 'FontSize',20);
    xlabel('v_{EAS}', 'FontSize',20);
    xlim([0 300]);
    ylim([0 1]);
    
    
    for n_plot = 1:numPlots/w
    
        autoplot_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_EAS_j_DEC, Ergebnisse_Flugleistung_1.S_G_vorh_j_DEC(n_plot*w,:), 'LineStyle', lineStyles{1,1}, 'Color', colors_1(n_plot*w, :)); %'LineStyle', lineStyles{1,n_plot});
        autoplot_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_EAS_j_DEC, Ergebnisse_Flugleistung_1.S_G_erf_j_DEC(n_plot*w,:), 'LineStyle', lineStyles{1,2}, 'Color', colors_1(n_plot*w, :)); %, 'LineStyle', lineStyles{1,n_plot});
        Legend{2*n_plot-1} =  '(S/G)_{erf} Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm' ;
        Legend{2*n_plot} =  '(S/G)_{vorh} Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm';
        legend(Legend(1:n_plot*2),Location='eastoutside', FontSize=10);
    
    end
    
    
    hold off

end
%% Aufgabe 2
if Disp_vec(1,2) == 1
    %% SET
    
    figure(2) % SET
    hold on 
    grid on
    title('Spezifischer Schubüberschuss G_{ICA}', 'FontSize',25)
    ylabel('SET','FontSize',20);
    xlabel('v_{TAS}','FontSize',20);
    xlim([0 350]);
    ylim([0 0.1]);
    
    
    for n_plot = 1:numPlots/w
        autoplot2_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot*w,:), Ergebnisse_Flugleistung_1.SET(n_plot*w,:), 'LineStyle', lineStyles{1,2}, 'Color', colors_1(n_plot*w, :)); 
        autoplot2_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot*w, Ergebnisse_Flugleistung_1.Hochpunkte.SET_x(n_plot*w,1)), Ergebnisse_Flugleistung_1.Hochpunkte.SET_y(n_plot*w,:), '*b');     
        Legend2{2*n_plot-1} =  'SET in rad, Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm' ;
        Legend2{2*n_plot} =  '' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm';
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
    ylim([0 15]);
    
    
    for n_plot = 1:numPlots/w
        autoplot3_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot*w, :), Ergebnisse_Flugleistung_1.SEP(n_plot*w, :), 'Color',colors_2(n_plot*w,:), 'LineStyle','-');
        autoplot3_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot*w, Ergebnisse_Flugleistung_1.Hochpunkte.SEP_x(n_plot*w,1)), Ergebnisse_Flugleistung_1.Hochpunkte.SEP_y(n_plot*w,1), '*b');
        Legend3{n_plot*2-1} = 'SEP in m/s, Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm' ;
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
    ylim([0 250]);
    
    for n_plot = 1:numPlots/w
        autoplot4_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot*w, :), Ergebnisse_Flugleistung_1.SR(n_plot*w, :), 'Color',colors_2(n_plot*w,:), 'LineStyle','-');
        autoplot4_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j(n_plot*w, Ergebnisse_Flugleistung_1.Hochpunkte.SR_x(n_plot*w,1)), Ergebnisse_Flugleistung_1.Hochpunkte.SR_y(n_plot*w,1), '*k');
        Legend4{n_plot*2-1} = 'SR in m/kg Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm' ;
        Legend4{n_plot*2} = '';
        legend(Legend4(1:n_plot*2), Location='eastoutside', FontSize=16)
    
    end
    
    hold off

end
%% Aufgabe 3
if Disp_vec(1,3) == 1
    %% Optimalgeschwindigkeiten 
    
    figure(6)
    hold on 
    grid on 
    
    title('Optimalgeschwindigkeiten bei G_{LDG}','FontSize',25)
    
    ylabel('H in m','FontSize',20);
    xlabel('v_{TAS} in m/s','FontSize',20);
    xlim([0 350]);
    
    
    plot(Ergebnisse_Flugleistung_1.TAS_SR_H_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SR_H_vec(:,3),'k');
    plot(Ergebnisse_Flugleistung_1.TAS_SEP_H_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SEP_H_vec(:,3),'b');
    plot(Ergebnisse_Flugleistung_1.TAS_SET_H_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SET_H_vec(:,3),'r');
    
    legend('SR', 'SEP', 'SET', 'FontSize',15)
    
    hold off

end
%% Aufgabe 4
if Disp_vec(1,4) == 1
    %% Spezifische Flugdauer
    
    figure(5) 
    hold on 
    grid on 
    title('Spezifische Flugdauer G_{ICA}','FontSize',25)
    ylabel('SE in s/kg','FontSize',20);
    xlabel('v_{TAS} in m/s','FontSize',20);
    xlim([0 350]);
    ylim([0 1.5]);
    
    for n_plot = 1:numPlots/w
    autoplot5_1(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j_DEC(n_plot*w, :) , Ergebnisse_Flugleistung_1.SE_DEC(n_plot*w,:), 'Color',colors_2(n_plot*w,:), 'LineStyle','-' );
    autoplot5_2(n_plot,1) = plot(Ergebnisse_Flugleistung_1.v_TAS_j_DEC(n_plot*w, Ergebnisse_Flugleistung_1.Hochpunkte.SE_DEC_x(n_plot*w,1)), Ergebnisse_Flugleistung_1.Hochpunkte.SE_DEC_y(n_plot*w,1), '*k');
    Legend5{n_plot*2-1} = ' in s/kg, Hoehe ' + sprintf("%d",Ergebnisse_Flugleistung_1.hoehe_m(n_plot*w)) + 'm' ;
    Legend5{n_plot*2} = '';
    legend(Legend5(1:n_plot*2), Location='eastoutside', FontSize=16)
    
    end
    hold off

end
%% Aufgabe 5
if Disp_vec(1,5) == 1
    
    %% Flugbereichsdiagramm
    
    
    
    hoehe_CR = round(unitsratio('m','ft')*(specs.flight_level*10^2));
    
    figure(7)
    hold on
    grid on
    
    title('Flugbereichsdiagramm Reiseflug', FontSize=25)
    xlabel('v_{TAS} in m/s', FontSize=20);
    ylabel('H in m', FontSize=20);
    xlim([0 450]);
    ylim([Ergebnisse_Flugleistung_1.hoehe_m(1,1), 12500]);

    
    pl(1) = plot(Ergebnisse_Flugleistung_1.v_s_1g, Ergebnisse_Flugleistung_1.hoehe_m,...
        '--r', 'DisplayName', 'v_s_{1g}');
    pl(2) = plot(Ergebnisse_Flugleistung_1.v_s_min, Ergebnisse_Flugleistung_1.hoehe_m,...
        'r', 'DisplayName', 'v_s_{min}');
    pl(3) = plot(Ergebnisse_Flugleistung_1.v_BO, Ergebnisse_Flugleistung_1.hoehe_m,...
        'm', 'DisplayName', 'Ma_{BO}');
    
    pl(4) = plot(Ergebnisse_Flugleistung_1.TAS_SR_H_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SR_H_vec(:,3),...
        '.-k', 'DisplayName', 'SR_{max}');
    pl(5) = plot(Ergebnisse_Flugleistung_1.TAS_SEP_H_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SEP_H_vec(:,3),...
        '.-b', 'DisplayName', 'SEP_{max}');
    pl(6) = plot(Ergebnisse_Flugleistung_1.TAS_SET_H_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SET_H_vec(:,3),...
        '.-r', 'DisplayName', 'SET_{max}');
    pl(7) = plot(Ergebnisse_Flugleistung_1.TAS_SE_H_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SE_H_vec(:,3),...
        '.-g', 'DisplayName', 'SE_{max}');
    
    pl(8) = plot(Ergebnisse_Flugleistung_1.v_MO, Ergebnisse_Flugleistung_1.hoehe_m,...
        '.-m', 'DisplayName', 'v_{MO}');
    
    % 
    pl(9) = plot(Ergebnisse_Flugleistung_1.v_min_HFD, Ergebnisse_Flugleistung_1.hoehe_m.',...
        'g', 'DisplayName', 'v_{min}');
    pl(10) = plot(Ergebnisse_Flugleistung_1.v_max_HFD, Ergebnisse_Flugleistung_1.hoehe_m,...
        '--g', 'DisplayName', 'v_{max}');
    
    pl(11) = plot((specs.Ma_CR * ISA.a(hoehe_CR)), hoehe_CR,...
        '*r', 'DisplayName', 'DP'); % Design point
    
    pl(12) = plot([0 500],...
        [Ergebnisse_Flugleistung_1.H_Dienstgipfel_CR, Ergebnisse_Flugleistung_1.H_Dienstgipfel_CR],...
        'Color', [0.3 0.3 0.3],'LineStyle', '--', 'DisplayName', 'H_{ops_{max}}');
    pl(13) = plot([0 500],...
        [Ergebnisse_Flugleistung_1.H_Kabienendruck_CR, Ergebnisse_Flugleistung_1.H_Kabienendruck_CR],...
        'k', 'DisplayName', 'H_{p,max}');
    
%     legend('v_s_{1g}', 'v_s_{min}', 'Ma_{BO}', 'SR_{max}', 'SEP_{max}', 'SET_{max}', 'SE_{max}', 'v_{MO}', 'v_{min}', 'v_{max}', 'DP', 'H_{ops_{max}}', 'H_{p,max}', 'FontSize',15 ,Location='eastoutside');
    legend('FontSize',15 ,Location='eastoutside')
    
    hold off
    
    
    figure(10) % Flugbereichsdiagramm LDG
    
    
    hold on
    grid on
    
    title('Flugbereichsdiagramm im Landefall', FontSize=25)
    xlabel('v_{TAS} in m/s', FontSize=20);
    ylabel('H in m', FontSize=20);
    xlim([0 450]);
    ylim([Ergebnisse_Flugleistung_1.hoehe_m(1,1), 12500]);

    pl_DEC(1) = plot(Ergebnisse_Flugleistung_1.v_s_1g_DEC, Ergebnisse_Flugleistung_1.hoehe_m,...
        '--r', 'DisplayName', 'v_s_{1g}');
    pl_DEC(2) = plot(Ergebnisse_Flugleistung_1.v_s_min_DEC, Ergebnisse_Flugleistung_1.hoehe_m,...
        'r', 'DisplayName', 'v_s_{min}');
    pl_DEC(3) = plot(Ergebnisse_Flugleistung_1.v_BO_DEC, Ergebnisse_Flugleistung_1.hoehe_m,...
        'm', 'DisplayName', 'Ma_{BO}');
    
    pl_DEC(4) = plot(Ergebnisse_Flugleistung_1.TAS_SR_H_DEC_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SR_H_DEC_vec(:,3),...
        '.-k', 'DisplayName', 'SR_{max}');
    pl_DEC(5) = plot(Ergebnisse_Flugleistung_1.TAS_SEP_H_DEC_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SEP_H_DEC_vec(:,3),...
        '.-b', 'DisplayName', 'SEP_{max}');
    pl_DEC(6) = plot(Ergebnisse_Flugleistung_1.TAS_SET_H_DEC_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SET_H_DEC_vec(:,3),...
        '.-r', 'DisplayName', 'SET_{max}');
    pl_DEC(7) = plot(Ergebnisse_Flugleistung_1.TAS_SE_H_DEC_vec(:,1), Ergebnisse_Flugleistung_1.TAS_SE_H_DEC_vec(:,3),...
        '.-g', 'DisplayName', 'SE_{max}');
    
    pl_DEC(8) = plot(Ergebnisse_Flugleistung_1.v_MO_DEC, Ergebnisse_Flugleistung_1.hoehe_m,...
        '.-m', 'DisplayName', 'v_{MO}');
    
    
    pl_DEC(9) = plot(Ergebnisse_Flugleistung_1.v_min_HFD_DEC, Ergebnisse_Flugleistung_1.hoehe_m.',...
        'g', 'DisplayName', 'v_{min}');
    pl_DEC(10) = plot(Ergebnisse_Flugleistung_1.v_max_HFD_DEC, Ergebnisse_Flugleistung_1.hoehe_m,...
        '--g', 'DisplayName', 'v_{max}');
    
    pl_DEC(11) = plot((specs.Ma_CR * ISA.a(hoehe_CR)), hoehe_CR,...
        '*r', 'DisplayName', 'DP'); % Design point
    
    pl_DEC(12) = plot([0 500],...
        [Ergebnisse_Flugleistung_1.H_Dienstgipfel_DEC, Ergebnisse_Flugleistung_1.H_Dienstgipfel_DEC],...
        'Color', [0.3 0.3 0.3],'LineStyle', '--', 'DisplayName', 'H_{ops_{max}}');
    pl_DEC(13) = plot([0 500],...
        [Ergebnisse_Flugleistung_1.H_Kabienendruck_DEC, Ergebnisse_Flugleistung_1.H_Kabienendruck_DEC],...
        'k', 'DisplayName', 'H_{p,max}');
    
%     legend('v_s_{1g}', 'v_s_{min}', 'Ma_{BO}', 'SR_{max}', 'SEP_{max}', 'SET_{max}', 'SE_{max}', 'v_{MO}', 'v_{min}', 'v_{max}', 'DP', 'H_{ops_{max}}', 'H_{p,max}', 'FontSize',15 ,Location='eastoutside');
    legend('FontSize',15 ,Location='eastoutside');

    hold off




end

