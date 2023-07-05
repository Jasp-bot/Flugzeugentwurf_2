clc
clear all
close all

Berechnung_FE2_PS4_Widerstand;


load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat
load Ergebnisse_Start_Landeanforderungen.mat;
load Ergebnisse_Widerstand_FE2.mat; 



%% Plot design

figure(1)

hold on
grid on

xlim([0, 0.05])
ylim([0, 1])

sz = size(Ergebnisse_Widerstand_FE2.x_vector);
numPlots = sz(1,1); %6;     % muss veraendert werden um off Design noch zu plotten


% Farbverlauf definieren
colorStart = [0, 0, 1];   % Startfarbe (RGB)
colorEnd = [0, 0, 0];     % Endfarbe (RGB)

% Farbwerte für jeden Plot berechnen
%colors = zeros(numPlots/2, 3);
for n_color = 1:(numPlots)
    colors(n_color, :) = colorStart + (n_color-1) * (colorEnd - colorStart) / ((numPlots)-1);
end

colors = vertcat(colors, colors);


% Linienarten definieren
lineStyles = {'--', ':', '-.', '-', '--', ':', '-.', '-', '--', ':', '-.', '-' }; % Gestrichelt, Gepunktet


% Plots erstellen
for n_plot = 1:numPlots
   autoplot(n_plot,1) = plot((Ergebnisse_Widerstand_FE2.x_vector_sum(n_plot,:)), Ergebnisse_Widerstand_FE2.c_A_F, 'LineStyle', lineStyles{1,n_plot}); % , 'Color', colors(n_plot, :), 'LineStyle', lineStyles{1,n_plot});
end
autoplot(numPlots+1, 1) = plot(Ergebnisse_Widerstand_FE2.x_vector_sum_off_D(numPlots,:), Ergebnisse_Widerstand_FE2.c_A_F_off_D,'red'); % Plot Off_Design


% Schnittpunkt Design / Off_design
schnittpunkt_des_off_D=InterX([Ergebnisse_Widerstand_FE2.x_vector_sum_off_D(numPlots,:); Ergebnisse_Widerstand_FE2.c_A_F_off_D], [(Ergebnisse_Widerstand_FE2.x_vector_sum(numPlots,:)); Ergebnisse_Widerstand_FE2.c_A_F]);
autoplot(numPlots+2, 1) = plot(schnittpunkt_des_off_D(1,1), schnittpunkt_des_off_D(2,1), '*b');


% Schnittpunkt Offdesign c_A_CR
schnittpunkt_off_D_c_A_CR = InterX([Ergebnisse_Widerstand_FE2.x_vector_sum_off_D(numPlots,:); Ergebnisse_Widerstand_FE2.c_A_F_off_D], [[0, 1]; [Ergebnisse_stat_Flaechenbelastung.C_A_CR, Ergebnisse_stat_Flaechenbelastung.C_A_CR]]);
autoplot(numPlots+3, 1) = plot(schnittpunkt_off_D_c_A_CR(1,1), schnittpunkt_off_D_c_A_CR(2,1), 'og');






% schoen machen des Plots 

legend(autoplot([1:n_plot+3]),{'+ SLW', '+ HLW', '+ Interferenz', '+ Rumpf', '+ Triebwerk', '+ Trimmung', '+ Abwind', '+ Profil', '+ ind. Widerstand', '+ Wellenwiderstand', '+Off-Design', 'Schnittpunkt von Design / Off-Design', 'C_{A,CR}'},...
    'Location','southeast','FontSize',18);

title('Kumulative Widerstandspolare')
xlabel('c_W');
ylabel('c_A');

hold off;

% % atomatisches Plotten
figure(2)

hold on
grid on
xlim([0, 0.05])
ylim([0, 1])

% Anzahl der Plots festlegen
numPlots = 10 %6;     % muss veraendert werden um off Design noch zu plotten

% Vector mit zu plottenden Werten


%x_vector = [c_W_SLW; c_W_HLW; c_w_int_fs; c_w_R; c_w_TW; c_w_trim; delta_c_w_H; c_w_p; c_w_ind; delta_c_WM; off_D];

% Farbverlauf definieren
colorStart = [0, 1, 0];   % Startfarbe (RGB)
colorEnd = [0, 0, 0];     % Endfarbe (RGB)

% Farbwerte für jeden Plot berechnen
%colors = zeros(numPlots/2, 3);
for n_color = 1:(numPlots/2)
    colors(n_color, :) = colorStart + (n_color-1) * (colorEnd - colorStart) / ((numPlots/2)-1);
end

colors = vertcat(colors, colors);



% Linienarten definieren
lineStyles = {'--', ':', '-.', '-', '--', ':', '-.', '-', '--', ':', '-.', '-' }; % Gestrichelt, Gepunktet

% % Vector mit zu plottenden Werten
% 
% x_vector = [c_W_SLW; c_W_HLW; c_w_int_fs; c_w_R; c_w_TW; c_w_trim; delta_c_w_H; c_w_p; c_w_ind; delta_c_WM];
% 
% for n_vec = 1:numPlots
%     if n_vec == 1
%         x_vector_sum(n_vec,:) = x_vector(n_vec,:);
%     else
%         x_vector_sum(n_vec,:) = x_vector_sum(n_vec-1,:) + x_vector(n_vec,:);
%     end
% end
% 
% % Offdesign Vector
% off_D_vector = [c_W_SLW_off_D; c_W_HLW_off_D; c_w_int_fs_off_D; c_w_R_off_D; c_w_TW_off_D; c_w_trim_off_D; delta_c_w_H_off_D; c_w_p_off_D; c_w_ind_off_D; delta_c_WM_off_D];
% 
% for n_vec_off_D = 1:numPlots
%     if n_vec_off_D == 1
%         x_vector_sum_off_D(n_vec_off_D,:) = off_D_vector(n_vec_off_D,:);
%     else
%         x_vector_sum_off_D(n_vec_off_D,:) = x_vector_sum_off_D(n_vec_off_D - 1,:) + off_D_vector(n_vec_off_D,:);
%     end
% end


% Plots erstellen

atoplots = cell(numPlots+1, 1);

% Plots erstellen
for n_plot = 1:numPlots
   autoplot(n_plot,1) = plot((Ergebnisse_Widerstand_FE2.x_vector_sum(n_plot,:)), Ergebnisse_Widerstand_FE2.c_A_F, 'LineStyle', lineStyles{1,n_plot}, 'Color', colors(n_plot, :), 'LineStyle', lineStyles{1,n_plot});
end
autoplot(numPlots+1, 1) = plot(Ergebnisse_Widerstand_FE2.x_vector_sum_off_D(numPlots,:), Ergebnisse_Widerstand_FE2.c_A_F_off_D,'red'); % Plot Off_Design


% Schnittpunkt Design / Off_design
schnittpunkt_des_off_D=InterX([Ergebnisse_Widerstand_FE2.x_vector_sum_off_D(numPlots,:); Ergebnisse_Widerstand_FE2.c_A_F_off_D], [(Ergebnisse_Widerstand_FE2.x_vector_sum(numPlots,:)); Ergebnisse_Widerstand_FE2.c_A_F]);
autoplot(numPlots+2, 1) = plot(schnittpunkt_des_off_D(1,1), schnittpunkt_des_off_D(2,1), '*b');


% Schnittpunkt Offdesign c_A_CR
schnittpunkt_off_D_c_A_CR = InterX([Ergebnisse_Widerstand_FE2.x_vector_sum_off_D(numPlots,:); Ergebnisse_Widerstand_FE2.c_A_F_off_D], [[0, 1]; [Ergebnisse_stat_Flaechenbelastung.C_A_CR, Ergebnisse_stat_Flaechenbelastung.C_A_CR]]);
autoplot(numPlots+3, 1) = plot(schnittpunkt_off_D_c_A_CR(1,1), schnittpunkt_off_D_c_A_CR(2,1), 'og');


% schoen machen des Plots 

legend(autoplot([1:11]),{'+ SLW', '+ HLW', '+ Interferenz', '+ Rumpf', '+ Triebwerk', '+ Trimmung', '+ Abwind', '+ Profil', '+ ind. Widerstand', '+ Wellenwiderstand', '+ Off-Design'},...
    'Location','southeast','FontSize',18);
title('Kumulative Widerstandspolare')
xlabel('c_W');
ylabel('c_A');

hold off;
% 
% 
% 
% 
% 
