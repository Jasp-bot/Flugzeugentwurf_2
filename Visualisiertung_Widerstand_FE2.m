clc
clear all
close all

Berechnung_FE2_PS4_Widerstand;


load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat
load Ergebnisse_Start_Landeanforderungen.mat;
load Ergebnisse_Widerstand_FE2.mat; 
load Getroffene_Annahmen_und_FUN.mat;


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
colors = zeros(numPlots/2, 3);
for n_color = 1:(numPlots)
    colors(n_color, :) = colorStart + (n_color-1) * (colorEnd - colorStart) / ((numPlots)-1);
end

% colors = vertcat(colors, colors);


% Linienarten definieren
lineStyles = {'--', '-', '-.', '-', '--', '-', '-.', '-', '--', '-', '-.', '-' }; % Gestrichelt, Gepunktet


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

legend(autoplot(1:n_plot+3),{'SLW', '+ HLW', '+ Interferenz', '+ Rumpf', '+ Triebwerk', '+ Trimmung', '+ Abwind', '+ Profil', '+ ind. Widerstand', '+ Wellenwiderstand', '+Off-Design', 'Schnittpunkt von Design / Off-Design', 'C_{A,CR}'},...
    'Location','southeast','FontSize',10)%15);

title('Kumulative Widerstandspolare','FontSize',25)
xlabel('c_W','FontSize',20);
ylabel('c_A','FontSize',20);

hold off;

% 
% figure(2)
% 
% hold on
% 
% Widerstandsanteile_CR_vec = Ergebnisse_Widerstand_FE2.x_vector(:,104);
% Gesamt_W=sum(Widerstandsanteile_CR_vec);
% 
% proz_Widerstandsanteile = Widerstandsanteile_CR_vec ./ Gesamt_W;
% labels = ["SLW", "HLW", "Interferenz", "Rumpf", "Triebwerk", "Trimmung", "Abwind", "Profil", "ind. Widerstand", "Wellenwiderstand"];
% pie(proz_Widerstandsanteile);     %, labels);
% 
% legend('+ SLW', '+ HLW', '+ Interferenz', '+ Rumpf', '+ Triebwerk', '+ Trimmung', '+ Abwind', '+ Profil', '+ ind. Widerstand', '+ Wellenwiderstand')
% 
% 
% hold off


%% Reziproke Gleitzahlpolare

[max_y, max_x] = max(Ergebnisse_Widerstand_FE2.Gleitverhaeltnis_Des);

figure(3)
hold on
grid on 

title('Reziproke Gleitzahlpolare','FontSize',25)
xlabel('c_A','FontSize',20);
ylabel('c_A / c_W','FontSize',20);
xlim([0 1])
ylim([0 25])


plot(Ergebnisse_Widerstand_FE2.c_A_F, Ergebnisse_Widerstand_FE2.Gleitverhaeltnis_Des);  % Designpolare
plot(Ergebnisse_Widerstand_FE2.c_A_F_off_D, Ergebnisse_Widerstand_FE2.Gleitverhaeltnis_off_D);      % Offdesignpolare
plot(Ergebnisse_Widerstand_FE2.c_A_F(1,max_x), max_y, '*k');        % Green Dot / optimales Gleitverhaeltnis

legend('Design', 'Offdesign', 'Maximales Gleitverhältnis', 'FontSize',25, Location='eastoutside')

hold off;



% %% Off D
% 
% figure(4)
% 
% hold on
% grid on
% 
% xlim([0, 0.05])
% ylim([0, 1])
% 
% sz = size(Ergebnisse_Widerstand_FE2.x_vector);
% numPlots = sz(1,1); %6;     % muss veraendert werden um off Design noch zu plotten
% 
% 
% % Farbverlauf definieren
% colorStart = [0, 0, 1];   % Startfarbe (RGB)
% colorEnd = [0, 0, 0];     % Endfarbe (RGB)
% 
% % Farbwerte für jeden Plot berechnen
% colors = zeros(numPlots/2, 3);
% for n_color = 1:(numPlots)
%     colors(n_color, :) = colorStart + (n_color-1) * (colorEnd - colorStart) / ((numPlots)-1);
% end
% 
% % colors = vertcat(colors, colors);
% 
% 
% % Linienarten definieren
% lineStyles = {'--', '-', '-.', '-', '--', '-', '-.', '-', '--', '-', '-.', '-' }; % Gestrichelt, Gepunktet
% 
% 
% % Plots erstellen
% for n_plot = 1:numPlots
%    autoplot_x(n_plot,1) = plot((Ergebnisse_Widerstand_FE2.x_vector_sum_off_D(n_plot,:)), Ergebnisse_Widerstand_FE2.c_A_F_off_D, 'LineStyle', lineStyles{1,n_plot}); % , 'Color', colors(n_plot, :), 'LineStyle', lineStyles{1,n_plot});
% end
% % autoplot(numPlots+1, 1) = plot(Ergebnisse_Widerstand_FE2.x_vector_sum_off_D(numPlots,:), Ergebnisse_Widerstand_FE2.c_A_F_off_D,'red'); % Plot Off_Design
% 
% 
% % % Schnittpunkt Design / Off_design
% % schnittpunkt_des_off_D=InterX([Ergebnisse_Widerstand_FE2.x_vector_sum_off_D(numPlots,:); Ergebnisse_Widerstand_FE2.c_A_F_off_D], [(Ergebnisse_Widerstand_FE2.x_vector_sum(numPlots,:)); Ergebnisse_Widerstand_FE2.c_A_F]);
% % autoplot(numPlots+2, 1) = plot(schnittpunkt_des_off_D(1,1), schnittpunkt_des_off_D(2,1), '*b');
% % 
% % 
% % % Schnittpunkt Offdesign c_A_CR
% % schnittpunkt_off_D_c_A_CR = InterX([Ergebnisse_Widerstand_FE2.x_vector_sum_off_D(numPlots,:); Ergebnisse_Widerstand_FE2.c_A_F_off_D], [[0, 1]; [Ergebnisse_stat_Flaechenbelastung.C_A_CR, Ergebnisse_stat_Flaechenbelastung.C_A_CR]]);
% % autoplot(numPlots+3, 1) = plot(schnittpunkt_off_D_c_A_CR(1,1), schnittpunkt_off_D_c_A_CR(2,1), 'og');
% % 
% % 
% % % schoen machen des Plots 
% % 
% % legend(autoplot(1:n_plot+3),{'SLW', '+ HLW', '+ Interferenz', '+ Rumpf', '+ Triebwerk', '+ Trimmung', '+ Abwind', '+ Profil', '+ ind. Widerstand', '+ Wellenwiderstand', '+Off-Design', 'Schnittpunkt von Design / Off-Design', 'C_{A,CR}'},...
% %     'Location','southeast','FontSize',10)%15);
% % 
% % title('Kumulative Widerstandspolare','FontSize',25)
% % xlabel('c_W','FontSize',20);
% % ylabel('c_A','FontSize',20);
% 
% hold off;


