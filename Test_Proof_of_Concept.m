% Test und Proof of concept
clc
clear all
close all
% 
% 
% load Projekt_specs.mat;
% load Ergebnisse_ISA_DATA.mat;
% load Ergebnisse_Widerstand_FE2.mat;
% load Getroffene_Annahmen_und_FUN.mat;
% load Ergebnisse_Fluegel_Tank_NP.mat;
% load Ergebnisse_Leitwerke.mat;
% 
% 
% addpath('Unterfunktionen Widerstand');
% 
% 
% hoehe_CR = round(unitsratio('m','ft')*specs.flight_level*10^2);
% 
% %v_air = linspace(50,300, 50).';
% 
% v_air = ones(50,1) .* (specs.Ma_CR * ISA.a(hoehe_CR)); 
% 
% %v_air_off_D = 
% Machzahl = v_air ./ ISA.a(hoehe_CR);
% c_A_F = linspace(0,1,50);
% 
% 
% Abwindfaktor = 1.75 * ((Annahmen.c_A_alpha_F)/(pi * Ergebnisse_Fluegel.streckung_phi25_max *...
%     (Ergebnisse_Fluegel.lambda * (HLW.r/(Ergebnisse_Fluegel.b/2))^0.25 *...
%     (1+ (abs(Annahmen.z_abstand))/((Ergebnisse_Fluegel.b/2))) )));
% 
% 
% 
% c_w_ind = Induzierter_W(c_A_F);
% c_w_int_fs = Interferenz_W(v_air);
% [c_A_H, c_A_ges, c_w_trim] = Leitwerke_Trim_W(c_A_F);
% [c_w_HLW_min, c_w_SLW_min] = Leitwerke_W(v_air); % muss noch integral geloest werden
% 
% % nur fuer des
% % for n_stuetz = 1:50
% % 
% % c_W_HLW(n_stuetz,:) = trapz(c_w_HLW_min(n_stuetz,:))*10^(-2);
% % c_W_SLW(n_stuetz,:) = trapz(c_w_SLW_min(n_stuetz,:))*10^(-2);
% % 
% % end
% 
% [c_w_p, c_w_p_min_Re] = Profilwiderstand(v_air);
% 
% 
% 
% [c_w_R, alpha_Rumpf_grad] = Rumpfwiderstand(Machzahl, Abwindfaktor, c_A_ges);
% [delta_c_WM, delta_Ma] = Transsonischer_W(Machzahl, c_A_F);
% [c_w_TW] = Triebwerkswiderstand(v_air, alpha_Rumpf_grad);
% [delta_c_w_H] = Zusatz_W(c_A_F, Abwindfaktor, c_A_H);
% 
% 
% % nur fuer des
% if length(v_air)> 1;
%     for n_stuetz = 1:50
%     
%     c_W_HLW(n_stuetz,:) = trapz(c_w_HLW_min(n_stuetz,:))*10^(-2);
%     c_W_SLW(n_stuetz,:) = trapz(c_w_SLW_min(n_stuetz,:))*10^(-2);
%     
%     end
% else length(v_air) == 1;
%     for n_stuetz = 1:50
%     
%     c_W_HLW(n_stuetz,:) = trapz(c_w_HLW_min)*10^(-2);
%     c_W_SLW(n_stuetz,:) = trapz(c_w_SLW_min)*10^(-2);
%     c_w_p(n_stuetz,:) = c_w_p(1,1);
%     c_w_int_fs(n_stuetz,:) = c_w_int_fs(1,1);
% 
%     end
% end
% 
% 
% %% atomatisches Plotten
% figure(2)
% 
% hold on
% grid on
% xlim([0, 0.05])
% ylim([0, 1])
% 
% % Anzahl der Plots festlegen
% numPlots = 10 %6;     % muss veraendert werden um off Design noch zu plotten
% 
% % Vector mit zu plottenden Werten
% 
% 
% %x_vector = [c_W_SLW.'; c_W_HLW.'; c_w_int_fs.'; c_w_R.'; c_w_TW.'; c_w_trim; delta_c_w_H; c_w_p.'; c_w_ind; delta_c_WM.']; % off_D];
% 
% % Farbverlauf definieren
% colorStart = [0, 1, 0];   % Startfarbe (RGB)
% colorEnd = [0, 0, 0];     % Endfarbe (RGB)
% 
% % Farbwerte für jeden Plot berechnen
% %colors = zeros(numPlots/2, 3);
% for n_color = 1:ceil(numPlots/2)
%     colors(n_color, :) = colorStart + (n_color-1) * (colorEnd - colorStart) / ((numPlots/2)-1);
% end
% 
% colors = vertcat(colors, colors);
% 
% 
% 
% % Linienarten definieren
% lineStyles = {'--', ':', '-.', '-', '--', ':', '-.', '-', '--', ':', '-.', '-' }; % Gestrichelt, Gepunktet
% 
% % Vector mit zu plottenden Werten
% % 
% x_vector = [c_W_SLW.'; c_W_HLW.'; c_w_int_fs.'; diag(c_w_R).'; c_w_TW.'; c_w_trim; delta_c_w_H; c_w_p.'; c_w_ind; delta_c_WM.']; % off_D];
% 
% for n_vec = 1:numPlots
%     if n_vec == 1
%         x_vector_sum(n_vec,:) = x_vector(n_vec,:);
%     else
%         x_vector_sum(n_vec,:) = x_vector_sum(n_vec-1,:) + x_vector(n_vec,:);
%     end
% end
% % 
% % % Offdesign Vector
% % off_D_vector = [c_W_SLW_off_D; c_W_HLW_off_D; c_w_int_fs_off_D; c_w_R_off_D; c_w_TW_off_D; c_w_trim_off_D; delta_c_w_H_off_D; c_w_p_off_D; c_w_ind_off_D; delta_c_WM_off_D];
% % 
% % for n_vec_off_D = 1:length(off_D_vector)
% %     if n_vec_off_D == 1
% %         x_vector_sum_off_D(n_vec_off_D,:) = off_D_vector(n_vec_off_D,:);
% %     else
% %         x_vector_sum_off_D(n_vec_off_D,:) = x_vector_sum_off_D(n_vec_off_D-1,:) + off_D_vector(n_vec_off_D,:);
% %     end
% % end
% % 
% % 
% % Plots erstellen
% 
% atoplots = cell(numPlots, 1);
% 
% for n_plot = 1:numPlots
%    autoplot(n_plot,1) = plot((x_vector_sum(n_plot,:)), c_A_F, 'Color', colors(n_plot, :), 'LineStyle', lineStyles{1,n_plot});
% end
% % plot(off_D_vector,c_A_F_off_D,'red')
% % 
% % 
% % % schoen machen des Plots 
% % 
% legend(autoplot([1:10]),{'+ SLW', '+ HLW', '+ Interferenz', '+ Rumpf', '+ Triebwerk', '+ Trimmung', '+ Abwind', '+ Profil', '+ ind. Widerstand', '+ Wellenwiderstand'},...
%     'Location','southeast','FontSize',18);
% title('Kumulative Widerstandspolare')
% xlabel('c_W');
% ylabel('c_A');
% % 
% % hold off;
% % 
% % 
% % 



% Y = [2 2 2]
% Z = [1 1.5 1]
% x = trapz(Y,2)

% Beispiel-Datenpunkte
x = [0, 1, 2, 3, 4];
y = [0, 1, 4, 9, 16];

% Berechnung der Fläche unter der Kurve mit trapz
area_under_curve = trapz(y)

hold on
plot(x)
plot(y)