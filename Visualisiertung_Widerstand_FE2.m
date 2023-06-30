clc
clear all
close all

load Ergebnisse_Start_Landeanforderungen.mat;
load Ergebnisse_Widerstand_FE2.mat; 

%% Plot design

% figure(1)
% hold on
% grid on
% xlim([0, 0.05])
% ylim([0, 1])
% 
% 
% 
% % plot Leitwerke SLW und HLW
% 
% p(1) = plot(c_W_SLW, n_iteration_vec, '.-b');   % SLW
% 
% p(2) = plot(c_W_HLW, n_iteration_vec, '.-c');   % HLW
% 
% % Interferenzwiderstand
% 
% p(3) =  plot(c_w_int_fs, n_iteration_vec, '.green');
% 
% % Plot Rumpfwiderstand
% p(4) = plot(c_w_R, n_iteration_vec, '-k');
% 
% 
% % Plot Widerstand Triebwerk
% p(5) = plot(c_w_TW,n_iteration_vec, '-m');
% 
% 
% % Plot Trimwiderstand HLW
% 
% p(6) = plot(c_w_trim, n_iteration_vec,'-blue');
% 
% % Abwindwiderstand
% p(7) = plot(delta_c_w_H, n_iteration_vec, 'c');
% 
% % Plot Profilwiderstand 
% 
% p(8) = plot(c_w_p, n_iteration_vec, '-.k');
% 
% % plot Induzierter Widerstand
% p(9) = plot(c_w_ind, n_iteration_vec, '-red');
% 
% % Plot Transsonischer Widersatnd
% p(10) = plot(delta_c_WM, n_iteration_vec, '-green');
%  
% 
% % Testplot Widerstände aufaddiert
% widerstaende_aufaddiert = c_w_ind + delta_c_WM + c_w_R + c_w_TW + c_w_trim + c_w_int_fs + delta_c_w_H + c_W_HLW + c_W_SLW;
% plot(widerstaende_aufaddiert, n_iteration_vec, '*r')
% 
% 
% legend(p([1:10]),{'+ SLW', '+ HLW', '+ Interferenz', '+ Rumpf', '+ Triebwerk', '+ Trimmung', '+ Abwind', '+ Profil', '+ ind. Widerstand', '+ Wellenwiderstand'},'Location','southeast','FontSize',18);
% title('Kumulative Widerstandspolare')
% xlabel('c_{W}');
% ylabel('c_A');


% atomatisches Plotten
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

% Vector mit zu plottenden Werten

x_vector = [c_W_SLW; c_W_HLW; c_w_int_fs; c_w_R; c_w_TW; c_w_trim; delta_c_w_H; c_w_p; c_w_ind; delta_c_WM];

for n_vec = 1:numPlots
    if n_vec == 1
        x_vector_sum(n_vec,:) = x_vector(n_vec,:);
    else
        x_vector_sum(n_vec,:) = x_vector_sum(n_vec-1,:) + x_vector(n_vec,:);
    end
end

% Offdesign Vector
off_D_vector = [c_W_SLW_off_D; c_W_HLW_off_D; c_w_int_fs_off_D; c_w_R_off_D; c_w_TW_off_D; c_w_trim_off_D; delta_c_w_H_off_D; c_w_p_off_D; c_w_ind_off_D; delta_c_WM_off_D];

for n_vec_off_D = 1:numPlots
    if n_vec_off_D == 1
        x_vector_sum_off_D(n_vec_off_D,:) = off_D_vector(n_vec_off_D,:);
    else
        x_vector_sum_off_D(n_vec_off_D,:) = x_vector_sum_off_D(n_vec_off_D - 1,:) + off_D_vector(n_vec_off_D,:);
    end
end


% Plots erstellen

atoplots = cell(numPlots, 1);

for n_plot = 1:numPlots
   autoplot(n_plot,1) = plot((x_vector_sum(n_plot,:)), n_iteration_vec, 'Color', colors(n_plot, :), 'LineStyle', lineStyles{1,n_plot});
end
autplot(numPlots+1,1)= plot(x_vector_sum_off_D(numPlots,:),c_A_F_off_D_vec,'red');


% schoen machen des Plots 

legend(autoplot([1:10]),{'+ SLW', '+ HLW', '+ Interferenz', '+ Rumpf', '+ Triebwerk', '+ Trimmung', '+ Abwind', '+ Profil', '+ ind. Widerstand', '+ Wellenwiderstand'},...
    'Location','southeast','FontSize',18);
title('Kumulative Widerstandspolare')
xlabel('c_W');
ylabel('c_A');

hold off;



intersections_1 = InterX([startschub.s1; startschub.S0_GTo_To(1,:)],[startschub.s1; startschub.S0_GTo_CL(3,:)]);

schnittpunkte=InterX([x_vector_sum_off_D(numPlots,:); c_A_F_off_D_vec], [(x_vector_sum(numPlots,:)); n_iteration_vec]);