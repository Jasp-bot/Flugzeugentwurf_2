% Visualisierung Fluegel Tank Neutralpunkte
clc
clear all
close all

%% Ausfuerung des Programms

Berechnung_PS8_Fluegel_Tank

%% Laden der Dateien

load Ergebnisse_Fluegel_Tank_NP.mat
load Projekt_specs.mat

%%%%%%%%%%%%%%% Visualisierung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %% Buffet Onset Diagramm
% 
% figure(1)
% hold on
% grid on
% %wert= 0:0.001:0.85;
% %f = @(x) (1459.9 .* x.^4 - 4203.5 .* x.^3 + 4483.4 .* x.^2 - 2100.9 .*x + 366.16);
% % Test_CA_BO = c_A_Bo >= c_A_ICA_MO; % Alle werte für Buffet onset sind unter c_A_ICA_MO
% p_bo(1) = plot(BO.wert, BO.f(BO.wert));
% p_bo(2) = plot([0.65 0.85],[BO.c_A_Bo_plot, BO.c_A_Bo_plot],'k');
% p_bo(3) = plot(BO.Ma_BO, BO.c_A_Bo_plot, 'or');
% ylim([0.5 1.2])
% xlim([0.65 0.85])
% title('Buffet Onset-Diagramm','FontSize',20)
% legend(p_bo([1:3]), {'Buffet Onset-Grenze', 'C_{A BO}','Schnittpunkt'},'FontSize',18);
% xlabel('Ma_{BO}','FontSize',16)
% ylabel('C_{A BO}','FontSize',16)
% hold off
% 
% 
% %% Streckung
% 
% vector_0_70_grad = 0:1:70;
% 
% figure(2)
% hold on 
% grid on 
% plot(Ergebnisse_Fluegel.f_streckung(vector_0_70_grad),'blue');
% plot(rad2deg(Ergebnisse_Fluegel.phi_25_max), Ergebnisse_Fluegel.streckung_phi25_max,'*r');
% xlabel('Pfeilung Phi_{25} in °','FontSize',16);
% ylabel('Streckung in 1','FontSize',16);
% ylim([0 40]);
% xlim([0 70]);
% title('Pitch Up-Limit','FontSize',20);
% lgd_PU =legend('Pitch Up-Limit','Gewählte Streckung','FontSize',18);
% ldg_PU.FontSize = 100;
% hold off
% 
% %% Visualisierung Flügel
% 
% 
% figure(3)
% hold on 
% grid on
% xlim([-1 (Ergebnisse_Fluegel.b/2 +5)])
% ylim([0 25])
% 
% % Neutralpunkte
% 
% 
% p(1) = plot((Ergebnisse_Fluegel.b/2) - DT.s_R - DT.s_I - NP.y_SP_A, (NP.x_NP_A),'xb');
% 
% p(2) = plot((Ergebnisse_Fluegel.b/2) - DT.s_R - NP.y_SP_I , DT.l_i_R + NP.versatz_HK - NP.x_NP_I,'*b');
% 
% p(3) = plot((Ergebnisse_Fluegel.b/2 -NP.y_SP_R),...
%     ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf))-NP.x_NP_R,'ob');
% 
% p(4) = plot((Ergebnisse_Fluegel.b/2) - NP.y_SP_ges,  NP.x_NP_ges,'*r'); %(versatz_HK + l_i_R) -
% 
% p(17)= plot((Ergebnisse_Fluegel.b/2) - NP.y_SP_ges, NP.x_SP_ges,'or');
% 
% % doppeltrapez
%  
% %DT.l_a
% p(5) = plot([0, 0],[0, DT.l_a],'black');
% %vorderkante 
% p(6) = plot([0 Ergebnisse_Fluegel.b/2-specs.R_rumpf],...
%     [DT.l_a, ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a)-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)]...
%     , 'black');
% %A_3_dt
% p(7) = plot([(Ergebnisse_Fluegel.b/2-specs.R_rumpf) Ergebnisse_Fluegel.b/2]...
%     , [((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)),...
%     ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf))]...
%     ,'black' );
% %hinterkante
% p(8) = plot([0 ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))]...
%     , [0 (tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink))))]...
%     ,'black');  % (tan(DT.phi_HK_dt)*Ergebnisse_Fluegel.b/2)-(tan(DT.phi_HK_dt) * (Ergebnisse_Fluegel.b/2) * DT.y_kink)
% % innenkante 
% p(9) = plot([Ergebnisse_Fluegel.b/2 Ergebnisse_Fluegel.b/2],...
%     [(tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))),...
%     ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf))],...
%     'black');
% %kink
% p(10) = plot([((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink))) Ergebnisse_Fluegel.b/2],...
%     [(tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))),...
%     (tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink))))],...
%     'black');
% % &Rumpf
% p(11) = plot([(Ergebnisse_Fluegel.b/2 - specs.R_rumpf), (Ergebnisse_Fluegel.b/2 - specs.R_rumpf)]...
%     ,[0 25], '--k');
% p(12) = plot([(Ergebnisse_Fluegel.b/2 + specs.R_rumpf), (Ergebnisse_Fluegel.b/2 + specs.R_rumpf)]...
%     ,[0 25], '--k');
% 
% % plot von l_i_A
% p(13) = plot([DT.s_A DT.s_A],[NP.versatz_HK, NP.versatz_HK+DT.l_i_A],'--k' );
% % plot phi_25_A
% p(14) = plot([0 DT.s_A],[DT.l_a*0.75 (NP.versatz_HK + DT.l_i_A*0.75 )],'--b');
% % plot phi_25_I
% p(15) = plot([DT.s_A (DT.s_A + DT.s_I)],[(NP.versatz_HK + DT.l_i_A*0.75 ) (NP.versatz_HK + DT.l_i_I*0.75 )],'--b');
% % plot phi_25_R
% p(16) = plot([(Ergebnisse_Fluegel.b/2 - specs.R_rumpf) (Ergebnisse_Fluegel.b/2)],...
%     [(NP.versatz_HK + DT.l_a_R*0.75) (NP.versatz_HK + DT.l_i_R*0.75 )],'--b');
% 
% legend(p([1:4 17 16]),{'Neutralpunkt F_A', 'Neutralpunkt F_I', 'Neutralpunkt F_R', 'Gesamtneutralpunkt','Gesamtschwerpunkt','Phi_{25}- Linie'},'Location','northwest','FontSize',18);
% title('Veranschaulichung Flügelansicht mit Neutralpunkten','FontSize',20);
% xlabel('Meter','FontSize',16)
% ylabel('Meter','FontSize',16)

% hinterkante veranschaulichung
%plot([0 Ergebnisse_Fluegel.b/2],[0 (tan(DT.phi_HK_dt)*Ergebnisse_Fluegel.b/2)],'-.g');

% 
% %% Trapezflügel
% figure(4)
% hold on 
% grid on
%  
% %l_a 
% plot([0 0],[0 ET.l_a],'b');
% %vorderkante
% plot([0 Ergebnisse_Fluegel.b/2], [ET.l_a ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2)+ET.l_a)], 'b');
% %hinterkante
% plot([0 Ergebnisse_Fluegel.b/2], [0 (tan(ET.phi_HK_a)*Ergebnisse_Fluegel.b/2)],'b');
% %innenkante 
% plot([Ergebnisse_Fluegel.b/2 Ergebnisse_Fluegel.b/2],...
%     [(tan(ET.phi_HK_a)*Ergebnisse_Fluegel.b/2) ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2)+ET.l_a)], 'b')
% 
% % &Rumpf
% pt(11) = plot([(Ergebnisse_Fluegel.b/2 - specs.R_rumpf) (Ergebnisse_Fluegel.b/2 - specs.R_rumpf)],[0 25], '--k');
% pt(12) = plot([(Ergebnisse_Fluegel.b/2 + specs.R_rumpf) (Ergebnisse_Fluegel.b/2 + specs.R_rumpf)],[0 25], '--k');
% 
% title('Einfachtrapezflügel','FontSize',20)
% xlabel('Meter','FontSize',16)
% ylabel('Meter','FontSize',16)
% xlim([-1 40])
% hold off

%% Plot Tank

% figure(5)
% hold on 
% grid on
% xlim([-1 (Ergebnisse_Fluegel.b/2+5)])
% ylim([0 25])
% 
% %doppeltrapez
%  
% %l_a_dt
% p2(1) = plot([0, 0],[0, DT.l_a],'black');
% %vorderkante 
% p2(2) = plot([0 Ergebnisse_Fluegel.b/2-specs.R_rumpf],...
%     [DT.l_a, ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a)-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)],...
%      'black');
% %A_3_dt
% p2(3) = plot([(Ergebnisse_Fluegel.b/2-specs.R_rumpf) Ergebnisse_Fluegel.b/2]...
%     , [((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)),...
%     ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf))]...
%     ,'black' );
% %hinterkante
% p2(4) = plot([0 ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))]...
%     , [0 (tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink))))]...
%     ,'black');  % (tan(DT.phi_HK_dt)*Ergebnisse_Fluegel.b/2)-(tan(DT.phi_HK_dt) * (Ergebnisse_Fluegel.b/2) * DT.y_kink)
% % innenkante 
% p2(5) = plot([Ergebnisse_Fluegel.b/2 Ergebnisse_Fluegel.b/2],...
%     [(tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))),...
%     ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf))],...
%     'black');
% %kink
% p2(6) = plot([((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink))) Ergebnisse_Fluegel.b/2],...
%     [(tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))),...
%     (tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink))))],...
%     'black');
% %Rumpf
% p2(7) = plot([(Ergebnisse_Fluegel.b/2 - specs.R_rumpf), (Ergebnisse_Fluegel.b/2 - specs.R_rumpf)]...
%     ,[0 25], '--k');
% p2(8) = plot([(Ergebnisse_Fluegel.b/2 + specs.R_rumpf), (Ergebnisse_Fluegel.b/2 + specs.R_rumpf)]...
%     ,[0 25], '--k');
% 
% 
% % plot Vorderer Holm
% p2(9) = plot([0.05*(Ergebnisse_Fluegel.b/2) Ergebnisse_Fluegel.b/2-specs.R_rumpf],...
%     [(tan(DT.phi_HK_dt)*0.05*(Ergebnisse_Fluegel.b/2))+Ergebnisse_Fluegel.Fluegeltiefen_eta(1,95)*(1-0.15),...
%     ((tan(DT.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a)-(tan(DT.phi_VK_max)*specs.R_rumpf)-(DT.l_i_R*0.15)], 'blue');
% 
% % plot Hinterkante Holm außen
% p2(10) = plot([0.05*(Ergebnisse_Fluegel.b/2), ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))],...
%     [(1-0.65)*Ergebnisse_Fluegel.Fluegeltiefen_eta(1,95)+(tan(DT.phi_HK_dt)*0.05*(Ergebnisse_Fluegel.b/2)),...
%     ((tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))))+(1-0.65)*DT.l_i_A],...
%     'blue'); 
% 
% % plot Hinterkante Holm innen
% p2(11) = plot([((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink))), (Ergebnisse_Fluegel.b/2)-specs.R_rumpf],...
%     [((tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))))+(1-0.65)*DT.l_i_A,...
%     (tan(DT.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(DT.phi_VK_max)*specs.R_rumpf)-(DT.l_i_R*0.65)],'blue');
% 
% % plot Vordekante Holm rumpf
% p2(12) = plot([(Ergebnisse_Fluegel.b/2-specs.R_rumpf) Ergebnisse_Fluegel.b/2],...
%     [((tan(DT.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a)-(tan(DT.phi_VK_max)*specs.R_rumpf)-(DT.l_i_R*0.15),...
%     ((tan(DT.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a)-(tan(DT.phi_VK_max)*specs.R_rumpf)-(DT.l_i_R*0.15)],...
%     'blue' );
% p2(13) = plot([(Ergebnisse_Fluegel.b/2-specs.R_rumpf) Ergebnisse_Fluegel.b/2],...
%     [((tan(DT.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a)-(tan(DT.phi_VK_max)*specs.R_rumpf)-(DT.l_i_R*0.65),...
%     ((tan(DT.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a)-(tan(DT.phi_VK_max)*specs.R_rumpf)-(DT.l_i_R*0.65)],...
%     'blue');
% 
% % plot außenkante Tank
% p2(14) = plot([0.05*(Ergebnisse_Fluegel.b/2) 0.05*(Ergebnisse_Fluegel.b/2)],...
%     [(tan(DT.phi_HK_dt)*0.05*(Ergebnisse_Fluegel.b/2))+Ergebnisse_Fluegel.Fluegeltiefen_eta(1,95)*(1-0.65),...
%     (tan(DT.phi_HK_dt)*0.05*(Ergebnisse_Fluegel.b/2))+Ergebnisse_Fluegel.Fluegeltiefen_eta(1,95)*(1-0.15)],...
%     'blue');
% 
% 
% title('Veranschaulichung Flügelansicht mit Tank','FontSize',20);
% xlabel('Meter','FontSize',16)
% ylabel('Meter','FontSize',16)




figure(6)
hold on 
grid on
xlim([-1 (Ergebnisse_Fluegel.b/2 +5)])
ylim([0 25])

% Neutralpunkte


p(1) = plot((Ergebnisse_Fluegel.b/2) - DT.s_R - DT.s_I - NP.y_SP_A, (NP.x_NP_A),'xb');

p(2) = plot((Ergebnisse_Fluegel.b/2) - DT.s_R - NP.y_SP_I , DT.l_i_R + NP.versatz_HK - NP.x_NP_I,'*b');

p(3) = plot((Ergebnisse_Fluegel.b/2 -NP.y_SP_R),...
    ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf))-NP.x_NP_R,'ob');

p(4) = plot((Ergebnisse_Fluegel.b/2) - NP.y_SP_ges,  NP.x_NP_ges,'*r'); %(versatz_HK + l_i_R) -

p(17)= plot((Ergebnisse_Fluegel.b/2) - NP.y_SP_ges, NP.x_SP_ges,'or');

% doppeltrapez
 
%DT.l_a
p(5) = plot([0, 0],[0, DT.l_a],'black');
%vorderkante 
p(6) = plot([0 Ergebnisse_Fluegel.b/2-specs.R_rumpf],...
    [DT.l_a, ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a)-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)]...
    , 'black');
zp(1) = plot([Ergebnisse_Fluegel.b/2-specs.R_rumpf, Ergebnisse_Fluegel.b/2],...
    [NP.versatz_HK + DT.l_i_R, (NP.versatz_HK + DT.l_i_R + tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)], 'green');
%A_3_dt
p(7) = plot([(Ergebnisse_Fluegel.b/2-specs.R_rumpf) Ergebnisse_Fluegel.b/2]...
    , [((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf)),...
    ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf))]...
    ,'black' );
%hinterkante
p(8) = plot([0 ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))]...
    , [0 (tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink))))]...
    ,'black');  % (tan(DT.phi_HK_dt)*Ergebnisse_Fluegel.b/2)-(tan(DT.phi_HK_dt) * (Ergebnisse_Fluegel.b/2) * DT.y_kink)
% innenkante 
p(9) = plot([Ergebnisse_Fluegel.b/2 Ergebnisse_Fluegel.b/2],...
    [(tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))),...
    ((tan(Ergebnisse_Fluegel.phi_VK_max) * Ergebnisse_Fluegel.b/2) + DT.l_a-(tan(Ergebnisse_Fluegel.phi_VK_max)*specs.R_rumpf))],...
    'black');
%kink
p(10) = plot([((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink))) Ergebnisse_Fluegel.b/2],...
    [(tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink)))),...
    (tan(DT.phi_HK_dt) * ((Ergebnisse_Fluegel.b/2)*(1-(DT.y_kink))))],...
    'black');
% &Rumpf
p(11) = plot([(Ergebnisse_Fluegel.b/2 - specs.R_rumpf), (Ergebnisse_Fluegel.b/2 - specs.R_rumpf)]...
    ,[0 25], '--k');
p(12) = plot([(Ergebnisse_Fluegel.b/2 + specs.R_rumpf), (Ergebnisse_Fluegel.b/2 + specs.R_rumpf)]...
    ,[0 25], '--k');

% plot von l_i_A
p(13) = plot([DT.s_A DT.s_A],[NP.versatz_HK, NP.versatz_HK+DT.l_i_A],'--k' );
% plot phi_25_A
p(14) = plot([0 DT.s_A],[DT.l_a*0.75 (NP.versatz_HK + DT.l_i_A*0.75 )],'--b');
% plot phi_25_I
p(15) = plot([DT.s_A (DT.s_A + DT.s_I)],[(NP.versatz_HK + DT.l_i_A*0.75 ) (NP.versatz_HK + DT.l_i_I*0.75 )],'--b');
% plot phi_25_R
p(16) = plot([(Ergebnisse_Fluegel.b/2 - specs.R_rumpf) (Ergebnisse_Fluegel.b/2)],...
    [(NP.versatz_HK + DT.l_a_R*0.75) (NP.versatz_HK + DT.l_i_R*0.75 )],'--b');

legend(p([1:4 17 16]),{'Neutralpunkt F_A', 'Neutralpunkt F_I', 'Neutralpunkt F_R', 'Gesamtneutralpunkt','Gesamtschwerpunkt','Phi_{25}- Linie'},'Location','northwest','FontSize',18);
title('Veranschaulichung Flügelansicht mit Neutralpunkten für Schwerpunktberechnung','FontSize',20);
xlabel('Meter','FontSize',16)
ylabel('Meter','FontSize',16)




