% Visualisierung Leitwerke

clc
clear all
close all

%% Ausfuerung des Programms
Berechnung_PS9_Leitwerke;

%% Laden der Dateien

load Projekt_specs.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Leitwerke.mat;

%% Visualisierung Hoehenleitwerk

kabinenversatz = specs.coanlaenge - specs.HLW_beginn;

figure(1)

hold on
grid on
% % plot Rumpf
% centerline
p_rumpf(1) = plot([0 0],[0 15],'k'); % 15 entspricht ca der Conelaenge am Flugzeugheck
% kabinenbegrinn
p_rumpf(2) = plot([0 specs.R_rumpf],[15 15],'k');
% aussenwand Cone
p_rumpf(3) = plot([0 specs.R_rumpf],[0 15],'k');



% innenkante außenfluegel
p_HLW(1) = plot([HLW.rumpfversatz_gross HLW.rumpfversatz_gross],...
    [kabinenversatz kabinenversatz-HLW.l_i],'blue');

% VK aussenfluegel
p_HLW(2) = plot([HLW.rumpfversatz_gross HLW.b/2],...
    [kabinenversatz (kabinenversatz-tan(HLW.phi_VK)*(HLW.s_A))],'blue');

% aussenkante aussenfluegel
p_HLW(3) = plot([HLW.b/2 HLW.b/2],[(kabinenversatz-tan(HLW.phi_VK)*(HLW.s_A)),...
    (kabinenversatz-tan(HLW.phi_VK)*(HLW.s_A))-HLW.l_a],'blue');

% HK aussenfluegel bis l_i_A_HLW
p_HLW(4) = plot([HLW.b/2 HLW.rumpfversatz_gross],...
    [(kabinenversatz-tan(HLW.phi_VK)*(HLW.s_A))-HLW.l_a,...
    (((kabinenversatz-tan(HLW.phi_VK)*(HLW.s_A))-HLW.l_a)...
    +tan(HLW.phi_HK)*(HLW.s_A))],...
    'blue'); %%kabinenversatz-HLW.l_i



% Rumpfanteil

% oberkante rumpfanteil
p_HLW(5) = plot([0 HLW.R_rumpf_oben],[kabinenversatz kabinenversatz],'blue');
% innenkante an symmetrieebene
p_HLW(6) = plot([0 0],[kabinenversatz HLW.schnittpunkt_y],'blue');
% unterkante rumpfteil
p_HLW(7) = plot([0 HLW.schnittpunkt_x],[HLW.schnittpunkt_y HLW.schnittpunkt_y],'blue');
% Hinterkante Dreiecksteil
p_HLW(8) = plot([HLW.R_rumpf_oben HLW.schnittpunkt_x],...
    [kabinenversatz-HLW.l_i, (kabinenversatz - HLW.l_i + tan(HLW.phi_HK)*(HLW.R_rumpf_oben - HLW.schnittpunkt_x))],'blue'); 


% 
% % plot phi_25_HLW linie
R_rumpf_oben_phi25 = ((kabinenversatz-HLW.l_i_R*0.25)/specs.coanlaenge)*specs.R_rumpf;
    % Rumfteil
p_HLW(9) = plot([0, R_rumpf_oben_phi25],[(kabinenversatz-HLW.l_i_R*0.25), (kabinenversatz-HLW.l_i_R*0.25)],'--k');
    % Fluegelteil
p_HLW(10) = plot([R_rumpf_oben_phi25, HLW.b/2],[(kabinenversatz-HLW.l_i_R*0.25),...
    (kabinenversatz-tan(HLW.phi_VK)*(HLW.s_A))-HLW.l_a*0.25],'--k');
% % plot Neutralpunkt
p_HLW(11) = plot(HLW.y_SP_ges/2, HLW.x_NP_ges/2,'*red');


% % plot von Holmen

% 15% Holm
R_rumpf_oben_Holm = ((kabinenversatz-HLW.l_i_R*0.15)/specs.coanlaenge)*specs.R_rumpf;
% Rumpfanteil oben
p_HLW(12) = plot([0, R_rumpf_oben_Holm],[(kabinenversatz-HLW.l_i_R*0.15), (kabinenversatz-HLW.l_i_R*0.15)],'g');
% Außenfluegel oben
p_HLW(13) = plot([R_rumpf_oben_Holm, HLW.b/2],[(kabinenversatz-HLW.l_i_R*0.15),...
    (kabinenversatz-tan(HLW.phi_VK)*(HLW.s_A))-HLW.l_a*0.15],'g');
% 65% Holm
R_rumpf_unten_Holm = ((kabinenversatz-HLW.l_i_R*0.65)/specs.coanlaenge)*specs.R_rumpf;
% Rumpfteil unten
p_HLW(14) = plot([0, R_rumpf_unten_Holm],[(kabinenversatz-HLW.l_i_R*0.65), (kabinenversatz-HLW.l_i_R*0.65)],'g');
% Außenteil unten
p_HLW(15) = plot([R_rumpf_unten_Holm, HLW.b/2],[(kabinenversatz-HLW.l_i_R*0.65),...
    (kabinenversatz-tan(HLW.phi_VK)*(HLW.s_A))-HLW.l_a*0.65],'g');


legend(p_HLW([11, 10, 12]),{'Neutralpunkt Höhenleitwerk', 'phi 25 HLW','Holme'},'Location','northeast','FontSize',18);
title('Veranschaulichung Höhenleitwerk mit Neutralpunkt','FontSize',20);
xlabel('Meter','FontSize',16)
ylabel('Meter','FontSize',16)
hold off


%% Visualisierung Seitenleitwerk

figure(2)

hold on
% Innenkante
p_SLW(1) = plot([0 SLW.l_i],[SLW.s_R SLW.s_R],'b');
% % VK
p_SLW(2) = plot([0 tan(SLW.phi_VK)*(SLW.s_A)],[SLW.s_R, (SLW.s_A + SLW.s_R)],'blue');
% % HK
p_SLW(3) = plot([SLW.l_i (SLW.l_i+tan(SLW.phi_HK_lokal)*(SLW.s_A))],[SLW.s_R (SLW.s_A+SLW.s_R)],'blue');
 % Aussenkante
p_SLW(4) = plot([tan(SLW.phi_VK)*(SLW.s_A), SLW.l_a+tan(SLW.phi_VK)*(SLW.s_A)],[(SLW.s_A+SLW.s_R) (SLW.s_A+SLW.s_R)],'blue');

% Plot Neutralpunkt
p_SLW(5) = plot(SLW.x_NP_ges, SLW.y_SP_ges-SLW.y_SP_R, '*red');
% Winkel phi25
p_SLW(6) = plot([SLW.l_i*0.25, tan(SLW.phi_VK)*(SLW.s_A)+(0.25*SLW.l_a)],[SLW.s_R (SLW.s_A+SLW.s_R)],'--k');
p_SLW(10) = plot([SLW.l_i*0.25 SLW.l_i*0.25],[0 SLW.s_R],'--k');

% Rumpfanteil
 % Vorderkante
p_SLW(7) = plot([0 0],[0 SLW.s_R],'blue');
 % Innenkante
p_SLW(8) = plot([0 SLW.l_i],[0 0],'blue');
 % hinterkante Rumpfteil
p_SLW(9) = plot([SLW.l_i SLW.l_i],[0 SLW.s_R],'blue');

% Holme
% 15% Holm 
% rumpfanteil
p_SLW(10) = plot([(SLW.l_i*0.15), (SLW.l_i*0.15)],[0 SLW.s_R],'green');
% außenanteil
p_SLW(11) = plot([(SLW.l_i*0.15), tan(SLW.phi_VK)*(SLW.s_A)+(0.15*SLW.l_a)],[SLW.s_R (SLW.s_A+SLW.s_R)],'green');
% 65% Holm
p_SLW(10) = plot([(SLW.l_i*0.65), (SLW.l_i*0.65)],[0 SLW.s_R],'green');
% außenanteil
p_SLW(11) = plot([(SLW.l_i*0.65), tan(SLW.phi_VK)*(SLW.s_A)+(0.65*SLW.l_a)],[SLW.s_R (SLW.s_A+SLW.s_R)],'green');

 legend(p_SLW([5:6 11]),{'Neutralpunkt Seitenleitwerk', 'phi 25 SLW','Holme'},'Location','northwest','FontSize',18);
title('Veranschaulichung Seitenleitwerk mit Neutralpunkt','FontSize',20);
grid on
xlabel('Meter','FontSize',16)
ylabel('Meter','FontSize',16)
axis equal
hold off





