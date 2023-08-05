% Visualisierung Start und Landeanforderung
clc
clear all
close all


%% Ausfuerung des Programms
% load Ergebnisse_Endwerte_Iteration_V1.mat
%  Berechnung_PS6_Startschub_Landeanforderung(Endwerte_Iteration,1);

%% Laden der Dateien
load Ergebnisse_Start_Landeanforderungen.mat

%% Visualisierung


hold on
plot(startschub.s1 ,startschub.S0_GTo_To(1,:),'b', ...
     startschub.s1 ,startschub.S0_GTo_To(2,:),'--b', ... 
     startschub.s1 ,startschub.S0_GTo_To(3,:),'-.b')  %,s1 ,S0_GTo_To(4,:),'b',s1 ,S0_GTo_To(5,:),'-.b',s1 ,S0_GTo_To(6,:),'--b'); % c_A_max

plot(startschub.s1, startschub.S0_GTo_CL(1,:),'-.k', ...
     startschub.s1, startschub.S0_GTo_CL(2,:),'--k', ...
     startschub.s1, startschub.S0_GTo_CL(3,:),'k')  %s1, S0_GTo_CL(4,:),'--k',s1, S0_GTo_CL(5,:),'-.k'); % Eta
grid on
xlim([1000 2600])
ylim([0 0.9])

% plot fuer Schnittpunkte der Graphen mit InterX
intersections_1 = InterX([startschub.s1; startschub.S0_GTo_To(1,:)],[startschub.s1; startschub.S0_GTo_CL(3,:)]);
intersections_2 = InterX([startschub.s1; startschub.S0_GTo_To(2,:)],[startschub.s1; startschub.S0_GTo_CL(2,:)]);
intersections_3 = InterX([startschub.s1; startschub.S0_GTo_To(3,:)],[startschub.s1; startschub.S0_GTo_CL(1,:)]);

plot(intersections_1(1,:), intersections_1(2,:),'or',...
    intersections_2(1,:), intersections_2(2,:),'or',...
    intersections_3(1,:), intersections_3(2,:),'or')


%plot(s1, S0_GTo_CL(1,:),s1, S0_GTo_CL(2,:),s1, S0_GTo_CL(3,:),s1, S0_GTo_CL(4,:),s1, S0_GTo_CL(5,:))

% plot(s1, S0_GTo_CR(1,1), s1, S0_GTo_CR(2,1),s1, S0_GTo_CR(3,1),s1, S0_GTo_CR(4,1))

% X(1,1) = 0;
% X(2,1) = 2800;
% Y(1,1) = S0_GTo_CR(1,1);
% Y(2,1) = S0_GTo_CR(1,1);
% Y2(1,1) = S0_GTo_CR(2,1);
% Y2(2,1) = S0_GTo_CR(2,1);
% Y3(1,1) = S0_GTo_CR(3,1);
% Y3(2,1) = S0_GTo_CR(3,1); 
% Y4(1,1) = S0_GTo_CR(4,1);
% Y4(2,1) = S0_GTo_CR(4,1);
% 
plot([0 2800],[schub_CR.S0_GTo_CR schub_CR.S0_GTo_CR],'--r')%,X,Y2,'--r',X,Y3, X,Y4);
%plot(1656,0.38877,'ored')
% Variabele Legende

C_A_string_1 = sprintf('C_A %.1f',startschub.c_A_max(1,1));
C_A_string_2 = sprintf('C_A %.1f',startschub.c_A_max(2,1));
C_A_string_3 = sprintf('C_A %.1f',startschub.c_A_max(3,1));

ETA_string_1 = sprintf('Eps %.1f', startschub.Eta_To_inv(1,1));
ETA_string_2 = sprintf('Eps %.1f', startschub.Eta_To_inv(2,1));
ETA_string_3 = sprintf('Eps %.1f', startschub.Eta_To_inv(3,1));

legend(C_A_string_1, C_A_string_2, C_A_string_3,...
    ETA_string_1 , ETA_string_2, ETA_string_3,...
    '','','', 'Schub-Gewichtsverhältnis CR','FontSize',18); %, 'Designpunkt')
xlabel('Startstrecke in m','FontSize',16);
ylabel('Schub-Gewichtsverhältnis','FontSize',16);
title('Thrustmatching','FontSize',20);


