clc
clear
close

load Ergebnisse_Hochauftrieb_1.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(HA1.alphas,HA1.CAs,'blue','LineWidth',1.5)

hold on
title("Aufgelöste Flügelpolare ohne Hochauftriebshilfen","FontSize",15)
ylabel("Auftriebsbeiwert des Flügels C_{A} in [-]","FontWeight","bold")
xlabel("Anstellwinkel \alpha_{F} in °","FontWeight","bold")
P1 = [0 0];
P2 = [-1 2];
plot(P1,P2,'black')
P3 = [-10 25];
P4 = [0 0];
plot(P3,P4,'black')

% Kritische Points
%Alpha MAX
plot(HA1.alpha_Ca_F_max, 0, 'xred')

%CA MAX
plot([0 20],[HA1.CA_F_max HA1.CA_F_max], 'red--')
plot(0, HA1.CA_F_max,"xred")

%Alpha 0
plot(HA1.alpha_MAC_0, 0,'xred')

plot([HA1.alpha_Ca_F_max  HA1.alpha_Ca_F_max-HA1.delta_alpha_CA_F_max],[HA1.CA_F_max HA1.CA_F_max],'xblue')

grid on

%% Auftriebsabfall Plotten

% Parameter der quadratischen Funktion
a = 0.0058; % Koeffizient von x^2
h = HA1.alpha_Ca_F_max; % x-Koordinate des Maximums
k = HA1.CA_F_max; % y-Koordinate des Maximums

% Bereich der x-Achse
x = linspace(h-5.5, h+7, 10); % Hier können Sie den Bereich anpassen

% Quadratische Funktion berechnen
y = -a * (x - h).^2 + k;

% Plot erstellen
plot(x, y,"blue--",'LineWidth',1.5)
plot(h, k, 'bx','LineWidth',1.5)
ylim([-0.2, 1.6])
xlim([-7, 20])



