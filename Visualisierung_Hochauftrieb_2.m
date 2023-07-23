clc
clear
close

load Ergebnisse_Hochauftrieb_1.mat
load Ergebnisse_Hochauftrieb_2.mat






%% Plotting
%Clean Polare
plot(HA1.alphas,HA1.CAs,'blue','LineWidth',1.5)
hold on

title("Aufgelöste Flügelpolare mit Hochauftriebshilfen","FontSize",15)
ylabel("Auftriebsbeiwert des Flügels C_{A} in [-]","FontWeight","bold")
xlabel("Anstellwinkel \alpha_{F} in °","FontWeight","bold")

%Landing Polare
alphas = -15:0.01:HA2.alpha_F_max_VFFK_deg-HA1.delta_alpha_CA_F_max; % normal Plotten bis alphamax - delta alpha

plot(alphas,HA2.CA_sl,'green','LineWidth',1.5)


%Takeoff Polare
alphas = -12:0.01:(HA2.alpha_F_max_VFFK_deg_TO-HA1.delta_alpha_CA_F_max); % normal Plotten bis alphamax - delta alpha

plot(alphas,HA2.CA_st,'red','LineWidth',1.5)

% Kritische Punkte -> CA_Max stimmt nicht perfekt mit ende der geraden
% überein ?!
grid on
plot([HA2.alpha_F_max_VFFK_deg_TO - HA1.delta_alpha_CA_F_max],[0],'redx','LineWidth',1.5)
plot([HA2.alpha_F_max_VFFK_deg - HA1.delta_alpha_CA_F_max],[0],'greenx','LineWidth',1.5)

plot([0],[HA2.CA_F_max_VFFK_TO],'redx','LineWidth',1.5)
plot([0],[HA2.CA_F_max_VFFK],'greenx','LineWidth',1.5)


plot(HA1.alpha_Ca_F_max, 0, 'xblue','LineWidth',1.5)
%CA MAX
plot(0, HA1.CA_F_max,"xblue",'LineWidth',1.5)

%% Abfall der Polaren Plotten

% Parameter der quadratischen Funktion
a = 0.008; % Koeffizient von x^2
h = HA2.alpha_F_max_VFFK_deg; % x-Koordinate des Maximums
k = HA2.CA_F_max_VFFK; % y-Koordinate des Maximums

% Bereich der x-Achse
x = linspace(h-6, h+7, 10); % Hier können Sie den Bereich anpassen

% Quadratische Funktion berechnen
y = -a * (x - h).^2 + k;

% Plot erstellen
plot(x, y,"green--",'LineWidth',1.5)
plot(h, k, 'gx','LineWidth',1.5)


% Parameter der quadratischen Funktion
a = 0.0072; % Koeffizient von x^2
h = HA2.alpha_F_max_VFFK_deg_TO; % x-Koordinate des Maximums
k = HA2.CA_F_max_VFFK_TO; % y-Koordinate des Maximums

% Bereich der x-Achse
x = linspace(h-7, h+7, 10); % Hier können Sie den Bereich anpassen

% Quadratische Funktion berechnen
y = -a * (x - h).^2 + k;

% Plot erstellen
plot(x, y,"red--",'LineWidth',1.5)
plot(h, k, 'rx','LineWidth',1.5)

% Parameter der quadratischen Funktion
a = 0.0058; % Koeffizient von x^2
h = HA1.alpha_Ca_F_max; % x-Koordinate des Maximums
k = HA1.CA_F_max; % y-Koordinate des Maximums

% Bereich der x-Achse
x = linspace(h-7, h+7, 10); % Hier können Sie den Bereich anpassen

% Quadratische Funktion berechnen
y = -a * (x - h).^2 + k;

% Plot erstellen
plot(x, y,"blue--",'LineWidth',1.5)
plot(h, k, 'bx','LineWidth',1.5)
% Achsen
P1 = [0 0];
P2 = [-1 3];
plot(P1,P2,'black','LineWidth',1.5)
P3 = [-15 35];
P4 = [0 0];
plot(P3,P4,'black','LineWidth',1.5)

ylim([-1, 3])

legend("Clean Konfiguration mit 0° Klappenausschlag","Landing mit 45° Klappenausschlag","Takeoff mit 20° Klappenausschlag",'','','','','Location', 'southeast')
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Widerstandspolaren















