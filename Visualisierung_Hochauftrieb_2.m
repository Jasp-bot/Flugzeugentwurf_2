clc
clear all
close all

load Ergebnisse_Hochauftrieb_1.mat
load Ergebnisse_Hochauftrieb_2.mat
load Ergebnisse_Start_Landeanforderungen.mat
load Ergebnisse_Widerstand_FE2.mat
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat




%% Plotting
%Clean Polare
plot(HA1.alphas,HA1.CAs,'blue','LineWidth',1.5)
hold on

title("Aufgelöste Flügelpolare mit Hochauftriebshilfen","FontSize",15)
ylabel("Auftriebsbeiwert des Flügels C_{A} in [-]","FontWeight","bold")
xlabel("Anstellwinkel \alpha_{F} in °","FontWeight","bold")

%Landing Polare
%alphas = -15:0.01:HA2.alpha_F_max_VFFK_deg - HA1.delta_alpha_CA_F_max; % normal Plotten bis alphamax - delta alpha

plot(HA2.alphas_sl,HA2.CA_sl,'green','LineWidth',1.5)


%Takeoff Polare
%alphas = -12:0.01:(HA2.alpha_F_max_VFFK_deg_TO - HA1.delta_alpha_CA_F_max); % normal Plotten bis alphamax - delta alpha

plot(HA2.alphas_st,HA2.CA_st,'red','LineWidth',1.5)

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
x = linspace(h-4, h+7, 10); % Hier können Sie den Bereich anpassen

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
x = linspace(h-4, h+7, 10); % Hier können Sie den Bereich anpassen

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

ylim([-0.2, 3])



%% Check ob genung CA TO / LDG

plot([0,25],[startschub.c_A_max_thrust_match startschub.c_A_max_thrust_match],'red--')

% Check ob lande CA erreicht
plot([0,25],[landeanvorderung.c_A_max_LDG landeanvorderung.c_A_max_LDG],'green--')


legend("Clean Konfiguration mit 0° Klappenausschlag","Landing mit 45° Klappenausschlag","Takeoff mit 20° Klappenausschlag",'','','','','','','','','','','','','','','Notwendiger Auftriebsbeiwert aus der Startanforderung','Notwendiger Auftriebsbeiwert aus der Landeanforderung','Location', 'southeast')
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Widerstandspolaren

%finde das Cruise CA im Array
targetValue = Ergebnisse_stat_Flaechenbelastung.C_A_CR;
absoluteDifferences = abs(Ergebnisse_Widerstand_FE2.c_A_ges - targetValue);
[minDifference, index] = min(absoluteDifferences);

figure(2)
hold on

plot(Ergebnisse_Widerstand_FE2.x_vector_sum(9,1:index),Ergebnisse_Widerstand_FE2.c_A_ges(1:index),'blue','LineWidth',1.5)
hold on

plot(HA2.TO_CW,HA2.c_A_F_TO,'red','LineWidth',1.5)
plot(HA2.TO_CW_FW,HA2.c_A_F_TO_FW,'red--','LineWidth',1.5)

plot(HA2.LDG_CW,HA2.c_A_F_LDG,'green','LineWidth',1.5)
plot(HA2.LDG_CW_FW,HA2.c_A_F_LDG_FW,'green--','LineWidth',1.5)



title("Lillienthalpolaren")
ylabel("Auftriebsbeiwert des Flügels C_{A} in [-]","FontWeight","bold")
xlabel("Widerstandsbeiwert des Flügels C_{W} in [-]","FontWeight","bold","FontWeight","bold")


grid on
legend("Clean Konfiguration im Cruise Zustand","Takeoff mit 20° Klappenausschlag","Takeoff mit 20° Klappenausschlag und Fahrwerk","Landing mit 45° Klappenausschlag","Landing mit 45° Klappenausschlag und Fahrwerk",'Location', 'southeast')
ylim([0, 3])
hold off






%% Reziproke Gleitzahl polaren E, nicht epsillon

figure(3)
grid on
hold on

plot(Ergebnisse_Widerstand_FE2.c_A_ges(1:index),(Ergebnisse_Widerstand_FE2.c_A_ges(1:index)./Ergebnisse_Widerstand_FE2.x_vector_sum(9,1:index)),'blue','LineWidth',1.5)

plot(HA2.c_A_F_TO,(HA2.c_A_F_TO./HA2.TO_CW),'red','LineWidth',1.5)
plot(HA2.c_A_F_TO_FW,(HA2.c_A_F_TO_FW./HA2.TO_CW_FW),'red--','LineWidth',1.5)

plot(HA2.c_A_F_LDG,(HA2.c_A_F_LDG./HA2.LDG_CW),'green','LineWidth',1.5)
plot(HA2.c_A_F_LDG_FW,(HA2.c_A_F_LDG_FW./HA2.LDG_CW_FW),'green--','LineWidth',1.5)


ylim([0, 20])
xlim([0,3])


title("Reziproke Gleitzahl")
ylabel("Gleitzahl E = C_{A}/C_W in [-]","FontWeight","bold")
xlabel("Auftriebsbeiwert des Flügels C_{A} in [-]","FontWeight","bold")

% PLotten der Gleitzahlen
plot(startschub.c_A_max_thrust_match,startschub.Eta_To_inv(3),'redo','LineWidth',2)

% Check ob lande CA erreicht
plot(landeanvorderung.c_A_max_LDG, (1/landeanvorderung.Eta_LDG),'greeno','LineWidth',2)



legend("Clean Konfiguration im Cruise Zustand","Takeoff mit 20° Klappenausschlag","Takeoff mit 20° Klappenausschlag und Fahrwerk","Landing mit 45° Klappenausschlag","Landing mit 45° Klappenausschlag und Fahrwerk","Gleitzahlvorgabe aus der Startanforderung","Gleitzahlvorgabe aus der Landeanforderung",'Location', 'southeast')

