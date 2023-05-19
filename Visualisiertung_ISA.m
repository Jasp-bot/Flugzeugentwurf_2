clear all                       %l�scht alle vorhandenen Eintr�ge im Workspace
close all                       %schlie�t alle ge�ffneten Figure()-Fenster
clc                             %l�scht alle Eintr�ge im Command-Window


%% Eingabe Parameter für Funktion ISA_Data

% Startwert der Höhe
H_Start = 1;

% Schrittweite zwischen H_Start und H_Ende
schrittweite = 1;

% max Höhe
H_Ende = 32000;

% delta der ISA
ISAPlus = 15;

%% Berechnung der Werte

ISA = Berechnung_ISA(H_Start,schrittweite ,H_Ende ,ISAPlus);

%atmosData(:,1) = Hoehe [m]
%atmosData(:,2) = Temperatur [K]
%atmosData(:,3) = Druck [N/m^2]
%atmosData(:,4) = Schallgeschwindigkeit [m/s]
%atmosData(:,5) = Dichte [kg/m^3]
%atmosData(:,6) = relativer Dichteverlauf
%atmosData(:,7) = dynamische Viskosität
%atmosData(:,8) = kinematische Viskosität


%% Plot der Werte / Visualisierung

%Muss noch erweitert werden
% figure()        % �ffnet neues MatLab-fenster in dem geplotet werden soll




% plot erm�glicht das zeichen von Graphen in einem Fenster, Achtung x und y m�ssen gleich lang sein
hold on
t=tiledlayout(1,1);


ax1 = axes(t);
hold on         % erm�glicht ploten von mehreren Graphen in ein Fenster
grid on         % blendet Gitternetzlinien in der Grafik ein
plot(ax1, ISA.T, ISA.H/1000, 'r', ...
    ISA.a, ISA.H/1000, '-.k')
ax1.XAxisLocation = 'Top';
%ax1.YAxisLocation = 'Left';
xlabel('Temperatur T in K   ||   Schallgeschwindigkeit a in m/s');
hold on
legend('T(H)', 'a(H)')



ax2 =axes(t);
% 
 plot(ax2, ISA.p/100000, ISA.H/1000, 'k',...
     ISA.rho*1, ISA.H/1000, 'g', ...
     ISA.rel_rho*1, ISA.H/1000, 'b')
ax2.XAxisLocation = 'Bottom';
ax2.YAxisLocation = 'Left';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
xlabel('Druck in Pa*1/(10^5) || Dichte rho in kg/(m^3) || relative Dichte Theta ')




% plot(atmosData(:,2), atmosData(:,1), 'r')
% plot(atmosData(:,3)/500, atmosData(:,1)/1000, '-.k')    % plot Druck
% plot(atmosData(:,4), atmosData(:,1)/1000, 'k')          % plot Schallgeschwindigkeit
% plot(atmosData(:,5)*100, atmosData(:,1)/1000, 'g')          % plot Dichte
% plot(atmosData(:,6)*100, atmosData(:,1)/1000, 'b')          % plot relative Dichte theta
% %plot(atmosData(:,7)*10^2, atmosData(:,1)/1000)               % plot dynamische Voskosität my
% %plot(atmosData(:,8)*10^5, atmosData(:,1)/1000)               % plot Kinematische Viskosität ny

% % Vektor f�r Lininen er Atmosp�hrenschichten
hold on
X(1,1) = 0;
X(2,1) = 1.5;
Y(1,1) = ISA.H_Tropo_ISA/1000;
Y(2,1) = ISA.H_Tropo_ISA/1000;
Z(1,1) = ISA.H_Str_u/1000;
Z(2,1) = ISA.H_Str_u/1000;
% 
% %ax3 =axes(t);
% %plot(ax3, X, Y, '-.b', X, Z, '-.b')
plot(X, Y, '-.b')                                                         %ploten der Linie f�r Tropopause nach ISA
% %text erm�glicht das Beschriften von einzelnen Graphen in einem Plot -> siehe Documentation f�r weitere Informationen
text(X(2,1)-0.4, Y(2,1)+0.90,  'Tropopause nach ISA', 'FontSize', 10)  
% axes hide
% 
% %ax3.Color = 'none';
% 
plot(X, Z, '-.b')
text(X(2,1)-0.4, Z(2,1)+0.90,  'zunahme Ozonschicht', 'FontSize', 10)    %ploten der Linie f�r Beginn der Zunahme der Ozonschicht nach ISA
% 
% 
% %Beschriftungen
% 
xlim([0 1.5])                      %erm�glicht die Vorgabe eines festen Bereiches der Dargestellt werden soll -> hier wird x-Achse nur von 160 bis 360 dargestell
ylim([0 35])                         %erm�glicht die Vorgabe eines festen Bereiches der Dargestellt werden soll -> hier wird y-Achse nur von 0 bis 32 dargestellt

legend('p(H)', 'rho(H)', 'theta(H)','','', 'Location','north')

%xlabel('Temperatur T in K   ||   Schallgeschwindigkeit a in m/s')       %Achsenbeschriftung der x-Achse
ylabel('Hoehe in Kilometer')                                             %Achsenbeschriftung der y-Achse
title('Verlauf der Atmosphaerendaten')                                   %Titel der Abbildung
%legend('T(H)', 'p(H)', 'a(H)', 'rho(H)', 'theta(H)')                        %, 'my(h)', 'ny(H)')                                                  %Erstellt eine Legende anhand der vorgegeben Daten, Achtung: muss immer soviele Eintr�ge enthalten wie Graphen im Plot zu sehen sind



















