clear
close
clc

G_to = 200000;
l_HFW = 90;
l_BFW = 20;
delta_z = 10;
l_BFW_max = 25;
l_BFW_min = 18;
d_reifen = 17;
l_HFW_min = 20;
load Projekt_specs.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bugfahrwerk
%min stat. Belastung
F_BFW_min = (G_to * ((l_BFW+l_HFW)-l_BFW_max))/(l_BFW + l_HFW);

%Max stat. Belastung
F_BFW_max = (G_to * ((l_BFW + l_HFW)-l_BFW_min))/(l_BFW + l_HFW);

%Max dyn. Belastung
F_BFW_dyn = F_BFW_max + ((10*(delta_z + 0.5 * d_reifen)*G_to) / (32.2 * (l_BFW + l_HFW)));

% Muss zwischen 6% und 20% des Abfluggewichts liegen

if ((F_BFW_min >= 0.06*G_to) && (F_BFW_max <= 0.2*G_to))
    disp("Within Limits 1 ohne 7%");
    % Wenn Lasten standhalten dann nochmal zusätzlicher Faktor von 7% und
    % nochmal vergleichen
    F_BFW_min = F_BFW_min + F_BFW_min * 0.07;
    F_BFW_max = F_BFW_max + F_BFW_max * 0.07;
    if ((F_BFW_min >= 0.06*G_to) && (F_BFW_max <= 0.2*G_to))
        disp("Within Limits 2 mit 7%");
    else
        disp("Außerhalb Limits 2 mit 7%");
    end
else
    disp("Außerhalb Limits 1 ohne 7%");
end
%weitere 25% können weggelassen werden weil Roskam/Torenbeek gemacht wurde
%% Hauptfahrwerk
n_FWB = 2; % Oder 3?

F_HFW_max = (G_to *((l_BFW + l_HFW) - l_HFW_min)) / (n_FWB * (l_BFW + l_HFW));


%% Entwurfslasten der Reifen    % PRO REIFEN !!
S_Reserve = 0.25;
S_FW = 1 + ((0.07 + S_Reserve)/100);

n_reifen_BFW = 2;
n_reifen_HFW = 12;

F_reifen_BFW_max = (F_BFW_max/n_reifen_BFW) * S_FW;
F_reifen_BFW_dyn = (F_BFW_dyn/n_reifen_BFW) * S_FW;

F_reifen_HFW_max = (F_HFW_max/n_reifen_HFW) * S_FW;

%% Reifengeschwindigkeiten max

v_s_LDG = 80;
v_s_TO = 100;

v_reifen_max_LDG = 1.3 * v_s_LDG;
v_reifen_max_TO = 1.2 * v_s_TO;

if (v_reifen_max_LDG > v_reifen_max_TO)
    v_reifen_krit = v_reifen_max_LDG;
else
    v_reifen_krit = v_reifen_max_TO;
end

%% Reifendruck 
% Bugfahrwerk High Pressure ab 160PSI/ Aber auch 220PSI
% Hauptfahrwerk Low-Pressure


%% BEstimmung LCN 

S_b =118; %Wheel Base von 1. bis 3. Rad
S_t = 20; % Spurbreite / Track width
tot_load = 9000; % Total Load on one undercarriage assembly
L = 40; %Relative Steifigkeit in inch

% Reifendatenbank
tires = readtable('Reifendaten_angepasst_Flaeche.xlsx');

% ESWL bestimmen
% Erster Schritt Sb/L
T1 = readmatrix('ReductionFactors1.csv');
T3 = readmatrix("LCN_new.csv");
T2 = readmatrix('reductionfactors2.csv');

for i=1:length(T2(:,1))
    for ii=1:length(T2(i,:))
        if isnan(T2(i,ii))
            T2(i,ii) = 0;
        end
    end

end


Sb_L = S_b / L; % in inch ?
St_L = S_t / L; % in Inch ?

rounding = (round(St_L*2)/2);

% Finde passende Spalte
[~,j]=find(T1(1,:)==(St_L*100));

plot(T1(:,j),T1(:,j+1));
hold on
z(:,1) = T1(:,j);
z(:,2) = T1(:,j+1);
z(1,1) = 0;
z(1,2) = 0;
z(2,1) = 0;
z(2,2) = 0;



%%PLotten zum Prüfen
punkt_x = [Sb_L,Sb_L];
punkt_y = [0,50];
punkt = [Sb_L Sb_L;0 50];
plot(punkt_x, punkt_y);


% interx Schnittpunkt finden

P = InterX(z.',punkt);


% Zweiter Schritt Reduction Factor


%Anzahl Räder an  einer Assembly
n_rad = 6;

figure(3);
for ii=1:2:19
loglog(T3(:,ii),T3(:,ii+1));
legend("25","30","35","40","50","60","70","80","90","100")
hold on

end


for i=1:length(tires.A_inch_2_)
    %z2 = zeros(1,1);
    A(i) = (tires.A_inch_2_(i)*n_rad)/(L^2); %Kontaktfläche aus Excel


rounding = (round(A(i)*5)/5);
if rounding > 100
    break
end

[i2,j2]=find(T2(1,:)==rounding);

%figure(2);
%plot(T2(:,j2),T2(:,j2+1));
%hold on
T2(:,(1:46)) = 0;
z2(:,1) = T2(:,j2);
z2(:,2) = T2(:,j2+1);
z2(1,1) = 0;
z2(1,2) = 0;
z2(2,1) = 0;
z2(2,2) = 0;


punkt_y= [P(2),P(2)];
punkt_x = [0,6];
punkt2 = [0 60;P(2) P(2)];
%plot(punkt_x,punkt_y)

P2 = InterX(z2.',punkt2);

reduc_fac(i) = P2(1);
ewsl(i) = (tot_load/reduc_fac(i))/1000;

% 3. Schritt ESWL bestimmen
% ESWL und Reifendruck nutzen um LCN zu ermitteln
loglog(tires.RatedInflation_PSI_(i),ewsl(i),'bx')

%legend(tires.RatedInflation_PSI_(i))
end
Legend= linspace(1,length(ewsl),length(ewsl));
legendStrings = "N = " + string(Legend);
%legend("25","30","35","40","50","60","70","80","90","100",legendStrings)


