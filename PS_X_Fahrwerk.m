clear
close  
clc

G_to = 2.6258e+06;


delta_z = 7.3;
l_BFW_max = 36.495;
l_BFW_min = 33.194;
%d_reifen = 17;
l_HFW_min = 2.793;
l_HFW = l_HFW_min;
l_BFW = l_BFW_max;

load Projekt_specs.mat

% Reifendatenbank
tires = readtable('Reifendaten_angepasst_Flaeche.xlsx');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Reifengeschwindigkeiten

v_s_LDG = 80; %-> Specs
v_s_TO = 100; %-> Specs

v_reifen_max_LDG = 1.3 * v_s_LDG;
v_reifen_max_TO = 1.2 * v_s_TO;

%Ist v_takeoff/Landing aus Spezifikation größer als diese Rechnung dann
%diese verwenden

if (v_reifen_max_LDG > v_reifen_max_TO)
    v_reifen_krit = v_reifen_max_LDG;
else
    v_reifen_krit = v_reifen_max_TO;
end

v_reifen_krit = 206.9667;%v_reifen_krit * 1.151078;

%% 1. Hauptfahrwerk Belastung ermitteln
S_Reserve = 25;
S_FW = 1 + ((7 + S_Reserve)/100);

n_FWB = 2; % Es handelt sich um Hauptfahrwerksbeine
n_reifen_HFW = 6; % Anzahl Reifen insgesamt

F_HFW_max = (G_to *((l_BFW + l_HFW) - l_HFW_min)) / (n_FWB * (l_BFW + l_HFW));

F_reifen_HFW_max = (F_HFW_max/n_reifen_HFW) * S_FW;
F_reifen_HFW_max = convforce(F_reifen_HFW_max,"N","lbf");


idx = (tires.RatedLoad_Lbs_ >= F_reifen_HFW_max & tires.RatedInflation_PSI_<=220 & tires.RatedSpeed_MPH_ >= v_reifen_krit);
tires_edit_HFW = tires(idx,:);





%% 2. LCN bestimmen -> beste Reifen wählen

%% BEstimmung LCN 

S_b = 50; %Wheel Base von 1. bis 3. Rad 
S_t = 30; % Spurbreite / Track width        -> Variabel für Iteration


tot_load = F_HFW_max; % Total Load on one undercarriage assembly
L = 40; %Relative Steifigkeit in inch
n_rad = 6; % Sechs Räder pro Undercarriage ASSEMBLY

% ESWL bestimmen

% Erster Schritt Sb/L
T1 = readmatrix('ReductionFactors1.csv');

Sb_L = S_b / L; % in inch ?
St_L = S_t / L; % in Inch ?
Y = roundtowardvec(St_L,[25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 110 120 130 140 150 160 170]);


% Finde passende Spalte
figure(1)
[~,j]=find(T1(1,:)==(Y));
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
 hold off

% interx Schnittpunkt finden

P = InterX(z.',punkt);



% Zweiter Schritt Reduction Factor und Contact Area
T2 = readmatrix('reductionfactors2.csv');
T2(:,(1:46)) = 0;

T3 = readmatrix("LCN_new.csv");

figure(3);
 for ii=1:2:19
 loglog(T3(:,ii),T3(:,ii+1));
 hold on
 end
 legend("25","30","35","40","50","60","70","80","90","100")


for i=1:length(tires_edit_HFW.A_inch_2_)
    A(i) = ((tires_edit_HFW.A_inch_2_(i)*n_rad)/L^2);% Hier eigentlich L^2 aber dann ergebnisse unrealistisch ; %Kontaktfläche aus Excel

A(i) = A(i)*100;

X=roundtowardvec(A(i),[5 8 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 110 120 130 140 150 160 170 180 190 200]);

% rounding = (round(A(i)*200)/200)*100;
% rounding = round(rounding);
% 
% % Wenn Wert lost zwischen zwei möglichen 
% if rounding <=6 
%    rounding =5;
% elseif (rounding >=6 && rounding <=7)
%     rounding = 8;
% elseif rounding ==9
%     rounding =10;
% end



[i2,j2]=find(T2(1,:)==X);
if isempty(i2)
    ewsl(i) = 0;
    continue
end
%figure(7);
%plot(T2(:,j2),T2(:,j2+1));
%hold on

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
if isempty(P2)
    ewsl(i) = 0;
    continue
end
reduc_fac(i) = P2(1);
ewsl(i) = (tot_load/reduc_fac(i))/10000;

% 3. Schritt ESWL bestimmen
% ESWL und Reifendruck nutzen um LCN zu ermitteln

%figure(3)
loglog(tires_edit_HFW.RatedInflation_PSI_(i),ewsl(i),'bx')

in(i) = inpolygon(tires_edit_HFW.RatedInflation_PSI_(i),ewsl(i),[71 88 218 238],[62 77 25 35]);
clear figure(2)
end
%% Check ob Reifen im LCN Limit liegt -> Über Polygons
big = 0;
big_idx = 0;
%for u=1:length(in)
%    if ewsl(u) > big && in(u) ==1
%        big = ewsl(u);
%        big_idx = u;
%    end
%end

%% 2. Bugfahrwerk ermitteln
%% Belastungen
%min stat. Belastung -> hintere Lage
F_BFW_min = (G_to * ((l_BFW+l_HFW)-l_BFW_max))/(l_BFW + l_HFW);

%Max stat. Belastung -> Vordere Lage
F_BFW_max = (G_to * ((l_BFW + l_HFW)-l_BFW_min))/(l_BFW + l_HFW);

%Max dyn. Bremslast -> Abhängig von durchmesser Reifen
%Outisde Diameter MAX verwenden?
d_reifen = tires_edit_HFW.OutsideDiameterMax(big_idx);
F_BFW_dyn(i) = F_BFW_max + ((10*(delta_z + 0.5 * d_reifen)*G_to) / (32.2 * (l_BFW + l_HFW)));


% Statische Belastung muss zwischen 6% und 20% des Abfluggewichts liegen
% NIEMALS Über/unterschreiten -> OPtimaler ist 8-15

if ((F_BFW_min >= 0.06*G_to) && (F_BFW_max <= 0.20*G_to))
    disp("Within Limits 1");
    if ((F_BFW_min >= 0.08*G_to) && (F_BFW_max <= 0.15*G_to))
        disp("Within Limits 2");
    else
        disp("Außerhalb Limits 2");
    end
else
    disp("Außerhalb Limits 1");
end

% Statische Belastung muss mit 7% Reserve erhöht werden 
% Auch 25% Reserve falls GTO sich erhöht, Aber A.Gobbin -> Wenn Roskam
% fertig reicht nur 7% Betrachtung !!Kristof fragen!!

%F_BFW_min = F_BFW_min + F_BFW_min * 0.07; -> Minimal nicht wichtig?

%F_BFW_max = F_BFW_max + F_BFW_max * 0.07;

%F_BFW_dyn = F_BFW_dyn + F_BFW_dyn .* 0.07; 

%weitere 25% können weggelassen werden weil Roskam/Torenbeek gemacht
%wurde??
S_Reserve = 0.25;
S_FW = 1 + ((0.07 + S_Reserve)/100);

n_reifen_BFW = 2;

F_reifen_BFW_max = (F_BFW_max/n_reifen_BFW) * S_FW;
F_reifen_BFW_dyn = (F_BFW_dyn./n_reifen_BFW) .* S_FW;
F_reifen_BFW_dyn_reduc = F_reifen_BFW_dyn./1.5;

if F_reifen_BFW_dyn_reduc > F_reifen_BFW_max
    reifencrit = F_reifen_BFW_dyn_reduc;
else
    reifencrit = F_reifen_BFW_max;
end

[~,j3]=find(tires.RatedInflation_PSI_>= 200 & tires.RatedLoad_Lbs_>=reifencrit & tires.RatedSpeed_MPH_>=v_reifen_krit);

tires.PartNumber(j3)

%% Auswahl Bugfahrwerksreifen






%% Reifendruck 
% Bugfahrwerk High Pressure ab 160PSI/ Aber auch 220PSI
% Hauptfahrwerk Low-Pressure





