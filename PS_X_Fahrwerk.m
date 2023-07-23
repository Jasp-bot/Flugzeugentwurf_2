clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Das mit Variablen aus CG ersetzen!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Projekt_specs.mat
load Ergebnisse_Auftrieb_Momente.mat
load Ergebnisse_Fluegel_Tank_NP.mat
load Ergebnisse_Auftrieb_Momente.mat
load Ergebnisse_Leitwerke.mat
load Ergebnisse_Start_Landeanforderungen.mat
load Zwischenergebnisse_PS5_Fluegelflaechen.mat
load Ergebnisse_Widerstand_FE2.mat
load Ergebnisse_Massen_FE2.mat
load Schwerpunkt.mat



% Reifendatenbank
tires = readtable('Reifendaten_angepasst_Flaeche.xlsx');

G_to = Ergebnisse_Massen_FE2.M_TO * 9.81; %2.6258e+06;
%G_to = 568793 * 9.81


delta_z = LS.delta_z; 
l_BFW_max = BFWL.l_BFW_max;
l_BFW_min = BFWL.l_BFW_min;

l_HFW_min = BFWL.l_HFW_min;

% Das Nochmal prüfen?!
l_HFW = l_HFW_min;
l_BFW = l_BFW_max;

% Finales HFW --> Durchmesser wählen 
d_reifen_aus_HFW = 52;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Reifengeschwindigkeiten

v_s_LDG = landeanvorderung.v_max_LDG; %-> Specs
v_s_TO = startschub.v_LOF; %-> Specs

v_reifen_max_LDG = 1.3 * v_s_LDG;
v_reifen_max_TO = 1.2 * v_s_TO;

%Ist v_takeoff/Landing aus Spezifikation größer als diese Rechnung dann
%diese verwenden

% if (v_reifen_max_LDG > v_reifen_max_TO)
%     v_reifen_krit = v_reifen_max_LDG;
% else
%     v_reifen_krit = v_reifen_max_TO;
% end

v_reifen_krit = v_reifen_max_LDG * 2.237; %220; % Einfach Festgelegt

%% 1. Hauptfahrwerk Belastung ermitteln

S_Reserve = 25; %Drehschraube Normalerweise 25%

S_FW = 1 + ((S_Reserve + 7)/100);

n_FWB = 2; % Es handelt sich um Hauptfahrwerksbeine
n_reifen_HFW = 8; % Anzahl Reifen insgesamt

F_HFW_max = (G_to *((l_BFW + l_HFW) - l_HFW_min)) / (n_FWB * (l_BFW + l_HFW));
F_HFW_max = F_HFW_max / 9.81;
F_HFW_max_pfund = convmass(F_HFW_max,'kg','lbm');
%F_HFW_max = convforce(F_HFW_max,"N","lbf");

F_reifen_HFW_max = (F_HFW_max/n_reifen_HFW) * S_FW;

F_reifen_HFW_max_pfund = convforce(F_reifen_HFW_max,'N','lbf');

%F_reifen_HFW_max = convmass(F_reifen_HFW_max,'kg','lbm');

%F_reifen_HFW_max = convforce(F_reifen_HFW_max,"N","lbf");


idx = ((tires.RatedLoad_Lbs_ >= F_reifen_HFW_max_pfund) & (tires.RatedInflation_PSI_<=240) & (tires.RatedSpeed_MPH_ >= v_reifen_krit) );
tires_edit_HFW = tires(idx,:);


%% 2. LCN bestimmen -> beste Reifen wählen

%% Bestimmung LCN 
tot_load = F_HFW_max_pfund; % Total Load on one undercarriage assembly
L = 40; %Relative Steifigkeit in inch
n_rad = 4; % Sechs Räder pro Undercarriage ASSEMBLY

% ESWL bestimmen

% Erster Schritt Sb/L
T1 = readmatrix('ReductionFactors1_new_2.csv');

for zz =1: length(tires_edit_HFW.A_inch_2_)
S_b = 2.3 * tires_edit_HFW.OutsideDiameterMin(zz); %Wheel Base von 1. bis 3. Rad 
S_t = 1.2 * tires_edit_HFW.OutsideDiameterMin(zz); % Spurbreite / Track width        -> Variabel für Iteration
      % Faktoren aus PS

Sb_L(zz) = S_b / L;
St_L(zz) = (S_t / L)*100; %-> Damit aus digitizer nutzbar
Y(zz) = roundtowardvec(St_L(zz),[25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 110 120 130 140 150 160 170]);

% Finde passende Spalte
%figure(1)
[~,j]=find(T1(1,:)==(Y(zz)));
%plot(T1(:,j),T1(:,j+1));
%hold on
z(:,1) = T1(:,j);
z(:,2) = T1(:,j+1);
z(1,1) = 0;
z(1,2) = 0;
z(2,1) = 0;
z(2,2) = 0;

%%PLotten zum Prüfen
 punkt_x = [Sb_L(zz),Sb_L(zz)];
 punkt_y = [0,50];
 punkt = [Sb_L(zz) Sb_L(zz);0 50];
 %plot(punkt_x, punkt_y);
 %hold off

% interx Schnittpunkt finden

P = InterX(z.',punkt);
if (isempty(P))
    Points(:,zz) = 0;
else
    Points(:,zz) = P(:,1);
end
end


% 2. Schritt Reduction Factor und Contact Area
T2 = readmatrix('reductionfactors2_new_2.csv');
T2(:,(1:46)) = 0;
T3 = readmatrix("LCN_new.csv");


%% Plot für LCN erstellen
figure(3);
 for ii=1:2:19
 loglog(T3(:,ii),T3(:,ii+1),'LineWidth',1.5);
 hold on
 end
 legend("25","30","35","40","50","60","70","80","90","100",'AutoUpdate','off','Location','northwest')


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 2 im Diagramm
for i=1:length(tires_edit_HFW.A_inch_2_)
    A(i) = (((tires_edit_HFW.A_inch_2_(i))*n_rad)/(L^2));% Hier eigentlich L^2 aber dann ergebnisse unrealistisch ; %Kontaktfläche aus Excel

X=roundtowardvec(A(i),[0.05 0.08 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.10 1.20 1.30 1.40 1.50 1.60 1.70 1.80 1.90 2.00]);
X = X*100;
%if X == 0
%continue
%end

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

punkt_y= [Points(2,i),Points(2,i)];
punkt_x = [0,20];
punkt2 = [0 20;Points(2,i) Points(2,i)];
%plot(punkt_x,punkt_y)
%hold off
P2 = InterX(z2.',punkt2);
if isempty(P2)
    ewsl(i) = 0;
    continue
end
reduc_fac(i) = P2(1);
ewsl(i) = (tot_load/reduc_fac(i))/1000;

% 3. Schritt ESWL bestimmen
% ESWL und Reifendruck nutzen um LCN zu ermitteln

figure(3)
loglog(tires_edit_HFW.RatedInflation_PSI_(i),ewsl(i),'bx','LineWidth',1.5)

if (ewsl(i) <= 0 || isinf(ewsl(i))== 1)

    tires_edit_HFW.RatedInflation_PSI_(i) = 0;
end 


if inpolygon(tires_edit_HFW.RatedInflation_PSI_(i),ewsl(i),[150 160 218 238],[62 77 25 35])
in(i) = 1;
else 
in(i) = 0;
end
end
title("Bestimmung der Reifen über LCG/LCN Methode","FontWeight","bold")
ylabel("Equivalent Single Wheel Load (ESWL) in 1000lb","FontWeight","bold")
xlabel("Druck in PSI","FontWeight","bold")
xlim([50 260])
grid on

%% Check ob Reifen im LCN Limit liegt -> Über Polygons
% big = 0;
% big_idx = 0;
% for u=1:length(in)
%     if (ewsl(u) > big) && (in(u) == 1)
%         big = ewsl(u);
%         big_idx = u;
%     end
% end
% tires_edit_HFW.PartNumber(big_idx)








%% 2. Bugfahrwerk ermitteln
%% Belastungen

%min stat. Belastung -> hintere Lage
F_BFW_min = (G_to * ((l_BFW+l_HFW)-l_BFW_max))/(l_BFW + l_HFW);
F_BFW_min = convforce(F_BFW_min,'N','lbf');

%Max stat. Belastung -> Vordere Lage
F_BFW_max = (G_to * ((l_BFW + l_HFW)-l_BFW_min))/(l_BFW + l_HFW);
F_BFW_max = convforce(F_BFW_max,'N','lbf');

%F_BFW_max = convmass(F_BFW_max,'kg','lbm')

%Max dyn. Bremslast -> Abhängig von durchmesser Reifen
%Outisde Diameter MIN verwenden!
d_reifen = convlength(d_reifen_aus_HFW,'in','m');%tires_edit_HFW.OutsideDiameterMin(big_idx);
F_BFW_dyn = convforce(F_BFW_max,'lbf','N') + ((10 * (delta_z + 0.5 * d_reifen)*G_to/9.81) / (32.2 * (l_BFW + l_HFW)));
F_BFW_dyn = convforce(F_BFW_dyn,'N','lbf');

F_BFW_Dyn_transf = F_BFW_dyn/1.5;


% Auf Bugfahrwerksreifen bezogen
n_BFW_reifen = 2;

F_reifen_BFW_max = S_FW * (F_BFW_max/n_BFW_reifen);

F_reifen_BFW_dyn = S_FW * (F_BFW_dyn/n_BFW_reifen); 



% Statische Belastung muss zwischen 6% und 20% des Abfluggewichts liegen
% NIEMALS Über/unterschreiten -> OPtimaler ist 8-15

if ((convforce(F_BFW_max,'lbf','N')/9.81 <= 0.20 * G_to/9.81) && (convforce(F_BFW_min,'lbf','N')/9.81 <= 0.06 * G_to/9.81))
    disp("Innerhalb Limits 1");
    
    if (convforce(F_BFW_min,'lbf','N')/9.81 >= 0.08 * G_to/9.81) && ((convforce(F_BFW_max,'lbf','N')/9.81 <= 0.15 * G_to/9.81))
        disp("Innerhalb Limits 2");
    else
        disp("Außerhalb Limits 2");
    end
else
    disp("Außerhalb Limits 1");
end                                 %& tires.RatedLoad_Lbs_ >= F_BFW_Dyn_transf & tires.MaximumBrakingLoad_Lbs_ >= F_BFW_dyn
idz  = (tires.RatedInflation_PSI_>=200  & (tires.RatedLoad_Lbs_ >= F_reifen_BFW_max) & tires.RatedSpeed_MPH_ >= v_reifen_krit & tires.MaximumBrakingLoad_Lbs_ >= F_reifen_BFW_dyn);
tires_bugfahrwerk = tires(idz,:)


%% Speichern

Gear.F_BFW_max = F_BFW_max;

Gear.F_BFW_min = F_BFW_min;

Gear.F_BFW_dyn = F_BFW_dyn;

Gear.F_reifen_BFW_max = F_reifen_BFW_max;

Gear.F_reifen_BFW_dyn = F_reifen_BFW_dyn;

Gear.F_HFW_max = F_HFW_max_pfund;

Gear.F_reifen_HFW_max = F_reifen_HFW_max_pfund;

Gear.Reifen_HFW = tires_edit_HFW.Size(49);

Gear.Reifen_BFW = tires_bugfahrwerk.Size(1);

Gear.durchmesser = tires_edit_HFW.OutsideDiameterMax(49);

Gear.breite = tires_edit_HFW.SectionWidthMax(49);

save Fahrwerk.mat Gear

%% Reifendruck 
% Bugfahrwerk High Pressure ab 160PSI/ Aber auch 220PSI
% Hauptfahrwerk Low-Pressure