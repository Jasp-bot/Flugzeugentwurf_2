%% Funktion zur Berechnung der ISA-Atmospherendaten

function [ISA] = Berechnung_ISA(H_Start, schrittweite, H_Ende, ISAPlus)

%atmosData(:,1) = Hoehe [m]
%atmosData(:,2) = Temperatur [K]
%atmosData(:,3) = Druck [N/m^2]
%atmosData(:,4) = Schallgeschwindigkeit [m/s]
%atmosData(:,5) = Dichte [kg/m^3]
%atmosData(:,6) = relativer Dichteverlauf
%atmosData(:,7) = dynamische Viskosität
%atmosData(:,8) = kinematische Viskosität

%% Konstanten

kappa       = 1.4;              %Isentropenexponent f�r Luft [-]
R           = 287.05287;        %spezifische Gaskonstante f�r Luft [J*kg^-1*K^-1]
S           = 110.4;            %Suteherland Konstante [K]
C1          = 1.458*10^-6;      %Konstante aus Sutherland-Gesetz [Kg/(m*s*K^1/2)]
T0          = 288.15;           %Temperatur auf Meeresh�he nach ISA [K]
g0          = 9.80665;          %Fallbeschleunigung nach ISA [m/s^2]
p0          = 101325;           %Druck auf Meeresh�he nach ISA [N/m^2]
T_Tro       = 216.65;           %Temperatur Tropopause nach ISA [K]
H_Tropo_ISA = 11000;            %H�he der Tropopause nach ISA [m]
H_Str_u     = 20000;            %Beginn der unteren Stratosph�re nach ISA [m]
H_Str_m     = 32000;            %Beginn der mittleren Stratosph�re nach ISA [m]
TGradTr     = -0.0065;          %Temperaturgradient der Tropopause nach ISA [K/m]
TGradStr    = 0.001;            %Temperaturgradient der unteren Stratosph�re nach ISA [K/m]
rho_0       = 1.225;            %Dichte Kg/m^3
beta_s      = 1.458e-6;         % Konstante aus Sutherland-Grenze [kg/(m*s*sqrt(K))]



% konstanten für Strukt

ISA.R = R;
ISA.T0 = T0;
ISA.T_Tro = T_Tro;
ISA.g0 = g0;
ISA.p0 = p0;
ISA.rho_0 = rho_0;
ISA.H_Tropo_ISA = H_Tropo_ISA;
ISA.H_Str_u = H_Str_u;
ISA.H_Str_m = H_Str_m;




%% Inputs

H = (H_Start:schrittweite:H_Ende);
ISA.atmosData(:,1) = H;
ISA.H = ISA.atmosData(:,1);


%% ideale Gasgleichung bei trockener Luft (p = rho * R * T)

rho_0 = p0 /(R * T0);           %Dichte auf Meeresh�he nach ISA [kg/m^3]
a0 = sqrt(kappa * R * T0);        %Schallgeschwindigkeit auf Meeresh�he nach ISA [m/s] 


%% Berechnung des Temperaturverlaufs

%H�he der Tropopaus
T_Gesamt = T0 + ISAPlus;            %reale Umgebungstemperatur auf Meeresh�he [K] T0+TISA
Delta_T = T_Gesamt - T_Tro;                   %reale Temperaturdifferenz zwischen Meeresh�he und Tropopause [K] TGesamt-T_Tro
H_Tropo = Delta_T / abs(TGradTr);            %H�he der Tropopause in Abh�ngigkeit von T_ISAPlus [m] DeltaT/abs(TGradTr)=H_Tropo


%Temperatur in Abh�ngigkeit der H�he
ISA.atmosData(1,2) = T_Gesamt; % Hoehe bei 1m


for n = 2:length(H)                     %Schleife f�r den Verlauf der Atmosph�rendaten �ber die H�he
              if ISA.H(n,1) <= H_Tropo     %wenn H kleiner H_Tropo
                ISA.atmosData(n,2) = ISA.atmosData(n-1,2) + (TGradTr * schrittweite);
              elseif ISA.H(n,1) <= H_Str_u
                  ISA.atmosData(n,2) = T_Tro;
              else 
                  ISA.atmosData(n,2) = ISA.atmosData(n-1,2) + (TGradStr * schrittweite);
              end
                  
                   %wenn H kleiner H_Str_u
                   %wenn kleiner H_Str_m
                   %sonst... disp        
end

ISA.T = ISA.atmosData(:,2);

%% Berechnung des Druckverlaufs

%Achtung Druck ist unabh�ngig von ISA-Plus Bedingung, da jeder H�he ein
%fester Druckwert zugeordnet ist -> Wichtig f�r H�henbestimmung im Flugzeug!!!


for j = 2:length(H)                     %Schleife f�r den Verlauf der Atmosph�rendaten �ber die H�he
              if ISA.H(j,1) <= H_Tropo_ISA     %wenn H kleiner H_Tropo_ISA da P unabhängig von Temperaturdiff ISAPlus 
                  ISA.atmosData(1,3) = p0;
                  ISA.atmosData(j,3) = p0 * (1+ (TGradTr*ISA.H(j,1))/(T0))^(-(g0/(R*TGradTr)));
              elseif ISA.H(j,1) <= H_Str_u
                  p_trop_o = p0 * (1+ (TGradTr*H_Tropo_ISA)/(T0))^(-(g0/(R*TGradTr)));
                  % H_tropo_oben =  11000m = H_Tropo_ISA
                  ISA.atmosData(j,3) = p_trop_o * exp(-(g0/(R*T_Tro))*(ISA.H(j,1) - H_Tropo_ISA ));
              else 
                  p_str_u = p_trop_o * exp(-(g0/(R*T_Tro)) * (H_Str_u - H_Tropo_ISA));
                  ISA.atmosData(j,3) = p_str_u * (1 + (TGradStr*(ISA.H(j,1)-H_Str_u))/...
                      (T_Tro))^(-(g0/(R*TGradStr)));
              end
                  
                   %wenn H kleiner H_Str_u
                   %wenn kleiner H_Str_m
                   %sonst... disp        
end

ISA.p = ISA.atmosData(:,3);

%% Verlauf der Schallgeschwindigkeit


% a(H) = sqrt(kappa*R*T(H))


ISA.atmosData(:,4) = sqrt (kappa .* R .* ISA.T);  %Schallgeschwindigkeit [m/s]

ISA.a = ISA.atmosData(:,4);

%% Berechnung Dichteverlauf


% rho(H) = p(H)/(R*T(H))

ISA.atmosData(:,5) = ISA.p ./ (R * ISA.T);

ISA.rho = ISA.atmosData(:,5);

%% Berechnung relativer Dichteverlauf


% theta = rho(H)/rho_0

ISA.atmosData(:,6) = ISA.rho / rho_0;

ISA.rel_rho = ISA.atmosData(:,6);


%% dynamische Viskosität

% my(H) = (beta_s * T(H)^(3/2))/(T(H)+S)

ISA.atmosData(:,7) = (beta_s * ISA.T.^(3/2))./(ISA.T + S);

ISA.dyn_visk = ISA.atmosData(:,7);


%% kinematische Viskosität

% ny(H) = my(H)/rho(H);

ISA.atmosData(:,8) = ISA.dyn_visk ./ ISA.rho;

ISA.kin_visk = ISA.atmosData(:,8);





%% Safe für schelleren zugriff



save Ergebnisse_ISA_DATA.mat ISA;
















