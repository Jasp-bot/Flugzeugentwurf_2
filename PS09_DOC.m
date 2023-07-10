%%PS09 DOC Kostenabschätzung
%Konstanten

clear all
close all
clc


load Ergebnisse_Massen_FE2.mat;
load Projekt_specs.mat
load Ergebnisse_ISA_DATA.mat


P_oe=1.235;      %Price per Kilogramm aus PS 09 [€/kg]
f_ir=0.05;       %interest rate [%]  
t_DEP=12;        %Depreciation Period [years]
f__RV=0.15;      %Residual fuel factor [%]
f_I=0.005;       %Insurance rate [%]
S_FA=50000;      %Average salery flight attendendent [€]
S_FC=[190000,237000];   %Langstrecke oder Ultra lang strecke Formel 3 PS09
n_crew=5;        %Anzahl an Crews und ersatz Crews
P_F=0.7;         %Preis je kg Sprit [€/kg]
P_LDG=0.001;     %landing gebühren [€/kg]
P_H = 0.1;       % handling fees [€/kg]
f_ATC=[1,0.7,0.6,0.5];  %Domestic Europe, transatlantic flights, far east flights, eu airport     costs of Range dependent ATC factort [€]
BT_avg=1.83;     %avarage block time supplement [h]
f_lr=50;         %Labor rate [€/h]
C_B=2;           %cost burden[1]
I_FR=0.2;        %Erlös [€/kg]
I_PAX=[550,400]; %income per seat long haul [3 KLassen, all eco]
n_pax=specs.n_pax;

%Variabeln
R= linspace(500, specs.max_range_basis_km, 20);     %Reichweite
%R=5000;


%Annahmen
R_STD=9000;

Flughoehe_CR = specs.flight_level * 10^2 ;     % in ft
hoehe_CR = round(unitsratio('m','ft')*Flughoehe_CR);
v=specs.Ma_CR*ISA.a(hoehe_CR)       %Geschwindigkeit
S_0=40;         %Standschub [t]
v_h=v/3.6;


%%Rechnung
alpha=f_ir*((1-f__RV*(1/(1+f_ir))^t_DEP)/((1-(1/(1+f_ir))^t_DEP)));

C_cap=P_oe*Ergebnisse_Massen_FE2.M_OE*(alpha+f_I);        %Capital costs

C_crew=n_crew*(S_FA*((specs.m_cargo+specs.m_pax)/5000)+S_FC(1)); %Kosten der Crew ca 3mio

C1=C_cap+C_crew    %C1 ist Routen unabhängige kosten



FC_pa=(6011/((R_STD/v_h)+BT_avg))          %PS09 Formel 10

FT_PA=6011/(1+(v_h*(BT_avg/R_STD)))    %yearly flight time 

FT=FT_PA/FC_pa

C_MRO_AF_MAT=(Ergebnisse_Massen_FE2.M_OE/1000)*(0.2*FT+13.7)+57.5; %Airframe Material maintenance costs (repair and replacement)
C_MRO_AF_PER=f_lr*(1+C_B)*((0.655+0.01*Ergebnisse_Massen_FE2.M_OE/1000))*FT+0.254+0.01*(Ergebnisse_Massen_FE2.M_OE/1000);   %Aiframe personal maintenance costs (inspection and repair)
C_MRO_ENG=specs.n_TW*(1.5*(S_0/specs.n_TW)+(30.5*FT)+10.6);

C_MRO=C_MRO_ENG+C_MRO_AF_PER+C_MRO_AF_MAT;

Fuel=fuel_range(R)

c2=FC_pa*(P_F*Fuel+((specs.m_cargo+specs.m_pax))*P_H+P_LDG*(Ergebnisse_Massen_FE2.M_TO/1000)+f_ATC(3)*R*sqrt((Ergebnisse_Massen_FE2.M_OE+Fuel)*1000)/50)+C_MRO; % PS09 Formel 4 Routen abhängige kosten

SKO=R*specs.n_pax

DOC=C1+c2           %Gesamten kosten 

SMC=DOC./SKO;

I_CAR=I_FR*(specs.m_cargo);          %%MUSS NOCHMAL ANGESEHEN WERDEN!!!! KA WAS DIE BEI FORMEL 13 DAMIT MEINEN 

n_PAX_CAR=I_CAR/I_PAX(1);

DOC_zuSKO_COR=(DOC./SKO)*(n_pax/(n_pax+n_PAX_CAR));

plot(R,SMC)



function [Fuel]=fuel_range(Range)

load Projekt_specs.mat;
load Ergebnisse_ISA_DATA.mat;
load Ergebnisse_Widerstand.mat;
load Ergebnisse_Start_Landeanforderungen.mat;
load Ergebnisse_Endwerte_Iteration_V1.mat;
load Ergebnisse_Basis_stat_m.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat;
load Ergebnisse_Fluegel_Tank_NP.mat;
load Ergebnisse_Leitwerke.mat;

% Laden der oben erstellten daten zur schnelleren benutzung muss angepasst 
% werden wenn genaueres zu den Hochauftriebshilfen bekannt ist 
load Data_PS1_High_lift_divice_Torenbeek.mat;   


%% Anfangswerte festlegen

M_TO_initial = Ergebnis_basis_m.m_To;
M_OE_initial = Ergebnis_basis_m.m_OE;
M_del_empty_initial = 140000; % Annahme die sehr nahe am ersten Schleifenduchlauf liegt
M_Zero_Fuel_initial = Ergebnis_basis_m.m_OE + specs.m_pax_all_eco + specs.m_cargo; % Annahme für all Eco Version nicht ACHTUNG
% Kappa_initial = 0.3064; % initiale annahme aus erstem Schleifenduchlauf
Zaehlvariabele = 0;
delta_M_to = 100;       % Anfangswert

% Unter der annahme, dass 50% der Strukturmasse aus CFK gefertigt werden und CFK 40% leichter ist als ALU
Technologiefaktor_ALU_CFK = 0.5 * 0.4;


   
    
    %% Berechnungen Fuel Fraction %%
    
    
    Flughoehe_CR = specs.flight_level * 10^2 ;     % in ft
    
    hoehe_CR = round(unitsratio('m','ft')*Flughoehe_CR);     % in m
    hoehe_CL = round(unitsratio('m','ft')*Flughoehe_CR*(2/3) );     % in m
    hoehe_CL_ALT = round(unitsratio('m','ft')*specs.flight_level_ALT * 10^2*(2/3) );     % in m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ?????
    hoehe_ALT = round(unitsratio('m','ft')*specs.flight_level_ALT* 10^2);     % in m
    hoehe_HLD = round(unitsratio('m', 'ft')*1500);
    
    %% Climbn_ult
    % Berechnungen für FuelFraction Methde nach Roskamp
    
    % Formel 9 PS1
    FF.S_S0_KF_CL = 0.9 * (ISA.rho(hoehe_CL)/ISA.rho_0) * ...
        exp(-0.35 * specs.Ma_CR*(2/3) * (ISA.p(hoehe_CL)/ISA.p0) * sqrt(specs.bypass));
    
    % Fromel 11 PS1
    FF.S_S0_CL = FF.S_S0_KF_CL * schub_CR.S_S0_E;
    
    % Formel 12 PS1
    FF.S_G_CL = FF.S_S0_CL * schub_CR.S0_GTo_CR *0.99;          %%%%%%%%%% schub_CR.S0_GTo_CR !!!!!!!!!!!!!!!!!! nicht sicher
    
    FF.eps_CL = 1.1 * (1/(Endwerte_Iteration.CA_CW_Clean));
    FF.v_CL = (2/3) * specs.Ma_CR * ISA.a(hoehe_CL);
    
    % Formel 8 PS1
    FF.t_CL = hoehe_CR / ((FF.S_G_CL - FF.eps_CL) * FF.v_CL);
    
    % SFC bestimmen für CL
    [FF.sfc_CL_daNh, FF.sfc_CL_1PERh, FF.sfc_CL_1PERs] = SFC(hoehe_CL, ((2/3) * specs.Ma_CR), specs.bypass);
    
    % Formel 7 PS1 (umgestellte Formel 6 PS1 um mf3 zu erhalten) 
    FF.mf3 = exp(-FF.t_CL * FF.sfc_CL_1PERs * FF.S_G_CL);
    
    %% Cruise
    
    FF.R_CR = Range .* 1000 - (FF.v_CL * FF.t_CL);
    FF.v_CR = specs.Ma_CR * ISA.a(hoehe_CR);
    [FF.sfc_CR_daNh, FF.sfc_CR_1PERh, FF.sfc_CR_1PERs] = SFC(hoehe_CR, (specs.Ma_CR), specs.bypass);
    
    FF.mf4 = 1./ exp((FF.R_CR .* FF.sfc_CR_1PERs .* (1./(Endwerte_Iteration.CA_CW_Clean))./(FF.v_CR)));
    
    FF.mf5 = 1;
    %% Diversion
    
    % Climb
    
    FF.S_S0_KF_CL_ALT = 0.9 * (ISA.rho(hoehe_CL_ALT)/ISA.rho_0) * ...
        exp(-0.35 * specs.Ma_CR*(2/3) * (ISA.p(hoehe_CL_ALT)/ISA.p0) * sqrt(specs.bypass));     %%%%%%%%%%% Nicht sicher ob richtige Machzahl
    
    FF.S_S0_CL_ALT = FF.S_S0_KF_CL_ALT * schub_CR.S_S0_E;
    
    FF.S_G_CL_ALT = FF.S_S0_CL_ALT * schub_CR.S0_GTo_CR *0.99;          %%%%%%%%%% schub_CR.S0_GTo_CR !!!!!!!!!!!!!!!!!! nicht sicher
    
    
    FF.Ma_CL_ALT = ((2/3) * specs.v_HLD) / ISA.a(hoehe_CL_ALT);
    FF.v_CL_ALT = ((2/3) * specs.v_HLD);
    
    FF.t_CL_ALT = hoehe_ALT / ((FF.S_G_CL_ALT - ((FF.eps_CL))) * FF.v_CL_ALT);  %% nicht sicher ob mit v_CL rechnen oder neuer geschwindigkeit
     
    
    [FF.sfc_CL_ALT_daNh, FF.sfc_CL_ALT_1PERh, FF.sfc_CL_ALT_1PERs] = SFC(hoehe_CL_ALT, FF.Ma_CL_ALT, specs.bypass);
    
    FF.mf6 = exp(-FF.t_CL_ALT * FF.sfc_CL_ALT_1PERs * FF.S_G_CL_ALT);
    
    % Diversion CR
    
    FF.R_ALT = specs.R_ALT ;% - (FF.v_CL_ALT * FF.t_CL_ALT); Es werden mit 200 nm im CR DIV gerechnet
    FF.v_ALT = specs.v_HLD;
    FF.Ma_CR_ALT = FF.v_ALT / ISA.a(hoehe_ALT);
    
    [FF.sfc_CR_ALT_daNh, FF.sfc_CR_ALT_1PERh, FF.sfc_CR_ALT_1PERs] = SFC(hoehe_ALT, FF.Ma_CR_ALT, specs.bypass);
    
    FF.mf7 = 1/exp((FF.R_ALT * FF.sfc_CR_ALT_1PERs * (1/(Endwerte_Iteration.CA_CW_Clean))/(FF.v_ALT)));
    
    FF.mf8 = 1;
    
    %% Holding
    
    FF.Ma_HLD = specs.v_HLD / ISA.a(hoehe_HLD);
    
    [FF.sfc_HLD_daNh, FF.sfc_HLD_1PERh, FF.sfc_HLD_1PERs] = SFC(hoehe_HLD, FF.Ma_HLD, specs.bypass);
    
    FF.mf9 = 1/(exp(specs.t_HLD * FF.sfc_HLD_1PERs * (1/Endwerte_Iteration.CA_CW_Clean)));
    
    
    %% Fuel Fraction ges
    
    
    % Nach Roskam
    FF.mf0 = 0.992; 
    FF.mf1 = 0.996;
    FF.mf2 = 0.996;
    FF.mf10 = 1; % hat kristof gesagt
    
    %%%%%%%%%%%%%%%%%%%%%%%%  mfi muss noch mal überprüft werden, ich habe
    %%%%%%%%%%%%%%%%%%%%%%%%  jetzt mit alles Massenanteilen gerechnet, von 0
    %%%%%%%%%%%%%%%%%%%%%%%%  bis 10 nicht von 2 bis 10. Muss eventuell nochmal
    %%%%%%%%%%%%%%%%%%%%%%%%  angepasst werden
    % Matrix für Produkt auf mf0 bis mf10
    FF.mfi = [FF.mf0, FF.mf1, FF.mf2, FF.mf3,...
        FF.mf4, FF.mf5, FF.mf6, FF.mf7,...
        FF.mf8, FF.mf9, FF.mf10]; 
    
    FF.Mff = prod(FF.mfi(3:11)); % Faktorprodukt  
    FF.mf_oC = (1-FF.Mff) * M_TO_initial; %%%%%%%%%%%%%%%%%%%%% Muss verändert werden für iteration
    
    FF.Mff_2_10 = prod(FF.mfi(3:11)); % Gesamtanteil Mff ab TO
    FF.Mff_2_5 = prod(FF.mfi(3:6)); % Nur Reisefluganteil ohne DIV
    
    % Kraftstoffmassenfaktor neu
    FF.Kappa_ges = 1- FF.Mff_2_10 + 0.05 * (1 - FF.Mff_2_5); 
    
    M_take_off_initial.M_fuel = M_TO_initial * FF.Kappa_ges;

    Fuel=M_take_off_initial.M_fuel
end






