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

if (F_BFW_min >= 0.06*G_to) || (F_BFW_max <= 0.2*G_to)
    disp("Within Limits");
end

%% Hauptfahrwerk

F_HFW_max = (G_to *((l_BFW + l_HFW) - l_HFW_min)) / (n_FWB * (l_BFW + l_HFW));


%% Entwurfslasten der Reifen    % PRO REIFEN !!
S_Reserve = 0.25
S_FW = 1 + ((0.07 + S_Reserve)/100);

F_reifen_BFW_max = (F_BFW_max/n_reifen_BFW) * S_FW;
F_reifen_BFW_dyn = (F_BFW_dyn/n_reifen_BFW) * S_FW;

F_reifen_HFW_max = (F_HFW_max/n_reifen_HFW) * S_FW;

%% Reifengeschwindigkeiten max

v_reifen_max_LDG = 1.3 * v_s_LDG;
v_reifen_max_TO = 1.2 * v_s_TO;

%% Reifendruck 



%% Reifenwahl


%% BEstimmung LCN 



% ESWL bestimmen

S_t = 


