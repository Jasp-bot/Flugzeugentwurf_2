% Berechnung PS9 Leitwerke

function Berechnung_PS9_Leitwerke

clc
clear all
close all



%% Laden der Werte

load Projekt_specs.mat;
load Ergebnisse_stat_Flaechenbelastung_Fluegelflaeche.mat
load Ergebnisse_Fluegel_Tank_NP.mat;
load SLW_data.mat;
load HLW_data.mat;

%%%%%%%%%%%%%%%%%%%%% Leitwerke %%%%%%%%%%%%%%%%%%%%%%%%%


%% Hoehenleitwerk

% nicht ueber regression bestimmen sondern aus werten ablesen
% % bestimmung con c_HLW (Hoehenleitwerksvolumenbeiwert)
% func_HLW = fit(HLW_data(:,1),HLW_data(:,2),'poly1');
%  %plot(func_HLW, HLW_data(:,1),HLW_data(:,2))
% % f_HLW(x) = p1*x + p2
% coef_func_HLW = coeffvalues(func_HLW);
% p1 = coef_func_HLW(1,1);
% p2 = coef_func_HLW(1,2);
% 
% c_HLW = p1 * Ma_CR + p2;



HLW.c = 0.95; % wert aus Datenblatt abgeschätzt

% betimmung der Hoehenleitwerksflaeche
HLW.r = specs.l_rumpf * 0.5;        %(0.5:0.01:0.55);

HLW.F = (HLW.c .* Ergebnisse_stat_Flaechenbelastung.F .* NP.l_mue_ges)./ HLW.r;


% Winkel phi25
HLW.phi_25 = Ergebnisse_Fluegel.phi_25_max + deg2rad(3);

% Strekung LAMBDA GROSS
% f_streckung_HLW = @(x) (-8.892 .* log(x) + 40.285);
% streckung_phi25_HLW_berechnung = f_streckung_HLW(rad2deg(phi_25_HLW));
HLW.streckung_phi25 = 4.3;
% Zuspitzung lambda (klein)

HLW.lambda = 0.4;

% pfeilung phi_VK_HLW

HLW.phi_VK = atan( tan(HLW.phi_25) - (4/HLW.streckung_phi25)*(-0.25) * ((1-HLW.lambda)/(1+HLW.lambda)));

% HLW Geometrie

HLW.b = sqrt(HLW.F .* HLW.streckung_phi25); % Lambda gross

HLW.l_i_et =  (2 .*HLW.F)./ (HLW.b .* (HLW.lambda+1));
HLW.l_a_et = HLW.lambda .* HLW.l_i_et;

HLW.phi_HK_et = atan( tan(HLW.phi_25) - (4/HLW.streckung_phi25)*(-0.25+1) * ((1-HLW.lambda)/(1+HLW.lambda)));

HLW.F_et = ((HLW.l_i_et + ((HLW.b/2)*tan(HLW.phi_HK_et)))*(HLW.b/2)...
    - ((tan(HLW.phi_VK)*(HLW.b/2)^2)/2) - (tan(HLW.phi_HK_et)*(HLW.b/2)^2)/2)*2; %% Trash rechnung %%%%%%%%%% ???? ich bin mir nicht sicher, was der kommentar bedeutet




% Itteration der Fluegelrumpfflaeche

l_i_HLW_high = 8;
l_i_HLW_low = 0.5;
dF_HLW = 1;
%phi_HK_HLW_lokal = phi_HK_HLW_et;
% specs.coanlaenge = 15;
% specs.HLW_beginn = 6.6;
kabinenversatz = specs.coanlaenge - specs.HLW_beginn;
HLW.R_rumpf_oben = (kabinenversatz/specs.coanlaenge)*specs.R_rumpf; % specs.R_rumpf - tan(atan(specs.HLW_beginn/specs.coanlaenge))*specs.HLW_beginn;
HLW.rumpfversatz_gross = HLW.R_rumpf_oben;
HLW.rumpfversatz_klein = 0.1;
HLW.s_A = HLW.b/2 - HLW.R_rumpf_oben;

%for zaehlvar = 1:10000
while abs(dF_HLW) > 0.0001
    % mittelwertsfestlegung
    l_i_HLW_try = (l_i_HLW_high + l_i_HLW_low)/(2);
    
    l_a_HLW_itt = HLW.lambda .* l_i_HLW_try;
   
    phi_HK_HLW_lokal = atan(((tan(HLW.phi_VK)*(HLW.s_A)) + l_a_HLW_itt - l_i_HLW_try)/(HLW.s_A));
    
    HLW.F_aussen = (l_i_HLW_try + ((HLW.s_A) * tan(phi_HK_HLW_lokal)))*(HLW.s_A)...
        - ((tan(HLW.phi_VK)*(HLW.s_A)^2)/2)...
        - ((tan(phi_HK_HLW_lokal)*(HLW.s_A)^2)/2); 
    
    % nebenrechnung um schnittpunkt von hinterkante mit rumpf zu betimmenrumpfversatz_gross
    % ausgehen von rumpfgeometrie
    % aufstellen von vectoren fuer interX function
    % vector geht von unterer aussenkante fluegel bis kurz vor mittelachse
    var1 = [(HLW.b/2) (HLW.rumpfversatz_klein);  
        ((kabinenversatz-tan(HLW.phi_VK)*(HLW.s_A))-l_a_HLW_itt),...
        ((kabinenversatz-tan(HLW.phi_VK)*(HLW.s_A))-l_a_HLW_itt ...
        + tan(phi_HK_HLW_lokal)*(HLW.b/2 -HLW.rumpfversatz_klein))];
    var2 = [0 specs.R_rumpf;0 specs.coanlaenge];
        

    HLW.variable_interx = InterX(var1,var2);
    HLW.schnittpunkt_x = HLW.variable_interx(1);
    HLW.schnittpunkt_y = HLW.variable_interx(2);

    % bestimmung der restlichen flaechen im Rumpf
    HLW.l_i_R = kabinenversatz - HLW.schnittpunkt_y;
    HLW.R_rumpf_unten = HLW.schnittpunkt_x;             % R_rumpf_oben = specs.R_rumpf - atan(specs.HLW_beginn/specs.coanlaenge);
    
    % Trapezflaeche Rumpf
    HLW.F_R = (1/2)*(HLW.R_rumpf_oben + HLW.R_rumpf_unten) * (HLW.l_i_R);
    
    % irreguläres dreieck außen F_dreieck
    HLW.teilflaeche_1_dreieck = (1/2) * (HLW.rumpfversatz_gross - HLW.R_rumpf_unten) * HLW.l_i_R;
    HLW.teilflaeche_2_dreieck = (tan(phi_HK_HLW_lokal) * (HLW.rumpfversatz_gross - HLW.R_rumpf_unten)^2)/(2);
    HLW.F_dreieck = HLW.teilflaeche_1_dreieck + HLW.teilflaeche_2_dreieck;


    F_HLW_itt = (HLW.F_aussen + HLW.F_dreieck + HLW.F_R)*2;
    dF_HLW = HLW.F - F_HLW_itt;
    
    if dF_HLW > 0 
        l_i_HLW_low = l_i_HLW_try;
    else 
        l_i_HLW_high = l_i_HLW_try;
    end


end

HLW.l_i = l_i_HLW_try;
HLW.l_a = l_a_HLW_itt;
HLW.phi_HK = phi_HK_HLW_lokal;





% HLW Neutralpunkte

HLW.l_mue_A = HLW.l_i * ((2/3) *((1+ HLW.lambda + HLW.lambda^2)/(1+HLW.lambda)));
% Als Annäherung wir die etscheidung getroffen das Innenstueck(Rumpfstück) nur als
% rechteck zu betrachten und um den neutralpunkt zu berechnen
HLW.lambda_R = 1;
HLW.l_mue_R = HLW.l_i_R * ((2/3) *((1+ HLW.lambda_R + HLW.lambda_R^2)/(1+HLW.lambda_R)));

HLW.y_SP_A = (HLW.b/2) * ((1/3)* ((1+2*HLW.lambda)/(1+HLW.lambda)));
HLW.y_SP_R = (HLW.b/2) * ((1/3)* ((1+2*HLW.lambda_R)/(1+HLW.lambda_R)));

HLW.streckung_phi25_A = 2*(HLW.s_A)^2 / HLW.F_aussen;
HLW.streckung_phi25_R = 2*(HLW.R_rumpf_oben)^2 / (HLW.R_rumpf_oben * HLW.l_i_R);
HLW.x_NP_A = HLW.l_i * ((1/4) + (HLW.streckung_phi25_A/12) * (1 + HLW.lambda * 2) * tan(HLW.phi_25));
HLW.x_NP_R = HLW.l_i_R * ((1/4) + (HLW.streckung_phi25_R/12) * (1 + HLW.lambda_R * 2) * tan(rad2deg(0)));


% % ges Neutalpunkte
% Alle annahmen getroffen unter der vorraussetzung, dass ein rechteck zwischen l_i_A und l_i_R exestiert

HLW.l_mue_ges = (HLW.l_mue_A*(HLW.F_aussen) + HLW.l_mue_R*(HLW.R_rumpf_oben * HLW.l_i_R)) /...
    (HLW.R_rumpf_oben * HLW.l_i_R + HLW.F_aussen);

HLW.y_SP_ges = ((HLW.y_SP_A + HLW.R_rumpf_oben )*(HLW.F_aussen) + HLW.y_SP_R*(HLW.R_rumpf_oben * HLW.l_i_R) )/...
    (HLW.F_aussen + HLW.R_rumpf_oben * HLW.l_i_R);

HLW.x_NP_ges = ((HLW.x_NP_A + kabinenversatz)*(HLW.F_aussen) +...
    (HLW.x_NP_R + kabinenversatz) * (HLW.R_rumpf_oben * HLW.l_i_R))/...
    (HLW.F_aussen + HLW.R_rumpf_oben * HLW.l_i_R);

% 
% HLW.versatz_HK = tan(HLW.phi_HK)*(HLW.b/2);

% berechchnung für abstand VK Flügel zu VK HLW
HLW.x_distanz_zu_fluegel = (HLW.r) - (kabinenversatz - HLW.x_NP_ges/2) + 8.52;

HLW.kabinenversatz = kabinenversatz;

% phi(50) HLW
HLW.phi_50 = tan((HLW.l_a/2 - HLW.l_i/2 + tan(HLW.phi_VK)*HLW.s_A)/(HLW.s_A));

HLW.l_phi50 = (HLW.b/2)/(cos(HLW.phi_50));



%% Seitenleitwerk


%% Seitenleitwerk

% [Ma   F_SLW/F     c_SLW]
% SLW_data = [0.711	0.143	0.0707;
% 0.825	0.161	0.0348;
% 0.74	0.174	0.0682;
% 0.87	0.106	0.0379;
% 0.86	0.132	0.0482;
% 0.89	0.192	0.081;
% 0.89	0.268	0.1117;
% 0.9	0.204	0.102;
% 0.875	0.191	0.0951];
% save SLW_data.mat SLW_data;


% interpolation von c_SLW
% func_SLW= fit(SLW_data(:,1),SLW_data(:,3),'poly1');
% % plot(func_SLW, SLW_data(:,1),SLW_data(:,2))
% % f_SLW(x) = p1_SLW*x + p2_SLW
% coef_func_SLW = coeffvalues(func_SLW);
% p1_SLW = coef_func_SLW(1,1);
% p2_SLW = coef_func_SLW(1,2);
% 
% c_SLW = p1_SLW * Ma_CR + p2_SLW;
SLW.c = 0.092;

% betimmung der Hoehenleitwerksflaeche

SLW.r= specs.l_rumpf .* 0.48; %(0.5:0.01:0.55);


%%% Achtung Trick um die reichtige Leitenleitwerksgrösse zu erreichen muss
%%% die Grundflaeche verdoppelt werden, da alle rehcnungen darauf ausgelegt
%%% sind mit Halber Flügelspannweite zu funktonieren
SLW.F = ((SLW.c .* Ergebnisse_stat_Flaechenbelastung.F .* Ergebnisse_Fluegel.b)./ SLW.r);
SLW.F_real = SLW.F; % Da es nur ein Seitenleitwerk gibt, die Rechnung aber für zwei Teilflügel ausgelegt ist
% Winkel phi25
SLW.phi_25 = HLW.phi_25 + deg2rad(3);

% Strekung LAMBDA GROSS
% f_streckung_SLW = @(x) (-8.892 .* log(x) + 40.285);
% streckung_phi25_SLW_berechnung = f_streckung_SLW(rad2deg(SLW.phi_25));
SLW.streckung_phi25 =3;
% Zuspitzung lambda klein

SLW.lambda = 0.5;

% pfeilung SLW.phi_VK

SLW.phi_VK = atan( tan(SLW.phi_25) - (4/SLW.streckung_phi25)*(0-0.25) * ((1-SLW.lambda)/(1+SLW.lambda)));
rad2deg(SLW.phi_VK);

% HLW Geometrie

SLW.b = sqrt(SLW.F .* SLW.streckung_phi25); % Lambda gross
% Da es nur ein Seitenleitwerk gibt, die Rechnung aber für zwei Teilflügel ausgelegt ist



SLW.l_i_et =  (2 .*SLW.F)./ (SLW.b .* (SLW.lambda+1));

SLW.l_a_et = SLW.lambda .* SLW.l_i_et;

SLW.phi_HK = atan( tan(SLW.phi_25) - (4/SLW.streckung_phi25)*(1-0.25) * ((1-SLW.lambda)/(1+SLW.lambda)));
% rad2deg(SLW.phi_HK);

SLW.F_test = ((SLW.l_i_et + ((SLW.b/2)*tan(SLW.phi_HK)))*(SLW.b/2) - ((tan(SLW.phi_VK)*(SLW.b/2)^2)/2) - (tan(SLW.phi_HK)*(SLW.b/2)^2)/2)*2; %% Trash rechnung



% Itteration der Fluegelrumpfflaeche

l_i_SLW_high = 100;
l_i_SLW_low = 1;
SLW.s_R = 2;
dF_SLW = 1;
SLW.SLW.phi_HK_lokal = SLW.phi_HK;
SLW.s_A = SLW.b - SLW.s_R;

%for zaehlvar = 1:10000
while abs(dF_SLW) > 0.001
    % mittelwertsfestlegung
    l_i_SLW_try = (l_i_SLW_high + l_i_SLW_low)/(2);
    
    % neue Berechnung

    l_a_SLW_itt = SLW.lambda .* l_i_SLW_try;
    
    
    gegenkathete_phi_HK = (tan(SLW.phi_VK) * SLW.s_A + l_a_SLW_itt) - l_i_SLW_try;

    SLW.phi_HK_lokal = atan(gegenkathete_phi_HK/SLW.s_A);
    

    SLW.F_aussen = (l_i_SLW_try + ((SLW.b/2)*tan(SLW.phi_HK_lokal)))*(SLW.b/2)...
        - ((tan(SLW.phi_HK_lokal)*(SLW.b/2)^2)/2)...
        - (tan(SLW.phi_HK_lokal)*(SLW.b/2)^2)/2; 
    SLW.F_innen = l_i_SLW_try * SLW.s_R;
    F_SLW_check = SLW.F_innen + SLW.F_aussen;
    dF_SLW = -F_SLW_check + SLW.F_real;

    if dF_SLW > 0 
        l_i_SLW_low = l_i_SLW_try;
    else 
        l_i_SLW_high = l_i_SLW_try;
    end
end
% rad2deg(SLW.phi_HK_lokal)
SLW.l_i = l_i_SLW_try;
SLW.l_a = l_a_SLW_itt;
% phi_HK_SLW_lokal_test = 0.4886
% F_SLW_aussen_test = (l_i_SLW_try + ((SLW.b/2)*tan(phi_HK_SLW_lokal_test)))*(SLW.b/2)...
%         - ((tan(phi_HK_SLW_lokal_test)*(SLW.b/2)^2)/2)...
%         - (tan(phi_HK_SLW_lokal_test)*(SLW.b/2)^2)/2; 
%     SLW.F_innen = l_i_SLW_try * SLW.s_R;
%     F_SLW_check_test = SLW.F_innen + SLW.F_aussen;



% SLW Neutralpunkte

SLW.lambda_R = 1;
SLW.lambda_lokal_A = SLW.l_a/SLW.l_i;
SLW.l_mue = SLW.l_i * ((2/3) *((1+ SLW.lambda_lokal_A + SLW.lambda_lokal_A^2)/(1+SLW.lambda_lokal_A)) );
SLW.l_mue_R = SLW.l_i * ((2/3) *((1+ SLW.lambda_R + SLW.lambda_R^2)/(1+SLW.lambda_R)) );

SLW.y_SP = (SLW.s_A) * ((1/3)* ((1+2*SLW.lambda_lokal_A)/(1+SLW.lambda_lokal_A)));
SLW.y_SP_R = (SLW.s_R) * ((1/3)* ((1+2*SLW.lambda_R)/(1+SLW.lambda_R)));

SLW.strechung_lokal_aussen = (SLW.s_A^2) /SLW.F_aussen;
SLW.x_NP = SLW.l_i * ((1/4) + (SLW.strechung_lokal_aussen /12) * (1 + SLW.lambda_lokal_A * 2) * tan(SLW.phi_25));

SLW.streckung_phi25_R =  (SLW.s_R^2) / (SLW.l_i*SLW.s_R);
SLW.phi_25_R = deg2rad(0);
SLW.x_NP_R = SLW.l_i * ((1/4)); %+ (SLW.streckung_phi25/12) * (1 + SLW.lambda_R * 2) * tan(SLW.phi_25))

% ges Neutalpunkte

SLW.l_mue_ges = (SLW.l_mue*(SLW.F_aussen) + SLW.l_mue_R*(SLW.F_innen))/(SLW.F_aussen+SLW.F_innen);
SLW.y_SP_ges = ((SLW.y_SP + SLW.s_R)*(SLW.F_aussen) + SLW.y_SP_R*(SLW.F_innen))/(SLW.F_aussen + SLW.F_innen);
SLW.x_NP_ges = (SLW.x_NP *(SLW.F_aussen) + SLW.x_NP_R*(SLW.F_innen))/(SLW.F_aussen + SLW.F_innen);

x_distanz_zu_fluegel_SLW = (SLW.r) - (SLW.x_NP_ges) + 8.52;


%% Speichern von Daten

save Ergebnisse_Leitwerke.mat HLW SLW