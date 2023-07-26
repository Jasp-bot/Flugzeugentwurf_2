%% Visulaisierung zur Schwerpunktsberechnung FE2

% Be- und Entladung 3-Klassen
% Be- und Entladung Maximale Nutzlast (AllEco)
% Be- und Entladung Überführung
% Vergleich Schwerpunktlage Rumpf, Flügel, Gesamt
% Einzelschwerpunkte Rumpf
% Einzelschwerpunkte Flügel

% siehe Berechnung_FE2_PS2_Schwerpunkt

load SchwerpunktPlot.mat

%% Plotten
% BELADUNG 3-KLASSEN
figure(1)
hold on 
grid on
xlim([-30 60])
ylim([Ergebnisse_Massen_FE2.M_OE-10000 Ergebnisse_Massen_FE2.M_TO+10000])

% Fuel
plot([Betankung.P1(1)*100 Betankung.P2(1)*100], [Betankung.P1(2) Betankung.P2(2)],"b-",'LineWidth', 2)
plot([Betankung.P2(1)*100 Betankung.P3(1)*100], [Betankung.P2(2) Betankung.P3(2)],"bo-",'LineWidth', 2)
% Pax Außen
plot(CG_Shift_Outer*100/NP.l_mue_ges, NewMassCounter,"rx-",'LineWidth', 2)
plot(BackwardsCG_Shift_Outer*100/NP.l_mue_ges, BackwardsNewMassCounter,"mx-",'LineWidth', 2)
% Pax Innen
plot(CG_Shift_Inner*100/NP.l_mue_ges, InnerMassCounter,"rx-",'LineWidth', 2)
plot(BackwardsCG_Shift_Inner*100/NP.l_mue_ges, BackwardsInnerMassCounter,"mx-",'LineWidth', 2)
% Fracht Vorn nach Hinten
plot([BackwardsCG_Startposition_Innen*100/NP.l_mue_ges, CG_Fracht.CG_BeladenVorn*100/NP.l_mue_ges],[BackwardsMass_Shift_Inner, CG_Fracht.Masse_BeladenVorn],"g-",'LineWidth', 2)
plot([CG_Fracht.CG_BeladenVorn*100/NP.l_mue_ges, CG_Fracht.CG_BeladenHinten*100/NP.l_mue_ges],[CG_Fracht.Masse_BeladenVorn, CG_Fracht.Masse_BeladenHinten],"go-",'LineWidth', 2)
% Fracht Hinten nach Vorn
plot([BackwardsCG_Startposition_Innen*100/NP.l_mue_ges, CG_Fracht.CG_BeladenHinten1*100/NP.l_mue_ges],[BackwardsMass_Shift_Inner, CG_Fracht.Masse_BeladenHinten1],"g--",'LineWidth', 2)
plot([CG_Fracht.CG_BeladenHinten1*100/NP.l_mue_ges, CG_Fracht.CG_BeladenVorn1*100/NP.l_mue_ges],[CG_Fracht.Masse_BeladenHinten1, CG_Fracht.Masse_BeladenVorn1],"go--",'LineWidth', 2)
%Grenzen
plot([BFWL.x_CG_BFW_Min_MAC_Prozent*100,BFWL.x_CG_BFW_Min_MAC_Prozent*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"-k",'LineWidth', 1.5)
plot(BFWL.x_CG_BFW_Max_MAC_Prozent*100,BFWL.MomentanMasse,"--k",'LineWidth', 1.5)
plot([LS.Laengsstabilitaet*100,LS.Laengsstabilitaet*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"--b",'LineWidth', 1.5)
plot([StatStab.CG_sigma_x*100,StatStab.CG_sigma_x*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],":b",'LineWidth', 1.5)
plot([-100,100],[Ergebnisse_Massen_FE2.M_OE,Ergebnisse_Massen_FE2.M_OE],'Color','#607B8B','LineWidth', 1.5)
plot([-100,100],[Ergebnisse_Massen_FE2.M_TO,Ergebnisse_Massen_FE2.M_TO],'Color',"#104E8B",'LineWidth', 1.5)

title('Beladung 3-Klassenbestuhlung','FontSize',20);
xlabel('X^{MAC}_{SP}/l_{\mu} [%]','FontSize',16)
ylabel('Masse [kg]','FontSize',16)
legend('Betankung Innentank','Betankung Außentank', 'Außenreihen von vorn','Außenreihen von hinten','Innenreihen von vorn','Innenreihen von hinten','Fracht vorn voll','Fracht hinten auffüllen','Fracht hinten','Fracht vorn voll','Minimale BFWL','Maximale BFWL','Längsstabilität','NP -5%','M_{OE}','M_{TO}','Location','southwest')

hold off

% %---------------------------
% ENTLADUNG 3-KLASSE
figure(2)
hold on 
grid on
xlim([-30 60])
ylim([Ergebnisse_Massen_FE2.M_OE-10000 Ergebnisse_Massen_FE2.M_TO+10000])

% Fuel
plot([Enttankung.P1(1)*100 Enttankung.P2(1)*100], [Enttankung.P1(2) Enttankung.P2(2)],"b-",'LineWidth', 2)
plot([Enttankung.P2(1)*100 Enttankung.P3(1)*100], [Enttankung.P2(2) Enttankung.P3(2)],"bo-",'LineWidth', 2)
% Fracht Vorn nach Hinten
plot([CG_Fracht.CG_EntladungStart*100/NP.l_mue_ges, CG_Fracht.CG_EntladenVorn*100/NP.l_mue_ges], [CG_Fracht.Masse_EntladungStart CG_Fracht.Masse_EntladungVorn],"g-",'LineWidth', 2)
plot([CG_Fracht.CG_EntladenVorn*100/NP.l_mue_ges, CG_Fracht.CG_EntladenHinten*100/NP.l_mue_ges], [CG_Fracht.Masse_EntladungVorn, CG_Fracht.Masse_EntladungHinten],"go-",'LineWidth', 2)
% Fracht Hinten nach Vorn
plot([CG_Fracht.CG_EntladungStart*100/NP.l_mue_ges, CG_Fracht.CG_EntladenHinten2*100/NP.l_mue_ges], [CG_Fracht.Masse_EntladungStart CG_Fracht.Masse_EntladungHinten2],"g--",'LineWidth', 2)
plot([CG_Fracht.CG_EntladenHinten2*100/NP.l_mue_ges, CG_Fracht.CG_EntladenVorn2*100/NP.l_mue_ges], [CG_Fracht.Masse_EntladungHinten2, CG_Fracht.Masse_EntladungVorn2],"go--",'LineWidth', 2)
% Pax Außen
plot(CG_PaxFBA.CG_Shift_Outer*100/NP.l_mue_ges, CG_PaxFBA.NewMassCounter,"rx-",'LineWidth', 2)
plot(CG_PaxBFA.BackwardsCG_Shift_Outer*100/NP.l_mue_ges, CG_PaxBFA.BackwardsNewMassCounter,"mx-",'LineWidth', 2)
% Pax Innen
plot(CG_PaxFBI.CG_Shift_Inner*100/NP.l_mue_ges, CG_PaxFBI.InnerMassCounter,"rx-",'LineWidth', 2)
plot(CG_PaxBFI.BackwardsCG_Shift_Inner*100/NP.l_mue_ges, CG_PaxBFI.BackwardsInnerMassCounter,"mx-",'LineWidth', 2)
%Grenzen
plot([BFWL.x_CG_BFW_Min_MAC_Prozent*100,BFWL.x_CG_BFW_Min_MAC_Prozent*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"-k",'LineWidth', 1.5)
plot(BFWL.x_CG_BFW_Max_MAC_Prozent*100,BFWL.MomentanMasse,"--k",'LineWidth', 1.5)
plot([LS.Laengsstabilitaet*100,LS.Laengsstabilitaet*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"--b",'LineWidth', 1.5)
plot([StatStab.CG_sigma_x*100,StatStab.CG_sigma_x*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],":b",'LineWidth', 1.5)
plot([-100,100],[Ergebnisse_Massen_FE2.M_OE,Ergebnisse_Massen_FE2.M_OE],'Color','#607B8B','LineWidth', 1.5)
plot([-100,100],[Ergebnisse_Massen_FE2.M_TO,Ergebnisse_Massen_FE2.M_TO],'Color',"#104E8B",'LineWidth', 1.5)

title('Entladung 3-Klassenbestuhlung','FontSize',20);
xlabel('X^{MAC}_{SP}/l_{\mu} [%]','FontSize',16)
ylabel('Masse [kg]','FontSize',16)
legend('Enttankung Innentank','Enttankung Außentank','Entladung Fracht vorn','Entladung Fracht hinten','Entladung Fracht hinten','Entladung Fracht vorn','Außenreihen von vorn','Außenreihen von hinten','Innenreihen von vorn','Innenreihen von hinten','Minimale BFWL','Maximale BFWL','Längsstabilität','NP -5%','M_{OE}','M_{TO}','Location','southwest')

hold off

%----------------------------------------------------------
% BELADUNG ALL-ECO
figure(3)
hold on 
grid on
xlim([-30 60])
ylim([Ergebnisse_Massen_FE2.M_OE-10000 Ergebnisse_Massen_FE2.M_TO+10000])

% Fuel
plot([BetankungAllEco.P1(1)*100 BetankungAllEco.P2(1)*100], [BetankungAllEco.P1(2) BetankungAllEco.P2(2)],"b-",'LineWidth', 2)
plot([BetankungAllEco.P2(1)*100 BetankungAllEco.P3(1)*100], [BetankungAllEco.P2(2) BetankungAllEco.P3(2)],"bo-",'LineWidth', 2)
% Pax Außen
plot(AllEcoPax.CG_Shift_Outer*100/NP.l_mue_ges, AllEcoPax.NewMassCounter,"rx-",'LineWidth', 2)
plot(AllEcoPax.BackwardsCG_Shift_Outer*100/NP.l_mue_ges, AllEcoPax.BackwardsNewMassCounter,"mx-",'LineWidth', 2)
% Pax Innen
plot(AllEcoPax.CG_Shift_Inner*100/NP.l_mue_ges, AllEcoPax.InnerMassCounter,"rx-",'LineWidth', 2)
plot(AllEcoPax.BackwardsCG_Shift_Inner*100/NP.l_mue_ges, AllEcoPax.BackwardsInnerMassCounter,"mx-",'LineWidth', 2)
% Fracht Vorn nach Hinten
plot([AllEcoPax.BackwardsCG_Startposition_Innen*100/NP.l_mue_ges, CG_Fracht_AllEco.CG_BeladenVorn*100/NP.l_mue_ges],[AllEcoPax.BackwardsMass_Shift_Inner, CG_Fracht_AllEco.Masse_BeladenVorn],"g-",'LineWidth', 2)
plot([CG_Fracht_AllEco.CG_BeladenVorn*100/NP.l_mue_ges, CG_Fracht_AllEco.CG_BeladenHinten*100/NP.l_mue_ges],[CG_Fracht_AllEco.Masse_BeladenVorn, CG_Fracht_AllEco.Masse_BeladenHinten],"go-",'LineWidth', 2)
% Fracht Hinten nach Vorn
plot([AllEcoPax.BackwardsCG_Startposition_Innen*100/NP.l_mue_ges, CG_Fracht_AllEco.CG_BeladenHinten1*100/NP.l_mue_ges],[AllEcoPax.BackwardsMass_Shift_Inner, CG_Fracht_AllEco.Masse_BeladenHinten1],"g--",'LineWidth', 2)
plot([CG_Fracht_AllEco.CG_BeladenHinten1*100/NP.l_mue_ges, CG_Fracht_AllEco.CG_BeladenVorn1*100/NP.l_mue_ges],[CG_Fracht_AllEco.Masse_BeladenHinten1, CG_Fracht_AllEco.Masse_BeladenVorn1],"go--",'LineWidth', 2)
%Grenzen
plot([BFWL.x_CG_BFW_Min_MAC_Prozent*100,BFWL.x_CG_BFW_Min_MAC_Prozent*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"-k",'LineWidth', 1.5)
plot(BFWL.x_CG_BFW_Max_MAC_Prozent*100,BFWL.MomentanMasse,"--k",'LineWidth', 1.5)
plot([LS.Laengsstabilitaet*100,LS.Laengsstabilitaet*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"--b",'LineWidth', 1.5)
plot([StatStab.CG_sigma_x*100,StatStab.CG_sigma_x*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],":b",'LineWidth', 1.5)
plot([-100,100],[Ergebnisse_Massen_FE2.M_OE,Ergebnisse_Massen_FE2.M_OE],'Color','#607B8B','LineWidth', 1.5)
plot([-100,100],[Ergebnisse_Massen_FE2.M_TO,Ergebnisse_Massen_FE2.M_TO],'Color',"#104E8B",'LineWidth', 1.5)

title('Beladung All-Eco','FontSize',20);
xlabel('X^{MAC}_{SP}/l_{\mu} [%]','FontSize',16)
ylabel('Masse [kg]','FontSize',16)
legend('Betankung Innentank','Betankung Außentank', 'Außenreihen von vorn','Außenreihen von hinten','Innenreihen von vorn','Innenreihen von hinten','Fracht vorn voll','Fracht hinten auffüllen','Fracht hinten','Fracht vorn voll','Minimale BFWL','Maximale BFWL','Längsstabilität','NP -5%','M_{OE}','M_{TO}','Location','southwest')

hold off

%---------------------------
% ENTLADUNG ALL-ECO
figure(4)
hold on 
grid on
xlim([-30 60])
ylim([Ergebnisse_Massen_FE2.M_OE-10000 Ergebnisse_Massen_FE2.M_TO+10000])

% Fuel
plot([Enttankung_AllEco.P1(1)*100 Enttankung_AllEco.P2(1)*100], [Enttankung_AllEco.P1(2) Enttankung_AllEco.P2(2)],"b-",'LineWidth', 2)
plot([Enttankung_AllEco.P2(1)*100 Enttankung_AllEco.P3(1)*100], [Enttankung_AllEco.P2(2) Enttankung_AllEco.P3(2)],"bo-",'LineWidth', 2)
% Fracht Vorn nach Hinten
plot([CG_Fracht_AllEco.CG_EntladungStart*100/NP.l_mue_ges, CG_Fracht_AllEco.CG_EntladenVorn*100/NP.l_mue_ges], [CG_Fracht_AllEco.Masse_EntladungStart CG_Fracht_AllEco.Masse_EntladungVorn],"g-",'LineWidth', 2)
plot([CG_Fracht_AllEco.CG_EntladenVorn*100/NP.l_mue_ges, CG_Fracht_AllEco.CG_EntladenHinten*100/NP.l_mue_ges], [CG_Fracht_AllEco.Masse_EntladungVorn, CG_Fracht_AllEco.Masse_EntladungHinten],"go-",'LineWidth', 2)
% Fracht Hinten nach Vorn
plot([CG_Fracht_AllEco.CG_EntladungStart*100/NP.l_mue_ges, CG_Fracht_AllEco.CG_EntladenHinten2*100/NP.l_mue_ges], [CG_Fracht_AllEco.Masse_EntladungStart CG_Fracht_AllEco.Masse_EntladungHinten2],"g--",'LineWidth', 2)
plot([CG_Fracht_AllEco.CG_EntladenHinten2*100/NP.l_mue_ges, CG_Fracht_AllEco.CG_EntladenVorn2*100/NP.l_mue_ges], [CG_Fracht_AllEco.Masse_EntladungHinten2, CG_Fracht_AllEco.Masse_EntladungVorn2],"go--",'LineWidth', 2)
% Pax Außen
plot(CG_PaxFBA_AllEco.CG_Shift_Outer*100/NP.l_mue_ges, CG_PaxFBA_AllEco.NewMassCounter,"rx-",'LineWidth', 2)
plot(CG_PaxBFA_AllEco.BackwardsCG_Shift_Outer*100/NP.l_mue_ges, CG_PaxBFA_AllEco.BackwardsNewMassCounter,"mx-",'LineWidth', 2)
% Pax Innen
plot(CG_PaxFBI_AllEco.CG_Shift_Inner*100/NP.l_mue_ges, CG_PaxFBI_AllEco.InnerMassCounter,"rx-",'LineWidth', 2)
plot(CG_PaxBFI_AllEco.BackwardsCG_Shift_Inner*100/NP.l_mue_ges, CG_PaxBFI_AllEco.BackwardsInnerMassCounter,"mx-",'LineWidth', 2)
%Grenzen
plot([BFWL.x_CG_BFW_Min_MAC_Prozent*100,BFWL.x_CG_BFW_Min_MAC_Prozent*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"-k",'LineWidth', 1.5)
plot(BFWL.x_CG_BFW_Max_MAC_Prozent*100,BFWL.MomentanMasse,"--k",'LineWidth', 1.5)
plot([LS.Laengsstabilitaet*100,LS.Laengsstabilitaet*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"--b",'LineWidth', 1.5)
plot([StatStab.CG_sigma_x*100,StatStab.CG_sigma_x*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],":b",'LineWidth', 1.5)
plot([-100,100],[Ergebnisse_Massen_FE2.M_OE,Ergebnisse_Massen_FE2.M_OE],'Color','#607B8B','LineWidth', 1.5)
plot([-100,100],[Ergebnisse_Massen_FE2.M_TO,Ergebnisse_Massen_FE2.M_TO],'Color',"#104E8B",'LineWidth', 1.5)

title('Entladung All-Eco','FontSize',20);
xlabel('X^{MAC}_{SP}/l_{\mu} [%]','FontSize',16)
ylabel('Masse [kg]','FontSize',16)
legend('Enttankung Innentank','Enttankung Außentank','Entladung Fracht vorn','Entladung Fracht hinten','Entladung Fracht hinten','Entladung Fracht vorn','Außenreihen von vorn','Außenreihen von hinten','Innenreihen von vorn','Innenreihen von hinten','Minimale BFWL','Maximale BFWL','Längsstabilität','NP -5%','M_{OE}','M_{TO}','Location','southwest')

hold off

% %-------------------------------
% Be- & Entladung ÜBERFÜHRUNG
figure(5)
hold on 
grid on
xlim([-30 60])
ylim([Ergebnisse_Massen_FE2.M_OE-10000 Ergebnisse_Massen_FE2.M_TO+10000])

% Fuel
plot([CG_Ueberfuehrung.P1(1)*100 CG_Ueberfuehrung.P2(1)*100], [CG_Ueberfuehrung.P1(2) CG_Ueberfuehrung.P2(2)],"b-",'LineWidth', 2)
plot([CG_Ueberfuehrung.P2(1)*100 CG_Ueberfuehrung.P3(1)*100], [CG_Ueberfuehrung.P2(2) CG_Ueberfuehrung.P3(2)],"bo-",'LineWidth', 2)
plot([Enttankung_Ueberfuehrung.P1(1)*100 Enttankung_Ueberfuehrung.P2(1)*100], [Enttankung_Ueberfuehrung.P1(2) Enttankung_Ueberfuehrung.P2(2)],"c-",'LineWidth', 2)
plot([Enttankung_Ueberfuehrung.P2(1)*100 Enttankung_Ueberfuehrung.P3(1)*100], [Enttankung_Ueberfuehrung.P2(2) Enttankung_Ueberfuehrung.P3(2)],"co-",'LineWidth', 2)
%Grenzen
plot([BFWL.x_CG_BFW_Min_MAC_Prozent*100,BFWL.x_CG_BFW_Min_MAC_Prozent*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"-k",'LineWidth', 1.5)
plot(BFWL.x_CG_BFW_Max_MAC_Prozent*100,BFWL.MomentanMasse,"--k",'LineWidth', 1.5)
plot([LS.Laengsstabilitaet*100,LS.Laengsstabilitaet*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],"--b",'LineWidth', 1.5)
plot([StatStab.CG_sigma_x*100,StatStab.CG_sigma_x*100],[Ergebnisse_Massen_FE2.M_OE Ergebnisse_Massen_FE2.M_TO+10000],":b",'LineWidth', 1.5)
plot([-100,100],[Ergebnisse_Massen_FE2.M_OE,Ergebnisse_Massen_FE2.M_OE],'Color','#607B8B','LineWidth', 1.5)
plot([-100,100],[Ergebnisse_Massen_FE2.M_TO,Ergebnisse_Massen_FE2.M_TO],'Color',"#104E8B",'LineWidth', 1.5)
plot([-100,100],[CG_Ueberfuehrung.m_TO,CG_Ueberfuehrung.m_TO],'Color',"#104E8B",'LineWidth', 1.5,'LineStyle','--')

title('Be- und Entladung Überführung','FontSize',20);
xlabel('X^{MAC}_{SP}/l_{\mu} [%]','FontSize',16)
ylabel('Masse [kg]','FontSize',16)
legend('Betankung Innentank','Betankung Außentank','Enttankung Innentank','Enttankung Außentank','Minimale BFWL','Maximale BFWL','Längsstabilität','NP -5%','M_{OE}','M_{TO,Design}','M_{TO,Überführung}','Location','southwest')

hold off

% %-------------------------------------
% PLOT SCHWERPUNKTLAGE
figure(6)
hold on
grid on
xlim([30 50])

% RUMPF
% Mittelpunkt der Ellipse
center_x = specs.l_rumpf/2;
center_y = 0;
% Erstellen der Ellipse
theta = linspace(0, 2*pi, 100); % Winkel für die Ellipsenpunkte
x = center_x + (specs.l_rumpf/2) * cos(theta);
y = center_y + (specs.D_rumpf/2) * sin(theta);

% ERSATZFLÜGEL
x1 = Wing_MAC.XMAC;
y1 = -24;
x2 = Wing_MAC.XMAC + NP.l_mue_ges;
y2 = 24;


% PLOTTEN
% Rumpf
plot(x, y, 'k', 'LineWidth', 1.5);
% Ersatzflügel
rectangle('Position', [x1, y1, x2-x1, y2-y1], 'EdgeColor', 'r', 'LineWidth', 1);
% Flügel
plot([Wing_Position1, Wing_Position1 + (DT.s_A+DT.s_I+DT.s_R)*tan(DT.phi_VK_max)], [0, (DT.s_A+DT.s_I+DT.s_R)], 'k', 'LineWidth', 1.5);
plot([Wing_Position1, Wing_Position1 + (DT.s_A+DT.s_I+DT.s_R)*tan(DT.phi_VK_max)], [0, -(DT.s_A+DT.s_I+DT.s_R)], 'k', 'LineWidth', 1.5);
plot([Wing_Position1+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max), Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max)],[DT.s_R+DT.s_I, DT.s_R+DT.s_I], 'k', 'LineWidth', 1.5);
plot([Wing_Position1+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max), Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max)],[-(DT.s_R+DT.s_I), -(DT.s_R+DT.s_I)], 'k', 'LineWidth', 1.5);
plot([Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max),DT.l_a+Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max)],[-(DT.s_R+DT.s_I+DT.s_A),-(DT.s_R+DT.s_I+DT.s_A)], 'k', 'LineWidth', 1.5);
plot([Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max),DT.l_a+Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max)],[DT.s_R+DT.s_I+DT.s_A,DT.s_R+DT.s_I+DT.s_A], 'k', 'LineWidth', 1.5);
plot([Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max),DT.l_a+Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max)],[-(DT.s_R+DT.s_I),-(DT.s_R+DT.s_I+DT.s_A)], 'k', 'LineWidth', 1.5);
plot([Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max),DT.l_a+Wing_Position1+(DT.s_R+DT.s_I+DT.s_A)*tan(DT.phi_VK_max)],[DT.s_R+DT.s_I,DT.s_R+DT.s_I+DT.s_A], 'k', 'LineWidth', 1.5);
plot([Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max),Wing_Position1+DT.l_i_A+(DT.s_R+DT.s_I)*tan(DT.phi_VK_max)],[DT.s_R+DT.s_I,-(DT.s_R+DT.s_I)], 'k', 'LineWidth', 1.5);
% Schwerpunkte Tank
% plot([CG_Fuel_X.Aussentrapez+(DT.s_I*tan(DT.phi_VK_max))+Wing_Position1+DT.s_R*tan(DT.phi_VK_max)+(0.15*DT.l_i_R),CG_Fuel_X.Aussentrapez+(DT.s_I*tan(DT.phi_VK_max))+Wing_Position1+DT.s_R*tan(DT.phi_VK_max)+(0.15*DT.l_i_R)],[-30,30], '-b');
%plot([Wing_MAC.XMAC+CG_Fuel_X.CG_Tank_MAC,Wing_MAC.XMAC+CG_Fuel_X.CG_Tank_MAC],[-30,30], '--m');
%plot([Wing_Position1+CG_Fuel_X.CG_Tank_FG,Wing_Position1+CG_Fuel_X.CG_Tank_FG],[-30,30], ':g');

%plot([Wing_MAC.XMAC+LS.x_MainGear_MAC*NP.l_mue_ges,Wing_MAC.XMAC+LS.x_MainGear_MAC*NP.l_mue_ges],[-30,30], '--m');
%plot([Wing_Position1+CG_Data_Wing.MainGear(2),Wing_Position1+CG_Data_Wing.MainGear(2)],[-30,30], ':g');

% CGs
plot(CG_Rumpf_X,0,'xb','LineWidth', 2);
plot(Wing_Position1+CG_Wing_X,0,'ob','LineWidth', 2);
plot(LS.x_MainGear_MAC*NP.l_mue_ges +Wing_MAC.XMAC, 0, 'xr','LineWidth', 2);
%plot(CG_Data.Bugfahrwerk(2)*specs.l_rumpf,0,'*k','LineWidth', 2);
plot(CG_Gesamt_x + Wing_MAC.XMAC,0,'og','LineWidth', 3)
plot([Wing_Position1+CG_Data_Wing.Fluegel(2),Wing_Position1+CG_Data_Wing.Fluegel(2)],[-30,30],'--r')
plot(Wing_MAC.XMAC+(X_NP_durch_l_mue*NP.l_mue_ges),0,'dm','LineWidth', 2);
plot(Wing_MAC.XMAC+(X_NP_OH_durch_l_mue*NP.l_mue_ges),0,'*m','LineWidth', 2);

axis equal;
legend('','','','','','','','','','','CG Rumpf','CG Flügelgruppe','CG HFW', 'CG Gesamt','','Gesamtneutralpunkt','NP oH')

% Achsenbeschriftungen und Titel
xlabel('x [m]');
ylabel('y [m]');
title('Schwerpunktlagen');
hold off

%----------------------------------------
% Plot CG Positionen Rumpf
img = imread('RumpfSide.png');
Rumpfplot.scale_x = 2617.5;
Rumpfplot.scale_z = 217/specs.D_rumpf;
Rumpfplot.delta_z = 144 + specs.D_rumpf*0.5*Rumpfplot.scale_z;
%Rumpfplot.CG_x = [CG_Data.Rumpf(2);CG_Data.HLW(2);CG_Data.SLW(2);CG_Data.Bugfahrwerk(2);CG_Data.APU(2);CG_Data.CockpitInstruments(2);CG_Data.Furnishing(2);CG_Data.CrewProvisions(2);CG_Data.PassengerCabinSupplies(2);CG_Data.WaterToiletChem(2);CG_Data.SafetyEq(2);CG_Data.Seating(2)];
Rumpfplot.CG_x = CG_DataMatrix(:,2);
%Rumpfplot.CG_z = [CG_Data.Rumpf(4);CG_Data.HLW(4);CG_Data.SLW(4);CG_Data.Bugfahrwerk(4);CG_Data.APU(4);CG_Data.CockpitInstruments(4);CG_Data.Furnishing(4);CG_Data.CrewProvisions(4);CG_Data.PassengerCabinSupplies(4);CG_Data.WaterToiletChem(4);CG_Data.SafetyEq(4);CG_Data.Seating(4)];
Rumpfplot.CG_z = CG_DataMatrix(:,4);

jetcustom=jet(length(Rumpfplot.CG_x)+2);
figure, imshow(flipud(img));
set(gca,'YDir','normal')
set(gca, 'ColorOrder', jetcustom , 'NextPlot', 'replacechildren');
hold on;
plot(Rumpfplot.CG_x(1)*Rumpfplot.scale_x, Rumpfplot.CG_z(1)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(2)*Rumpfplot.scale_x, Rumpfplot.CG_z(2)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(3)*Rumpfplot.scale_x, Rumpfplot.CG_z(3)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(4)*Rumpfplot.scale_x, Rumpfplot.CG_z(4)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(5)*Rumpfplot.scale_x, Rumpfplot.CG_z(5)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(6)*Rumpfplot.scale_x, Rumpfplot.CG_z(6)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(7)*Rumpfplot.scale_x, Rumpfplot.CG_z(7)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(8)*Rumpfplot.scale_x, Rumpfplot.CG_z(8)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(9)*Rumpfplot.scale_x, Rumpfplot.CG_z(9)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(10)*Rumpfplot.scale_x, Rumpfplot.CG_z(10)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(11)*Rumpfplot.scale_x, Rumpfplot.CG_z(11)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(12)*Rumpfplot.scale_x, Rumpfplot.CG_z(12)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(13)*Rumpfplot.scale_x, Rumpfplot.CG_z(13)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(14)*Rumpfplot.scale_x, Rumpfplot.CG_z(14)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(15)*Rumpfplot.scale_x, Rumpfplot.CG_z(15)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(16)*Rumpfplot.scale_x, Rumpfplot.CG_z(16)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(Rumpfplot.CG_x(17)*Rumpfplot.scale_x, Rumpfplot.CG_z(17)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(CG_Fracht.X_vorn_Prozent*Rumpfplot.scale_x, (-1.7)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(CG_Fracht.X_hinten_Prozent*Rumpfplot.scale_x, (-1.7)*Rumpfplot.scale_z+Rumpfplot.delta_z,'Marker','x','MarkerSize',8,'LineWidth',2)

legend('Rumpf','HLW','SLW','Bugfahrwerk','APU','Cockpit Instruments','Hydraulics Front','Hydraulics Back','Furnishing','Anti Ice & AC','Misc','Crew Provisions','Passenger Cabin Supplies','Water & Toilet Chems','Safety Equipment','Seating','Residual Fuel','Fracht vorn','Fracht hinten','Location','eastoutside')
hold off;

%---------------------------------------------
% Plot CG Positionen Flügel

img2 = imread('WingGroup.png');
Wingplot.scale_x = 581/(DT.s_R+DT.s_I+DT.s_A);
Wingplot.scale_z = 201/DT.l_i_I;

Wingplot.CG_x = CG_DataMatrix_Wing(:,2);
Wingplot.CG_z = CG_DataMatrix_Wing(:,4);

jetcustom=jet(length(Wingplot.CG_x)+3);
figure, imshow(flipud(img2));
%figure, imshow(img2);
set(gca,'YDir','normal')
set(gca, 'ColorOrder', jetcustom , 'NextPlot', 'replacechildren');
hold on;
plot(581, Wingplot.CG_x(1)*Wingplot.scale_x,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(581, Wingplot.CG_x(2)*Wingplot.scale_x,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(581, Wingplot.CG_x(3)*Wingplot.scale_x,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(581, Wingplot.CG_x(4)*Wingplot.scale_x,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(581, Wingplot.CG_x(5)*Wingplot.scale_x,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(560, CG_Fuel_X.Rumpf_FG*Wingplot.scale_x,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(465, CG_Fuel_X.Innentrapez_FG*Wingplot.scale_x,'Marker','x','MarkerSize',8,'LineWidth',2)
plot(230, CG_Fuel_X.Aussentrapez_FG*Wingplot.scale_x,'Marker','x','MarkerSize',8,'LineWidth',2)

legend('Flügel','Main Gear','Surface Controls','Engine','Nacelle','Tank Rumpf','Tank Innen','Tank Außen','Location','eastoutside')
hold off;