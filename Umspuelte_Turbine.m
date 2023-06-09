%%
%%   Quelle: E. Torenbeek: "Synthesis of Subsonic Airplane Design", Appendix B S446-449
%%Umströmte Fläche für EEEEEEIIIIINNNNNEEEEEEE Gondel

clc
clear all
load Projekt_specs.mat;

%%
feat=3.28084;     %%Umrechnungsfactor
ln=4*feat;        %%länge fan cowling
lh=2.3*feat;      %%länge zur größten Stelle an der urbine
Dn=3.3*feat;      %%größter Durchmesser Engine
Dh=3*feat;        %%Durchmesser Eingang Engine
Gesamtlaenge = 4.49*feat;
beta_c = ln/Gesamtlaenge;   %%Gesamtlänge des Triebwerkes

A_cowling = ln*Dn*(2+0.35*beta_c*(Dh/Dn)+1.15*(1-beta_c));

%%
%%Gas Generator section
%%

lg=0.3*feat;     %%Länge section
Dg=1.3*feat;     %%größter Durchmesser section
Deg=0.9*feat;    %%KLeinster Durchmesser section

A_generator = pi*lg*Dg*(1-((1/3)*(1-(Deg/Dg)))*1-0.18*((Dg/lg)^(5/3)));



%%
%%Plug
%%
lp=0.19*feat;    %%Länge section
Dp=0.5*feat;     %%Anfangs Durchmesser

A_plug=0.7*pi*lp*Dp;


%%
%%Alles zsm

A_Turbine=(A_plug+A_cowling+A_plug);

%%
%%Gewicht berechnen nach:
%%Commercinal airplane design principals
%%  Refined Weight and Balance Estimate  s.317

W_e= 2.20462 * specs.m_TW;       %%Gewicht Triebwerk
S_p=10.3*feat^2;              %%umspülte fläche des pilon
n_ult=3.2;                    %%Danach soll itteriert werden sagt jasper


W_n=35.45*2*((2.33*(1.1*W_e)*A_Turbine)/(10000))^(0.59);    %%Ausrechnen in pounds
W_p=24.11*2*S_p^(0.381)*((1.46*n_ult*1.1*W_e*ln*Dn)/(10^6*cos(deg2rad(80)))^(0.952));

Gewicht=W_n/2.20462
Gewicht=W_p/2.20462

