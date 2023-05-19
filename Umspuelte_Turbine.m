%%
%%   Quelle: E. Torenbeek: "Synthesis of Subsonic Airplane Design", Appendix B S446-449
%%Umströmte Fläche für EEEEEEIIIIINNNNNEEEEEEE Gondel

clc
clear all
%%
feat=3.28084;      %%Umrechnungsfactor
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

A_Turbine=(A_plug+A_cowling+A_plug)/(feat^2)