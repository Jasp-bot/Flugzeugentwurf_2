%%
%%   Quelle: E. Torenbeek: "Synthesis of Subsonic Airplane Design", Appendix B S446-449
%%Umströmte Fläche Turbine
%%

ln=4;        %%länge fan cowling
lh=2.3;      %%länge zur größten Stelle an der urbine
Dn=3.3;      %%größter Durchmesser Engine
Dh=3;        %%Durchmesser Eingang Engine
Gesamtlaenge= 4,49;
beta_c=ln/Gesamtlaenge;   %%Gesamtlänge des Triebwerkes

A_cowling=ln*Dn*(2+0.35*beta*(Dh/Dn)+1.15*(1-beta))

%%
%%Gas Generator section
%%

lg=0.3;     %%Länge section
Dg=1.3;     %%größter Durchmesser section
Deg=0.9;    %%KLeinster Durchmesser section
beta_g=lg/Gesamtlaenge;

A_generator=pi*lg*Dg*(1-((1/3)*(1-(Deg/Dg)))*1-0.18*((Dg/lg)^(5/3)))

%%
%%Plug
%%
lp=0.19;    %%Länge section
Dp=0.5;     %%Anfangs Durchmesser

A_plug=0.7*pi*lp*Dp