% Gruppe nm177 
% 369337 SINGH Jalaj
% 414809 PETERSCHMITT Ariane Clara
% 468310 DOROZALLA Jill
% 468311 BLEI Saskia Maxi

function [T,S] = integral(f,a,b,n)

% Ziel : Integral von a nach b von f(x)dx berechnen

% Eingabewerte: Funktion f, a untere Intervallsgrenze, b obere
% Intervallsgrenze, n Anzahl der Stützstellen
% Intervalllänge h : h=(b-a)/(n-1)
% 1. Rückgabewert T : summierte Trapezenregel
% 2. Rückgabewert S : summierte Simpson-Regel

x=linspace(a,b,n);
h=(b-a)/(n-1);
x(1)=a;
for j=2:n
    x(j)=a+(j-1)*h;
end
x(n)=b;

summeT=0;
for j=2:n-1
    summeT=summeT+f(x(j));
end

T=h*(1/2*(f(a)+f(b))+summeT);

summeS=0;
for j=1:n-1
    summeS=summeS+h/6*((f(x(j))+4*f((x(j)+x(j+1))/2)+f(x(j+1))));
end

S=summeS;

% T in Zeile 30 & S in Zeile 37 mit Semikolon versehen
% sonst T & S doppelt ausgegeben : von Zeilen 30, 37 + [T,S] Zeile 7

% Test für die Verbesserung
% f=@(x)sin(x)
% a=-1
% b=1
% n=100
% [T,S] = integral(f,a,b,n)
% das zu erwarten Ergebnis : T=S≈0 