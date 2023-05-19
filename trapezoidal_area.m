function [area] = trapezoidal_area(x, y)
% x: Vektor mit x-Datenpunkten
% y: Vektor mit y-Datenpunkten
% area: die berechnete Fläche

% Berechnung der Teilintervallbreite
h = x(2) - x(1);

% Berechnung der Fläche mit der Trapezregel
area = h/2 * (y(1) + 2*sum(y(2:end-1)) + y(end));

end