% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function res = SpringDesign(x)
% SPRINGDESIGN is a two-objective, constrained, 3-dimensional test function.
%   Design variables are of type integer and float. They specify no. of
%   spring coils, (discrete) wire diameter and coil diameter.
%   [3 1 1.0] <= x <= [20 42 3.0]
%
% LITERATURE:
%   K. Deb, a. Pratap and S. Moitra (2000) Mechanical Component Design for
%   Multiple Ojectives Using Elitist Non-dominated Sorting GA, PPSN

if isempty(x)
    res={2, 8, 0};
    return
end

diameters = [0.0090 0.0095 0.0104 0.0118 0.0128 0.0132 0.0140 0.0150 0.0162 0.0173 0.0180 0.0200 0.0230 0.0250 0.0280 0.0320 0.0350 0.0410 0.0470 0.0540 0.0630 0.0720 0.0800 0.0920 0.1050 0.1200 0.1350 0.1480 0.1620 0.1770 0.1920 0.2070 0.2250 0.2440 0.2630 0.2830 0.3070 0.3310 0.3620 0.3940 0.4375 0.5000];

x(2) = diameters(x(2));

G = 11.5e6;
S = 189000;
P = 300;
P_max = 1000;
l_max = 14.0;
d_min = 0.2;
D_max = 3.0;
delta_w = 1.25;
delta_pm = 6.0;
V_max = 30;

C = x(3) / x(2);
K = (4 * C-1) / (4*C-4) + 0.615*x(2)/x(3);
k = G*x(2)^4/(8*x(1)*x(3)^3);
delta_p = P/k;

y(1) = 0.25 * pi^2 * x(2)^2 * x(3) * (x(1) + 2);
y(2) = 8.0 * K * P_max * x(3) / (pi * x(2)^3);

% g>0 is feasible 
g(1) = l_max - P_max / k;
g(2) = x(2) - d_min;
g(3) = D_max - (x(2) + x(3));
g(4) = C-3.0;
g(5) = delta_pm - delta_p;
g(6) = (P_max - P) / k - delta_w;
g(7) = S - 8.0 * K * P_max * x(3) / (pi * x(2)^3);
g(8) = V_max - 0.25 * pi^2 * x(2)^2 * x(3) * (x(1) + 2);

res={y, g, []}; 
end