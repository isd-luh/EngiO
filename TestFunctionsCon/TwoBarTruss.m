% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function res = TwoBarTruss(x)
% TWOBARTRUSS is a two-objective, constrained, 3-dimensional test function.
%   The design variables specify vertical distance and two vertical
%   distances.
%   [1 0 0] <= x <= [3 0.01 0.01]
%
% LITERATURE:
%   K. Deb, a. Pratap and S. Moitra (2000) Mechanical Component Design for
%   Multiple Ojectives Using Elitist Non-dominated Sorting GA, PPSN

if isempty(x)
    res={2, 1, 0};
    return
end

y  = x(1);
x1 = x(2);
x2 = x(3);

sig_ac = 20 * sqrt(16 + y^2) / (y * x1);
sig_bc = 80 * sqrt( 1 + y^2) / (y * x2);

y(1) = x1 * sqrt(16 + y^2) + x2 * sqrt(1+y^2);
y(2) = max(sig_ac, sig_bc);

% g>0 is feasible  
g(1) = 1e5-max(sig_ac, sig_bc)

res={y, g}; 
end