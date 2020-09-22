% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Poloni(x)
% POLONI is a two-objective, unconstrained, 2-dimensional test function.
%   [-pi -pi] <= x <= [pi pi]
%
% LITERATURE:
%   C. Poloni (1995) Hybrid GA for multi objective aerodynamic shape
%   optimization, Genetic Algorithms in Engineering and Computer Science

if isempty(x)
    y=2;
    return
end

A_1 = 0.5 * sin(1) - 2 * cos(1) + sin(2) - 1.5 * cos(2);
A_2 = 1.5 * sin(1) - cos(1) + 2 * sin(2) - 0.5 * cos(2);
B_1 = 0.5 * sin(x(1)) - 2 * cos(x(1)) + sin(x(2)) - 1.5 * cos(x(2));
B_2 = 1.5 * sin(x(1)) - cos(x(1)) + 2 * sin(x(2)) - 0.5 * cos(x(2));

y(1) = 1 + (A_1 - B_1)^2 + (A_2 - B_2)^2;
y(2) = (x(1) + 3)^2 + (x(2) + 1)^2;

end