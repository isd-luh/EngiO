% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function res = Simionescu(x)
% SIMIONESCU is a single-objective, constrained, 2-dimensional test function.
%   Global minima are at
%   f(-/+0.84852813,-/+0.84852813)=-0.072.
%   [-1.25 -1.25] <= x <= [1.25 ... 1.25]
%
% LITERATURE:
%   P.A. Simionescu (2014) Computer Aided Graphing and Simulation Tools for
%   AutoCAD Users

if isempty(x)
    res={1, 1, 0};
    return
end

y = 0.1 * x(1) * x(2);
g = (1.0 + 0.2 * cos(8.0 * atan2(x(1),x(2) ) ) )^2 - (x(1)^2+x(2)^2);

res={y, g};
end

