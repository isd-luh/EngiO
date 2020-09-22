% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Binh(x)
% BINH and Korn is a two-objective, unconstrained, 2-dimensional test function.
%   [0 0] <= x <= [5 3]
%
% LITERATURE:
%   T.T. Binh and U. Korn (1997) MOBES: A Multiobjective Evolution Strategy
%   for Constrained Optimization Problems, 3rd International Conference on
%   Genetic Algorithms

if isempty(x)
    y=2;
    return
end

y(1)=4*(x(1)^2+x(2)^2);
y(2)=(x(1)-5)^2+(x(2)-5)^2;
end