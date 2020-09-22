% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Viennet2(x)
% VIENNET2 is a tree-objective, unconstrained, 2-dimensional test function.
%   Objective function on p. 256
%   [-2 -2] <= x <= [2 2]
%
% LITERATURE:
%   R. Viennet, C.Fonteix, I. Marc (1996) Multicriteria optimization using
%   a genetic algorithm for determining a Pareto set

if isempty(x)
    y=3;
    return
end

y(1) = x(1)^2 + (x(2) - 1)^2;
y(2) = x(1)^2 + (x(2)+1)^2 + 1;
y(3) = (x(1)-1)^2 + x(2)^2 + 2;

end
