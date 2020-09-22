% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Schaffer1(x)
% SCHAFFER1 is a two-objective, unconstrained, 1-dimensional test function.
%   [-10] <= x <= [10]
%
% LITERATURE:
%   J. Schaffer (1985) Multiple Objective Optimization with Vector 
%   Evaluated Genetic Algorithms, Proceedings of the First Int. Conference
%   on Genetic Algortihms

if isempty(x)
    y=2;
    return
end

y = [0 0];
y(1) = x.^2;
y(2) = (x-2).^2;

end