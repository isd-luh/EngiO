% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = ZDT1(x)
% ZDT1 is a two-objective, unconstrained, 30-dimensional test function.
%   The Pareto optimal front is formed with g(x)=1.
%   [0 ... 0] <= x <= [1 ... 1]
%
% LITERATURE:
%   E. Zitzler, K. Deb and L. Thiele (2000) Comparison of multiobjective
%   evolutionary algorithms: empirical results, Evolutionary computation 

if isempty(x)
    y=2;
    return
end

y(1)=x(1);
g=1+9*sum(x(2:length(x)))/(length(x)-1);
y(2)=g*(1-sqrt(x(1)/g));
end