% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Beale(x)
% BEALE is a single-objective, unconstrained, 2-dimensional test function.
%   Global minimum is at f(3,0.5)=0.
%   [-4.5 -4.5] <= x <= [4.5 4.5]
%
% LITERATURE:
%   M. Jamil and X.-S. Yang (2013) A Literature Survey of Benchmark
%   Functions For Global Optimization Problems, Int. Journal of
%   Mathematical Modelling and Numerical Optimisation

if isempty(x)
    y=1;
    return
end
y = (1.5-x(1)+x(1)*x(2)).^2+(2.25-x(1)+x(1)*x(2).^2).^2+(2.625-x(1)+x(1)*x(2).^3).^2;
end