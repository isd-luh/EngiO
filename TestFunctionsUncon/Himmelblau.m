% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Himmelblau(x)
% HIMMELBLAU is a single-objective, unconstrained, 2-dimensional test function.
%   Global minima are at f(3,2)=0, f(-2.805118,3.131312)=0, f(-3.779310,-3.283186)=0
%	and f(3.584428,-1.848126)=0.
%   [-5 -5] <= x <= [5 5]
%
% LITERATURE:
%   D. M. Himmelblau (1972) Applied Nonlinear Programming, McGraw-Hill
%
%   M. Jamil and X.-S. Yang (2013) A Literature Survey of Benchmark
%   Functions For Global Optimization Problems, Int. Journal of
%   Mathematical Modelling and Numerical Optimisation

if isempty(x)
    y=1;
    return
end

y=(x(1)^2+x(2)-11)^2+(x(1)+x(2)^2-7)^2;
end