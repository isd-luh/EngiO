% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Ackley(x)
% ACKLEY is a single-objective, unconstrained, 2-dimensional test function.
%   Global minimum is at f(0,0)=0.
%   [-35 -35] <= x <= [35 35]
%
% LITERATURE:
%   D. H. Ackley (1987) Stochastic iterated genetic hillclimbing (Diss.)
%
%   M. Jamil and X.-S. Yang (2013) A Literature Survey of Benchmark
%   Functions For Global Optimization Problems, Int. Journal of
%   Mathematical Modelling and Numerical Optimisation

if isempty(x)
    y=1;
    return
end

y=-20*exp(-0.2*sqrt(0.5*(x(1)^2+x(2)^2)))-exp(0.5*(cos(2*pi*x(1))+cos(2*pi*x(2))))+20+exp(1);
end