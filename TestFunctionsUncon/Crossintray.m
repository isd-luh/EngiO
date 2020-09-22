% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Crossintray(x)
% CROSSINTRAY is a single-objective, unconstrained, 2-dimensional test function.
%   Global minima are at
%   f(-/+1.349406685353340,-/+1.349406608602084)=-2.06261218.
%   [-10 -10] <= x <= [10 10]
%
% LITERATURE:
%   S. K. Mishra (2006) Global Optimization by Differential Evolution 
%	and Particle Swarm Methods: Evaluation on Some Benchmark Functions,
%	SSRN Electronic Journal
%
%   M. Jamil and X.-S. Yang (2013) A Literature Survey of Benchmark
%   Functions For Global Optimization Problems, Int. Journal of
%   Mathematical Modelling and Numerical Optimisation

if isempty(x)
    y=1;
    return
end

y =-0.0001*(abs(sin(x(1))*sin(x(2))*exp(abs(100.0-sqrt(x(1)^2+x(2)^2)/pi)))+1.0)^0.1;
end

