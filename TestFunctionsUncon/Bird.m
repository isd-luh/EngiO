% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Bird(x)
% BIRD is a single-objective, unconstrained, 2-dimensional test function.
%   Global minimum is at f(4.70104,3.15294)=-106.764537.
%   [-2pi -2pi] <= x <= [2pi 2pi]
%
% LITERATURE:
%	S. K. Mishra (2006) Global Optimization by Differential Evolution 
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

y = sin(x(1)).*exp((1-cos(x(2))).^2)+cos(x(2)).*exp((1-sin(x(1))).^2)+(x(1)-x(2)).^2;
end