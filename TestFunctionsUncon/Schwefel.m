% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Schwefel(x)
% SCHWEFEL is a single-objective, unconstrained, i-dimensional test function.
%   Global minimum for i=... is at f(420.9687,...,420.9687)=-418.983.
%   [-500 ... -500] <= x <= [500 ... 500]
%
% LITERATURE:
%   H.-P. Schwefel (1981) Numerical optimization of computer models, John
%   Wiley & Sons, Inc.
%
%   M. Jamil and X.-S. Yang (2013) A Literature Survey of Benchmark
%   Functions For Global Optimization Problems, Int. Journal of
%   Mathematical Modelling and Numerical Optimisation

if isempty(x)
    y=1;
    return
end

nDims=numel(x);
sum=0;

for i = 1:nDims
    sum=sum+x(i)*sin(sqrt(abs(x(i))));
end

y = -1/nDims*sum;
end