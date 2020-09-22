% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Rosenbrock(x)
% ROSENBROCK is a single-objective, unconstrained, i-dimensional test function.
%   Global minimum for i=... is at f(1,...,1)=0.
%   [-30 ... -30] <= x <= [30 ... 30]
%
% LITERATURE:
%   H. H. Rosenbrock (1960) An automatic method for finding the greatest or
%   least value of a function, The Computer Journal
%
%   M. Jamil and X.-S. Yang (2013) A Literature Survey of Benchmark
%   Functions For Global Optimization Problems, Int. Journal of
%   Mathematical Modelling and Numerical Optimisation

if isempty(x)
    y=1;
    return
end

y=0;
nDims=numel(x);

for i=1:nDims-1
    y=y+100*(x(i+1)-x(i)^2)^2+(1-x(i))^2;
end

end