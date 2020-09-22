% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Rastrigin(x)
% RASTRIGIN is a single-objective, unconstrained, i-dimensional test function.
%   Global minimum for i=... is at f(0,...,0)=0.
%   [-5.15 ... -5.15] <= x <= [5.15 ... 5.15]
%
% LITERATURE:
%   L. A. Rastrigin (1974) Systems of extremal control

if isempty(x)
    y=1;
    return
end

nDims=numel(x);
sum=0;
A=10;

for i=1:nDims
    sum=sum+(x(i)^2-10*cos(2*pi*x(i)));
end

y = A*nDims+sum;
end