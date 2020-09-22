% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Camel6(x)
% CAMEL6 is a single-objective, unconstrained, 2-dimensional test function.
%   Global minima are at f(-0.0898,0.7126)=-1.0316 and
%   f(0.0898,-0.7126)=-1.0316. 
%   [-5 -5] <= x <= [5 5]
%
% LITERATURE:
%   M. Molga and C. Smutnicki (2005) Test functions for optimization needs,
%   <http://new.zsd.iiar.pwr.wroc.pl/files/docs/functions.pdf>
%
%   M. Jamil and X.-S. Yang (2013) A Literature Survey of Benchmark
%   Functions For Global Optimization Problems, Int. Journal of
%   Mathematical Modelling and Numerical Optimisation

if isempty(x)
    y=1;
    return
end

y =(4-2.1*x(1)^2+(x(1)^4)/3)*x(1)^2+x(1)*x(2)+(4*x(2)^2-4)*x(2)^2;
end