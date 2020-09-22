% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = FourBarTruss(x)
% FOURBARTRUSS is a two-objective, unconstrained, 1-dimensional test function.
%   [1, sqrt(2), sqrt(2), 1] <= x <= [3, 3, 3, 3]
%
% LITERATURE:
%   F Cheng (1999) GENERALIZED CENTER METHOD FOR MULTIOBJECTIVE ENGINEERING
%   OPTIMIZATION, Engineering Optimization

% sigma=10 has been factored out...

if isempty(x)
    y=2;
    return
end

F=10;
E=2e5;
L=200;

y = [0 0];
y(1) = L * ( 2 * x(1) +sqrt(2) * x(2) + sqrt(x(3)) + x(4));
y(2) = F*L/E *( 2/x(1) + 2 * sqrt(2)/x(2) - 2*sqrt(2)/x(3) + 2/x(4));

end