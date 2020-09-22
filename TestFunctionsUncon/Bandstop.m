% EngiO - Engineering Optimization framework 
% Copyright C 2020 Institute of Structural Analysis, Leibniz University Hannover 
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or at your option any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

function y = Bandstop(x)
% BANDSTOP is a digital filter optimization task.
% The coefficients of a FIR filter are optimized for a low-pass bandstop
% task. Frequencies below (above) unity shall pass with unity (zero)
% amplitude.

if isempty(x)
    y=1;
    return
end

[h, w] = freqz(x, 1, 70);
thres = 1.0;

y = sum( (1.0 - abs(h(find(w < thres)))).^2) + sum(abs(h(find(w >= thres))).^2);
end