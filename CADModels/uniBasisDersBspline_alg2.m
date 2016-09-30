%{
Copyright (C) <2014-2016>  <Khanh Chau-Nguyen, Hung Nguyen-Xuan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

clear
close all
clc

KntVect = [0 0 0 1 2 3 4 5 6 7 8 9 10 10 10];
KntVect = KntVect ./ max(KntVect);
p = 2; % order of basis functions
NCtrlPts = numel(KntVect) - p - 1;

xi = linspace(0, 1, 1001);
Idx = FindSpan(NCtrlPts, p, xi, KntVect);

N = DersBasisFuns(Idx, xi, p, 1, KntVect);
N0 = squeeze(N(:, :, 1));

figure 
hold on
set(gcf, 'color', 'white');
plot(xi, N0');




