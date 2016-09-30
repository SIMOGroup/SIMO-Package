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

t = 0 : pi / 50 : 20 * pi;
x = exp(-0.05 * t).*sin(t);
y = exp(-0.05 * t).*cos(t);

figure
hold on
set(gcf, 'color', 'white')
set(gca, 'XTick', -1:1:1)
set(gca, 'YTick', -1:1:1)
axis([-1.4 1.4 -1.4 1.4])
daspect([1 1 1])

plot(x,y, 'LineWidth', 1.5);

