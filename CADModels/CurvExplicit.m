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

close all
clear
clc

x = -2 : 0.01 : 2;
y = 2 * x.^2 - 3;

figure
hold on
set(gcf,'color','white')
set(gca,'XTick', -4:1:4)
set(gca,'YTick', -5:1:6)
axis([-4.4 4.4 -5.4 6.4])
daspect([1 1 1])
plot(x, y, '', 'LineWidth', 1.5);


