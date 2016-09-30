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

figure
hold on
set(gcf,'color','white')
set(gca,'XTick', -1:1:1)
set(gca,'YTick', -1:1:1)
axis([-1.4 1.4 -1.4 1.4])
daspect([1 1 1])
h = ezplot('2*x^2+y^2/1.5-1');
title([]);
set(h,'LineWidth', 1.5, 'Color', 'b');
