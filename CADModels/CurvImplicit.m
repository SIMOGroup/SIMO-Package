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
