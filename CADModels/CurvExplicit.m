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


