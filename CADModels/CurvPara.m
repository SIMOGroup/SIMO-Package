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

