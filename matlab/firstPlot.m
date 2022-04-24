figure 

flow = xlsread("resflows.xlsx")

subplot(2, 2, 1)
t=0:0.1:2*pi;
plot(sin(t), '.b')
title("cool plot")
plot(flow(:, 2))

subplot(2, 2, 3)
plot(flow(:, 3))

subplot(2, 2, 2)
boxplot(log(flow(:, 2)));