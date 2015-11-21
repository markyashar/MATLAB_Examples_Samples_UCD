
% Homework problem #2 -- linear plot
par.g = 1; % This is A_2
par.h = 1; % This is A_3
par.i = 1; % This is A_4  

x = linspace(0.1,1); % define a vector x to be the x coordinates (ranging from 0.1 to 1)
y2 = g2(x,par); % define vector of y coordinates
y3 = g3(x,par); % define another vector of y coordinates
y4 = g4(x,par); % define another vector of y coordinates

plot(x,y2,'k-');  % plots curve for g_2(x) as a black solid line
hold on % this command causes the next plot to appear in the same figure
plot(x,y3,'k:'); % plots curve for g_3(x) as a black dotted line 
plot(x,y4,'k-.'); % plots curve for g_4(x_ as a dashdot line
hold off

xlim([0,1]); % set axis limits
ylim([0, 100]);

xlabel('x')
ylabel('g(x)')
title('Linear Axes Plot for Problem #2 of Hw #1 in Phys 262')
legend('g_2(x)','g_3(x)','g_4(x)')