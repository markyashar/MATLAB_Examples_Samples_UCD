% Mark Yashar
% Modified top level sample Matlab program for Phys 262
% Homework #1 ; Problem #1
% Sample program has been modified so that all the constants
% (in the "par" structure) are double their original values
% Corresponding plot is generated

% the data structure "par" contains all the fixed paramters
% It is reccomended that we use data structures so that all paramters can easily be
% passed to functions.
par.m = 2*1; % Note: If you do not put a semicolon at the end of a line the result will appear on the command line
par.b = 2*4;  

par.a = 2*3;
par.c = 2*-1;

par.A = 2*1; %matlab is case sensitive
par.f = 2*2*pi; % Note that pi is predifined in Matlab;

x = linspace(0,1) % define a vector x to be the x coordinates
y1 = f1(x,par); % define vector of y coordinates
y2 = f2(x,par); % define another vector of y coordinates
y3 = fSin(x,par); % define another vector of y coordinates

plot(x,y1,'k-');  % plots y1(x) curve as a black solid line
hold on % this command causes the next plot to appear in the same figure
plot(x,y2,'k:'); % plots y2(x) curve as a black dotted line
plot(x,y3,'k-.'); % plots y3(x) curve as a black dash-dot curve
hold off

xlim([0,1]); % set axis limits
ylim([-3, 16]);

xlabel('x axis')
ylabel('y axis')
title('Plot for Homework problem #1')
legend('linear','qadratic','sine')
