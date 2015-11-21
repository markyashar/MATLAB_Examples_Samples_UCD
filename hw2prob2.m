% Homework #2, Problem #2, Physics 262

%Define parameters
par.rhoc0 =  4.46e3; % present critical density   
par.rhom0 = 0.28 * par.rhoc0; % present mass density
par.rhor0 = 0.260; % present radiation density  
par.rhoL0 = 0.72 * par.rhoc0; % present cosmological constant density

% Plotting rho(a) (problem 2.2)
% note that the various rho_i functions are defined in rhe rhom.m, rhor.m, and
% rhol.m files (attached)
a = logspace(-6,0); % creates "a" values for the x axis
y1=rhor(a,par);
y2=rhom(a,par);
y3=rhol(a,par);
loglog(a,y1,'k-');
hold on
loglog(a,y2,'k:');
loglog(a,y3,'k-.');
hold off
%plot values
xlabel('log(a)')
ylabel('log(\rho )')
title('Log-Log Plot for Problem #2 of Hw #2 in Phys 262')
legend('\rho_r','\rho_m','\rho_{\Lambda}')




