% Homework #2, problem #4, Physics 262

%Define parameters
par.rhoc0 =  4.46e3; % present critical density   
par.rhom0 = 0.28 * par.rhoc0; % present mass density
par.rhor0 = 0.260; % present radiation density  
par.rhoL0 = 0.72 * par.rhoc0; % present cosmological constant density

% Plotting Omega(a) 
a = logspace(-6,0);
y1=omegar(a,par);
y2=omegam(a,par);
y3=omegal(a,par);
semilogx(a,y1,'k-');
hold on %
semilogx(a,y2,'k:');
semilogx(a,y3,'k-.');
hold off
ylim([-0.1,1]);
xlabel('log(a)')
ylabel('\Omega')
title('Semi-Log plot for Problem #4 of Hw #2 in Phys 262')
legend('\Omega_r','\Omega_m','\Omega_{\Lambda}')




