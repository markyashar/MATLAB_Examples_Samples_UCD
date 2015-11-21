% Homework #3, Problem #4, Physics 262

%Define parameters
par.rhoc0 =  4.46e3; % present critical density (units: Mev/m^3)   
par.rhom0 = 0.28 * par.rhoc0; % present matter density 
par.rhor0 = 0.260; % present radiation density  (Mev/m^3)
par.rhoL0 = 0.72 * par.rhoc0; % present cosmological constant density
par.rhoc0100 = 1.053e4; %  value critical density would have if H = 100 km/s/Mpc (units: Mev/m^3)
par.rhom0100 = 0.28 * par.rhoc0100; % this is rho_{m,0} as in HW #2 except the critical density
% is being calculated as if H = 100 km/s/Mpc
par.rhor0100 = 0.260; % this is the present radiation density rho_{lambda,0} as in HW #2 except the critical density
% is being calculated as if H = 100 km/s/Mpc
par.rhoL0100 = 0.72 * par.rhoc0100; % rho_{lambda,0} as in HW #2 except the critical density
% is being calculated as if H = 100 km/s/Mpc

% Plotting rho_i_100_(a) (problem 3.?) (using value of critical density as if H = 100 km/s/Mpc 
% Note that the various rho_i_100(a) functions are defined in rhe rhom100.m, rhor100.m, and
% rhol100.m files (attached)
a = logspace(-6,0); % creates "a" values for the x axis
y1=rhor100(a,par);
y2=rhom100(a,par);
y3=rhol100(a,par);
loglog(a,y1,'k-');
hold on
loglog(a,y2,'k:');
loglog(a,y3,'k-.');
hold off
%plot values
xlabel('log(a)')
ylabel('log(\rho(a) )')
title('Log-Log Plot for Problem #? of Hw #3 in Phys 262 using critical denstiy evaluated at H = 100km/s/Mpc')
legend('\rho_r_{(100)}','\rho_m_{(100)}','\rho_{\Lambda}_{(100)}')




