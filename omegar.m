function y = omegar(a,par)
  y = par.rhor0*a.^(-4)./(par.rhor0*a.^(-4) + par.rhom0*a.^(-3) + par.rhoL0*a.^(0));
