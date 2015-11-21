chisquared2d= Chisq2(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% chisquared3d= Chisq3(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
cl_abcfree = gammainc(chisquared3d/2,3/2);
cl_abfree = gammainc(chisquared2d/2,1);

figure(16)
plot(chisquared2d,chisquared3d);
xlim([-1.0 max(chisquared2d)]);
ylim([0 max(chisquared3d)]);
xlabel('\chi^2 values for case of two free parameters, a, b'); ylabel('\chi^2 values for case of three free parameters, a, b, c')
title(['Graph of \chi^2 (a,b,c, free parameters) vs.\chi^2 (a,b, free parameters for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_chi2abcfree_chi2abfree;

figure(17)
plot(cl_abfree,cl_abcfree);
xlim([-1.1 1.1]);
ylim([0 1.1]);
xlabel('C.L. values for case of two free parameters, a, b'); ylabel('C.L. values for case of three free parameters, a, b, c')
title(['Graph of C.L. (a,b,c, free parameters) vs. C.L.(a,b, free parameters for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_CLabcfree_CLabfree;






