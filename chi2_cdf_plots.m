x1 = 1.0:25.5; v = 1:15;
[X1,V] = meshgrid(x1,v);
probability = chi2cdf(X1(:),V(:))
probability = reshape(probability,length(v),length(x1));
surf(x1,v,probability)
xlabel('chi^{2}'); ylabel('degrees of freedom'); zlabel('Chi^{2} cumulative distribution function (cdf)')
% print -deps plot_chi2_cdf_dof

figure 

x1 = linspace(0.0,25.5,200);
prob1 = chi2cdf(x1,1);
prob2 = chi2cdf(x1,2);
prob3 = chi2cdf(x1,3);
prob4 = chi2cdf(x1,4);
prob5 = chi2cdf(x1,5);
plot(x1,prob1,':*', x1,prob2,'b-',x1,prob3,'b--', x1, prob4,'b:+',x1, prob5,':ks')
grid on
legend('1 d.o.f.','2 d.o.f.','3 d.o.f.','4 d.o.f.','5 d.o.f.',0);
% probability = reshape(probability,length(v),length(x1));
% plot3(x1,v,probability)
% plot(x1,probability)
xlabel('chi^{2}'); ylabel('Chi^{2} cumulative distribution function (cdf)')
% print -deps plot2D_chi2_cdf_dof


% >> as = 1.0:0.005:1.2;
% >> bs = 1.4:0.005:1.6;
% >> cs = 1.7:0.005:1.9;
% >> [X1,X2] = meshgrid(as,bs);
% >> f = subs(chisq);
% >> f = reshape(f,length(bs),length(as));