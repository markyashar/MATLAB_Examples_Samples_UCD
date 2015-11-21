clear all
% Main program for calculating different chi^2 values, etc...

a = 1.1000;
b = 1.5000;
c = 5.8000;
% sigmaa = 0.001:0.001:0.1; % standard deviations of parameters chosen by me
% sigmab = 0.001:0.001:0.1;
% sigmac = 0.001:0.001:0.1;
sigmax  = 0.0030;
sigmay  = 0.3000;
sigmaz  = 0.5000;
as = 1.0975:0.00005:1.1025;
bs = 1.4975:0.00005:1.5025;
cs = 5.7975:0.00005:5.8025;
chisquared3d= Chisq3(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% [X1,X2] = meshgrid(as,bs);
% chisquare3d = chisquared3d(X1(:),X2(:));
% chisquare3d=reshape(chisquare3d,length(bs),length(as));
plot3(as,bs,chisquared3d)
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of three free parameters, a, b, c');
grid on
title({'Graph of \chi^2 vs. a_{data} and b_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
    '3 free parameters and 3 observables'});
print -deps plot2_chi2_3dof;

figure(2)
plot3(as,cs,chisquared3d)
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of three free parameters, a, b, c');
grid on
title({'Graph of \chi^2 vs. a_{data} and c_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
    '3 free parameters and 3 observables'});
print -deps plot3_a_c_chi2_3dof;


cl = gammainc((chisquared3d)/2,3/2);
figure(3)
plot3(as,bs,cl);
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Confidence level values for case of three free parameters, a, b, c')
title({'Graph of confidence level vs. a_{data} and b_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
   '3 free parameters and 3 observables'});
grid on
print -deps plot2_cl_3dof;

figure(4)
plot(chisquared3d,cl);
xlabel('Chi^{2} values for case of three free parameters, a, b, c'); ylabel('Confidence level values for case of three free parameters, a, b, c')
title({'Graph of confidence level vs.\chi^2 for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
    '3 free parameters and 3 observables'});
grid on
print -deps plot_cl_chi2_3dof;

figure(5)
plot3(as,cs,cl);
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of three free parameters, a, b, c')
title({'Graph of confidence level vs. a_{data} and c_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
   '3 free parameters and 3 observables'});
grid on
print -deps plot3_a_c_cl_3dof;

figure(6)
plot3(bs,cs,chisquared3d)
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of three free parameters, a, b, c');
grid on
title({'Graph of \chi^2 vs. b_{data} and c_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
    '3 free parameters and 3 observables'});
print -deps plot4_b_c_chi2_3dof;

figure(7)
plot3(bs,cs,cl);
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of three free parameters, a, b, c')
title({'Graph of confidence level vs. b_{data} and c_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
   '3 free parameters and 3 observables'});
grid on
print -deps plot4_b_c_cl_3dof;


chisquared2d= Chisq2(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
figure(8)
plot3(as,bs,chisquared2d)
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of two free parameters, a, b');
grid on
title({'Graph of \chi^2 vs. a_{data} and b_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
    '2 free parameters and 3 observables'});
print -deps plot2_chi2_2dof;

cl = gammainc((chisquared2d)/2,2/2);
figure(9)
plot3(as,bs,cl);
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Confidence level values for case of two free parameters, a, b')
title({'Graph of confidence level vs. a_{data} and b_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
   '2 free parameters and 3 observables'});
grid on
print -deps plot2_cl_2dof;

figure(10)
plot(chisquared2d,cl);
xlabel('Chi^{2} values for case of two free parameters, a, b, with c fixed'); ylabel('Confidence level values for case of two free parameters, a, b, and c fixed')
title({'Graph of confidence level vs.\chi^2 for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
    '2 free parameters and 3 observables'});
grid on
print -deps plot_cl_chi2_2dof;


figure(11)

cl1 = gammainc((chisquared3d)/2,3/2);
cl2 = gammainc((chisquared2d)/2,2/2);
plot3(as,bs,cl1, ':*', as,bs,cl2,'b-')
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Confidence level values for case of two and three free parameters')
title({'Graph of confidence level vs. a_{data} and b_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
   '3 parameters and 3 observables'});
grid on
legend('3 free parameters (a,b,c)','2 free parameters (a,b)',0);
print -deps plot2_cl_2and3dof;

figure (12)
plot3(as,bs,chisquared2d, ':*', as,bs,chisquared3d,'b-')
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of two and three free parameters');
grid on
legend('3 free parameters (a,b,c)','2 free parameters (a,b)',0);
title({'Graph of \chi^2 vs. a_{data} and b_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
    '2 and three free parameters and 3 observables'});
print -deps plot2_chi2_2and3dof;


