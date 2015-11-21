clear all
syms a b c as bs cs as_draw bs_draw cs_draw sigma_as sigma_bs sigma_cs sigmaa... 
    sigmab sigmac sigmax sigmay sigmaz real;
% x = a*b*c; y = a^2 + b^2 + c^2; z = a + b + c;
x = a*b*c;   % three observables here
y = a^2 + b^2 + c^2; 
z = a + b + c;
Fdata = [1/(sigmax)^2 0 0; 0 1/(sigmay)^2 0; 0 0 1/(sigmaz)^2]; % data fisher matrix
M = jacobian([x;y;z],[a,b,c]); % Jacobian matrix for transforming fisher matrix
F = M'*Fdata*M; % transformed fisher matrix
% a = 0.1; b = 5.9; c = 0.8; sigmax = 0.005; sigmay = 0.009; sigmaz = 0.003;
% a = 2.9; b = 5.9; c = 0.8; sigmax = 0.005; sigmay = 0.009; sigmaz = 0.003;
% a = 2.9 
% b = 5.9 
% c = 0.8 
% as = 2.8985
% bs = 5.8984
% cs = 0.7985
% sigmax = 0.02
% sigmay = 0.08
% sigmaz = 0.01
a = 1.16; % 'theoretical' or central values of parameters (chosen by me)'
b = 1.63;
c = 1.18;
sigmaa = 0.001; % standard deviations of parameters chosed by me
sigmab = 0.002;
sigmac = 0.003;
% as = 2.8985
% bs = 5.8984
% cs = 0.7985
% as = 2.8985
% bs = 5.8984
% cs = 0.7985
% as = 1.116
% bs = 1.538
% as = 1.110
% bs = 1.532
% cs = 1.799
% sigmax = sqrt(sigmaa^2*(b*c)^2 +sigmab^2*(a*c)^2 + sigmac^2*(a*b)^2);
% standard deviation of observables calculated via error propagation
% approximation equations
% sigmay = sqrt(sigmaa^2*(2*a)^2 + (sigmab^2)*(2*b)^2 + sigmac^2*(2*c)^2);
% sigmaz = sqrt(sigmaa^2 + sigmab^2 + sigmac^2);
sigmax = 0.004;
sigmay = 0.005;
sigmaz = 0.006;
Fdata_subs = subs(Fdata);
Fsubs = subs(F);
eigsFsubs = eig(Fsubs);
eig1 = eigsFsubs(1);
eig2 = eigsFsubs(2);
eig3 = eigsFsubs(3);
sigma1 = 1/sqrt(eigsFsubs(1));
sigma2 = 1/sqrt(eigsFsubs(2));
sigma3 = 1/sqrt(eigsFsubs(3));
Msubs = subs(M);
fom = sqrt(det(Fsubs));
xsubs = subs(x);
ysubs = subs(y);
zsubs = subs(z);
precisionx = (sigmax/xsubs)*100;
precisiony = (sigmay/ysubs)*100;
precisionz = (sigmaz/zsubs)*100;
% n = 1
% for i = 1:n
% xs= normrnd(xsubs,sigmax,[1]);
% ys= normrnd(ysubs,sigmay,[1]);
% zs= normrnd(zsubs,sigmaz,[1]);
% as = normrnd(a,sigma3,[1])
% bs = normrnd(b,sigma2,[1]);
% as = 1.0900
% bs = 1.4800
% cs = 1.7700
% cs = 1.7999
% cs=mean(normrnd(1.8,sigma2+0.015,[1 1000000]))
% cs = mean(normrnd(c,sigma3,[1 100]))
% pt = [xsubs;ysubs;zsubs];
j=3500;
for i = 1:j;
% i=1:j;    
as_draw=normrnd(a,sigmaa,[1]);
bs_draw=normrnd(b,sigmab,[1]);
cs_draw=normrnd(c,sigmac,[1]);
as=as_draw;
result_as(i) = as;
bs=bs_draw;
result_bs(i) = bs;
cs=cs_draw;
result_cs(i) = cs;
% sigma_as =(std(as_draw))/1000;
% sigma_bs = (std(bs_draw))/1000;
% sigma_cs = (std(cs_draw))/1000;
pt = [a;b;c];
% ps = [xs;ys;zs];
ps = [as;bs;cs];
chisq = (pt-ps)'*Fsubs*(pt-ps);
result_chisq(i)=chisq;
% chisq_data(i) = (pt-ps)'*Fdata_subs*(pt-ps);
% end
% xsmean = mean(xs)
% ysmean = mean(ys)
% zsmean = mean(zs)
% pt = [xsubs;ysubs;zsubs]
% psmean = [xsmean;ysmean;zsmean]
% mean_chisq = mean(chisq)
% stdev_chisq = std(chisq)
% stdev_mean_chisq = std(chisq)/sqrt(n)
% prob_dens = ((2*pi)^-(3/2))*sqrt(det(Fsubs))*exp(-mean_chisq/2)
prob_dens = ((2*pi)^-(3/2))*sqrt(det(Fsubs))*exp(-chisq/2);
result_prob_dens(i) = prob_dens;
% mean_chisq_data = mean(chisq_data)
% stdev_chisq_data = std(chisq_data)
% stdev_mean_chisq_data = std(chisq_data)/sqrt(n)
% prob_dens_data = ((2*pi)^-(3/2))*sqrt(det(Fsubs))*exp(-mean_chisq_data/2)
n=3;
cl = gammainc((chisq)/2,n/2);
result_cl(i) =cl;
chisq_red = chisq/n;
result_chisq_red(i) = chisq_red;
alpha = 1-cl;
result_alpha(i) = alpha;
cl_percent = 100*cl;
result_cl_percent(i) = cl_percent;
% alpha_percent = (1-cl)*100;
alpha_percent = 100 - cl_percent;
result_alpha_percent(i) = alpha_percent;
% end
sigma_as = (std(result_as));
sigma_bs = (std(result_bs));
sigma_cs = (std(result_cs));
std_as = (std(result_as))/j;
std_bs = (std(result_bs))/j;
std_cs = (std(result_cs))/j;
pd = mean(result_prob_dens);
con_level = mean(result_cl);
chisq = mean(result_chisq);
chisq_reduced = mean(result_chisq_red);
alpha_sig = mean(result_alpha);
conlevel_percent = mean(result_cl_percent);
alpha_sig_percent = mean(result_alpha_percent);
std_pd = (std(result_prob_dens))/j;
std_cl = (std(result_cl))/j;
sigma_chisq = std(result_chisq);
sigma_chisq_red = std(result_chisq_red);
std_chisq = (std(result_chisq))/j;
std_chisq_red = (std(result_chisq_red))/j;
std_alpha_sig = (std(result_alpha))/j;
sigma_alpha_sig = std(result_alpha);
std_cl_percent = (std(result_cl_percent))/j;
sigma_cl_percent = std(result_cl_percent);
std_alpha_percent = (std(result_alpha_percent))/j;
sigma_alpha_percent = std(result_alpha_percent);
end
% end
% fomsub12 = sqrt(eig1*eig2)
% fomsub13 = sqrt(eig1*eig3)
% fomsub23 = sqrt(eig2*eig3)
fprintf('                                         \n')
fprintf('                                         \n')
fprintf('                                         \n')
fprintf(' 3 observables, 3 free parameters \n')
fprintf('a                                      %g \n',a)
fprintf('b                                      %g \n', b)
fprintf('c                                      %g \n', c)
fprintf('sigma_1                                %g \n', sigma1)
fprintf('sigma_2                                %g \n', sigma2)
fprintf('sigma_3                                %g \n', sigma3)
fprintf('sigma_a                                %g \n', sigmaa)
fprintf('sigma_b                                %g \n', sigmab)
fprintf('sigma_c                                %g \n', sigmac)
fprintf('x                                      %g \n', xsubs)
fprintf('y                                      %g \n', ysubs)
fprintf('z                                      %g \n', zsubs)
fprintf('sigma_x                                %g \n', sigmax)
fprintf('sigma_y                                %g \n', sigmay)
fprintf('sigma_z                                %g \n', sigmaz)
fprintf('sigma_as                               %g \n', sigma_as)
fprintf('sigma_bs                               %g \n', sigma_bs)
fprintf('sigma_cs                               %g \n', sigma_cs)
fprintf('stddev_mean_as                               %g \n', std_as)
fprintf('stddev_mean_bs                               %g \n', std_bs)
fprintf('stddev_mean_cs                               %g \n', std_cs)
fprintf('probability                            %g \n', pd)
fprintf('standard dev. of mean of probability   %g \n', std_pd)
fprintf('CONFIDENCE LEVEL                      %g \n', conlevel_percent)
fprintf('standard deviation of mean of c.l.     %g \n', std_cl_percent)
fprintf('SIGMA_CL_PERCENT                      %g \n', sigma_cl_percent)
fprintf('significance (alpha)                   %g \n',alpha_sig_percent)
fprintf('standard deviation of mean of alpha    %g \n',std_alpha_percent)
fprintf('sigma_alpha_percent                    %g \n', sigma_alpha_percent)
fprintf('chi^2                                  %g \n', chisq)
fprintf('standard deviation of mean of chi^2    %g \n', std_chisq)
fprintf('sigma_chisq                  %g \n', sigma_chisq)
fprintf('chi^2_reduced                          %g \n', chisq_reduced)
fprintf('standard dev. of mean of chi^2_red %g \n',std_chisq_red)
fprintf('sigma_chisq_red                  %g \n',sigma_chisq_red)

% % a = 1.1000;
% % b = 1.5000;
% % c = 5.8000;
% % sigmaa = 0.001:0.001:0.1; % standard deviations of parameters chosen by me
% % sigmab = 0.001:0.001:0.1;
% % sigmac = 0.001:0.001:0.1;
% % sigmax  = 0.0030;
% % sigmay  = 0.3000;
% % sigmaz  = 0.5000;
% as = 1.0975:0.00005:1.1025;
as = linspace(min(result_as),max(result_as),3500);
bs = linspace(min(result_bs),max(result_bs),3500);
cs = linspace(min(result_cs),max(result_cs),3500);
% % as =min(result_as):max(result_as);
% % bs =min(result_bs):max(result_bs);
% % cs =min(result_cs):max(result_cs);
% bs = 1.4975:0.00005:1.5025;
% cs = 5.7975:0.00005:5.8025;
chisquared3d= Chisq3(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% % [X1,X2] = meshgrid(as,bs);
% % chisquare3d = chisquared3d(X1(:),X2(:));
% % chisquare3d=reshape(chisquare3d,length(bs),length(as));
figure(1)
plot3(as,bs,chisquared3d)
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of three free parameters, a, b, c');
grid on
title(['Graph of \chi^2 vs. a_{data} and b_{data} for the case of a =',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ';  3 free parameters & 3 observables']);
print -deps plot2_as_bs_chi2_3dof;

figure(2)
plot3(as,cs,chisquared3d)
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of three free parameters, a, b, c');
grid on
title(['Graph of \chi^2 vs. a_{data} and c_{data} for the case of a =',num2str(a),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ';  3 free parameters & 3 observables']);
print -deps plot3_as_cs_chi2_3dof;


cl = gammainc(chisquared3d/2,3/2);
figure(3)
plot3(as,bs,cl);
% xlim([0 max(as)]);
% ylim([0 max(bs)]);
% zlim([0 1.5]);
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Confidence level values for case of three free parameters, a, b, c')
title(['Graph of confidence level vs. a_{data} and b_{data} for the case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ';  3 free parameters & 3 observables']);
grid on
print -deps plot2_as_bs_cl_3dof;

figure(4)
plot(chisquared3d,cl);
xlim([-1000.0 max(chisquared3d)]);
ylim([0 2]);
xlabel('Chi^{2} values for case of three free parameters, a, b, c'); ylabel('Confidence level values for case of three free parameters, a, b, c')
title(['Graph of confidence level vs.\chi^2 for the case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ';  3 free parameters & 3 observables']);
grid on
print -deps plot_cl_chi2_3dof;

figure(5)
plot3(as,cs,cl);
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of three free parameters, a, b, c')
title(['Graph of confidence level vs. a_{data} and c_{data} for the case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ';  3 free parameters & 3 observables']);
grid on
print -deps plot3_as_cs_cl_3dof;

figure(6)
plot3(bs,cs,chisquared3d)
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of three free parameters, a, b, c');
grid on
title(['Graph of \chi^2 vs. b_{data} and c_{data} for the case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ';  3 free parameters & 3 observables']);
print -deps plot4_bs_cs_chi2_3dof;

figure(7)
plot3(bs,cs,cl);
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of three free parameters, a, b, c')
title(['Graph of confidence level vs. b_{data} and c_{data} for the case of a =',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ';  3 free parameters & 3 observables']);
grid on
print -deps plot4_bs_cs_cl_3dof;


% chisquared2d= Chisq2(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% figure(8)
% plot3(as,bs,chisquared2d)
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of two free parameters, a, b');
% grid on
% title({'Graph of \chi^2 vs. a_{data} and b_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
%     '2 free parameters and 3 observables'});
% print -deps plot2_chi2_2dof;
% 
% cl = gammainc((chisquared2d)/2,2/2);
% figure(9)
% plot3(as,bs,cl);
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Confidence level values for case of two free parameters, a, b')
% title({'Graph of confidence level vs. a_{data} and b_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
%    '2 free parameters and 3 observables'});
% grid on
% print -deps plot2_cl_2dof;
% 
% figure(10)
% plot(chisquared2d,cl);
% xlabel('Chi^{2} values for case of two free parameters, a, b, with c fixed'); ylabel('Confidence level values for case of two free parameters, a, b, and c fixed')
% title({'Graph of confidence level vs.\chi^2 for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
%     '2 free parameters and 3 observables'});
% grid on
% print -deps plot_cl_chi2_2dof;
% 
% 
% figure(11)
% 
% cl1 = gammainc((chisquared3d)/2,3/2);
% cl2 = gammainc((chisquared2d)/2,2/2);
% plot3(as,bs,cl1, ':*', as,bs,cl2,'b-')
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Confidence level values for case of two and three free parameters')
% title({'Graph of confidence level vs. a_{data} and b_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
%    '3 parameters and 3 observables'});
% grid on
% legend('3 free parameters (a,b,c)','2 free parameters (a,b)',0);
% print -deps plot2_cl_2and3dof;
% 
% figure (12)
% plot3(as,bs,chisquared2d, ':*', as,bs,chisquared3d,'b-')
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of two and three free parameters');
% grid on
% legend('3 free parameters (a,b,c)','2 free parameters (a,b)',0);
% title({'Graph of \chi^2 vs. a_{data} and b_{data} for the case of a = 1.1, b =1.5, c =5.8, \sigma_{x} = 0.03, \sigma_{y} = 0.3, \sigma_{z} =0.5';...
%     '2 and three free parameters and 3 observables'});
% print -deps plot2_chi2_2and3dof;
% 
% 
% 
% 
