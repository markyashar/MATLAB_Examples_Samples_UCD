clear all
syms a b c as bs cs as_draw bs_draw cs_draw sigma_as sigma_bs sigma_cs sigmaa... 
    sigmab sigmac sigmax sigmay sigmaz real;
x = a*b*c; y = a^2 + b^2 + c^2; z = a + b + c;
Fdata = [1/(sigmax)^2 0 0; 0 1/(sigmay)^2 0; 0 0 1/(sigmaz)^2];
M = jacobian([x;y;z],[a,b,c]);
F = M'*Fdata*M;
% a = 2.9; b = 5.9; c = 0.8; sigmax = 0.005; sigmay = 0.009; sigmaz = 0.003;
% a = 2.9 
% b = 5.9 
% c = 0.8
% a = 21.96;
a = 1.16;
b = 1.63;
c = 1.18;
% sigmaa = 0.2;
% sigmaa = 0.00001;
sigmaa = 0.2;
sigmab = 0.002;
sigmac = 0.003;
% as = 2.8999
% bs = 5.8988
% cs = 0.7999
% as = 2.8985
% bs = 5.8984
% cs = 0.7985
% as = 1.110
% bs = 1.532
% as = 1.0998
% bs = 1.4999
% cs = 1.7999
% as=mean(normrnd(1.1,sigmarem3_3,[1 100]))
% bs=mean(normrnd(1.5,sigmarem3_2,[1 100]))
% cs=mean(normrnd(1.8,sigmarem3_1,[1 100]))
% sigmax = sqrt(sigmaa^2*(b*c)^2 +sigmab^2*(a*c)^2 + sigmac^2*(a*b)^2);
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
Fsubsrem3 = Fsubs(1:2,1:2);
fom_rem3 = sqrt(det(Fsubsrem3));
eigsFsubsrem3 = eig(Fsubsrem3);
eigrem3_1 = eigsFsubsrem3(1);
eigrem3_2 = eigsFsubsrem3(2);
% eigrem3_3 = eigsFsubsrem3(3)
eigrem3_3 = eigsFsubs(3);
sigmarem3_1 = 1/sqrt(eigsFsubsrem3(1));
sigmarem3_2 = 1/sqrt(eigsFsubsrem3(2));
sigmarem3_3 = 1/sqrt(eigsFsubs(3));
% as=normrnd(a,sigmaa,[1])
% bs=normrnd(b,sigmab,[1])
% cs=normrnd(c,sigmac,[1])
% as = 1.0900
% bs = 1.4800
% cs = 1.7700
% n = 5000
% for i = 1:n
% xs= normrnd(xsubs,sigmax,[1]);
% ys= normrnd(ysubs,sigmay,[1]);
% zs= normrnd(zsubs,sigmaz,[1]);
% as = normrnd(a,sigmarem3_2,[1]);
% bs = normrnd(b,sigmarem3_1,[1]);
% cs = normrnd(c,sigma1,[1]);
% pt = [xsubs;ysubs;zsubs];
% j=3500;
j=100;
for i = 1:j;
% as_draw=normrnd(a,sigmaa,[1]);
% bs_draw=normrnd(b,sigmab,[1]);
% cs_draw=normrnd(c,sigmac,[1]);
% as=as_draw;
as=normrnd(a,sigmaa,[1]);
bs=normrnd(b,sigmab,[1]);
cs=normrnd(c,sigmac,[1]);
result_as(i) = as;
% bs=bs_draw;
result_bs(i) = bs;
% cs=cs_draw;
result_cs(i) = cs;
% sigma_as =(std(as_draw))/1000
% sigma_bs = (std(bs_draw))/1000
% sigma_cs = (std(cs_draw))/1000
pt = [a;b];
ptall =[a;b;c];
% ps = [xs;ys;zs];
% ps = [as;bs];
ps = [as;bs];
% psall =  [as;bs;cs];
psall = [as;bs;cs];
chisq_rem3 = (pt-ps)'*Fsubsrem3*(pt-ps);
chisq = (ptall-psall)'*Fsubs*(ptall-psall);
result_chisq_rem3(i)=chisq_rem3;
result_chisq(i)=chisq;
% end
% xsmean = mean(xs)
% ysmean = mean(ys)
% zsmean = mean(zs)
% pt = [xsubs;ysubs;zsubs]
% psmean = [xsmean;ysmean;zsmean]
% mean_chisq_rem3 = mean(chisq_rem3)
% stdev_chisq_rem3 = std(chisq_rem3)
% stdev_mean_chisq_rem3 = std(chisq_rem3)/sqrt(n)
% prob_dens = ((2*pi)^-(2/2))*sqrt(det(Fsubsrem3))*exp(-mean_chisq_rem3/2)
prob_dens = ((2*pi)^-(2/2))*sqrt(det(Fsubsrem3))*exp(-(chisq_rem3)/2);
prob_dens_all = ((2*pi)^-(3/2))*sqrt(det(Fsubs))*exp(-chisq/2);
result_prob_dens(i) = prob_dens;
result_prob_dens_all(i) = prob_dens_all;
% n=2;
cl = gammainc(chisq_rem3/2,1);
cl_all = gammainc((chisq)/2,3/2);
result_cl(i) =cl;
result_cl_all(i) = cl_all;
chisq_rem3_red = chisq_rem3/2;
chisq_red(i) = chisq/3;
result_chisq_rem3_red(i) = chisq_rem3_red;
result_chisq_red = chisq_red;
% alpha(i) = 1-cl;
% alpha_all(i) = 1 - cl_all;
result_alpha(i) = 1-cl;
result_alpha_all(i) = 1-cl_all;
%cl_percent(i) = 100*cl;
%cl_percent_all(i) = 100*cl_all;
result_cl_percent(i) = 100*cl;
result_cl_percent_all(i) = 100*cl_all;
result_alpha_percent(i) = (1-cl)*100;
result_alpha_percent_all(i) = (1-cl_all)*100;
% result_alpha_percent(i) = alpha_percent;
% result_alpha_percent_all(i) = alpha_percent_all;
% end
a;
b;
c;
sigmaa;
sigmab;
sigmac;
xsubs;
ysubs;
zsubs;
sigmax;
sigmay;
sigmaz;
sigma_as = (std(result_as));
sigma_bs = (std(result_bs));
sigma_cs = (std(result_cs));
std_as = (std(result_as))/sqrt(j);
std_bs = (std(result_bs))/sqrt(j);
std_cs = (std(result_cs))/sqrt(j);
pd = mean(result_prob_dens);
pd_all = mean(result_prob_dens_all);
con_level = mean(result_cl);
con_level_all = mean(result_cl_all);
chisq_rem3_full = mean(result_chisq_rem3);
chisq_full_all = mean(result_chisq);
chisq_rem3_reduced = mean(result_chisq_rem3_red);
chisq_all_reduced = mean(result_chisq_red);
alpha_sig = mean(result_alpha);
alpha_sig_all = mean(result_alpha_all);
conlevel_percent = mean(result_cl_percent);
conlevel_percent_all = mean(result_cl_percent_all);
alpha_sig_percent = mean(result_alpha_percent);
alpha_sig_percent_all = mean(result_alpha_percent_all);
std_pd = (std(result_prob_dens))/sqrt(j);
std_pd_all = (std(result_prob_dens_all))/sqrt(j);
std_cl = (std(result_cl))/sqrt(j);
std_cl_all = (std(result_cl_all))/sqrt(j);
sigma_cl = std(result_cl);
sigma_cl_all = std(result_cl_all);
std_chisq_rem3_full = (std(result_chisq_rem3))/sqrt(j);
std_chisq_all = (std(result_chisq))/sqrt(j);
sigma_chisq_rem3_full = std(result_chisq_rem3);
sigma_chisq_all = std(result_chisq);
sigma_chisq_rem3_red = std(result_chisq_rem3_red);
sigma_chisq_red = std(result_chisq_red);
std_chisq_rem3_red = (std(result_chisq_rem3_red))/sqrt(j);
std_chisq_red = (std(result_chisq_red))/sqrt(j);
std_alpha_sig = (std(result_alpha))/sqrt(j);
std_alpha_sig_all = (std(result_alpha_all))/sqrt(j);
sigma_alpha_sig = std(result_alpha);
sigma_alpha_sig_all = std(result_alpha_all);
std_cl_percent = (std(result_cl_percent))/sqrt(j);
std_cl = (std(result_cl))/sqrt(j);
std_cl_percent_all = (std(result_cl_percent_all))/sqrt(j);
std_cl_all = (std(result_cl_all))/sqrt(j);
sigma_cl_percent = std(result_cl_percent);
sigma_cl = std(result_cl);
sigma_cl_percent_all = std(result_cl_percent_all);
sigma_cl_all = std(result_cl_all);
std_alpha_percent = (std(result_alpha_percent))/sqrt(j);
std_alpha_percent_all = (std(result_alpha_percent_all))/sqrt(j);
sigma_alpha_percent = std(result_alpha_percent);
sigma_alpha_percent_all = std(result_alpha_percent_all);

end
fprintf('                                         \n')
fprintf('                                         \n')
fprintf('                                         \n')
fprintf(' 3 observables, 2 free parameters and 1 parameter held fixed \n')
fprintf('a                                      %g \n',a)
fprintf('b                                      %g \n', b)
fprintf('c                                      %g \n', c)
fprintf('sigmarem3_1                            %g \n', sigmarem3_1)
fprintf('sigmarem3_2                            %g \n', sigmarem3_2)
fprintf('sigmarem3_3                            %g \n', sigmarem3_3)
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
fprintf('stddev_mean_as                         %g \n', std_as)
fprintf('stddev_mean_bs                         %g \n', std_bs)
fprintf('stddev_mean_cs                         %g \n', std_cs)
fprintf('probability                            %g \n', pd)
fprintf('standard dev. of mean of probability   %g \n', std_pd)
fprintf('CONFIDENCE LEVEL (percent)             %g \n', conlevel_percent)
fprintf('CONFIDENCE LEVEL                       %g \n', con_level)
fprintf('standard deviation of mean of c.l(per) %g \n', std_cl_percent)
fprintf('standard deviation of mean of c.l.     %g \n', std_cl)
fprintf('SIGMA_CL_PERCENT (percent)             %g \n', sigma_cl_percent)
fprintf('standard deviation of c.l.             %g \n', sigma_cl)
fprintf('significance (alpha)                   %g \n',alpha_sig_percent)
fprintf('standard deviation of mean of alpha    %g \n',std_alpha_percent)
fprintf('sigma_alpha_percent                    %g \n', sigma_alpha_percent)
fprintf('chi^2                                  %g \n', chisq_rem3_full)
fprintf('standard deviation of mean of chi^2    %g \n', std_chisq_rem3_full)
fprintf('sigma_chisq_rem3_full                  %g \n', sigma_chisq_rem3_full)
fprintf('chi^2_reduced                          %g \n', chisq_rem3_reduced)
fprintf('standard dev. of mean of chi^2_red     %g \n',std_chisq_rem3_red)
fprintf('sigma_chisq_rem3_red                   %g \n',sigma_chisq_rem3_red)

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
fprintf('stddev_mean_as                         %g \n', std_as)
fprintf('stddev_mean_bs                         %g \n', std_bs)
fprintf('stddev_mean_cs                         %g \n', std_cs)
fprintf('probability                            %g \n', pd_all)
fprintf('standard dev. of mean of probability   %g \n', std_pd_all)
fprintf('CONFIDENCE LEVEL (percent)             %g \n', conlevel_percent_all)
fprintf('CONFIDENCE LEVEL                       %g \n', con_level_all)
fprintf('standard deviation of mean of c.l(per) %g \n', std_cl_percent_all)
fprintf('standard deviation of mean of c.l.     %g \n', std_cl_all)
fprintf('SIGMA_CL_PERCENT                       %g \n', sigma_cl_percent_all)
fprintf('standard deviation of c.l.             %g \n', sigma_cl_all)
fprintf('significance (alpha)                   %g \n',alpha_sig_percent_all)
fprintf('standard deviation of mean of alpha    %g \n',std_alpha_percent_all)
fprintf('sigma_alpha_percent                    %g \n', sigma_alpha_percent_all)
fprintf('chi^2                                  %g \n', chisq_full_all)
fprintf('standard deviation of mean of chi^2    %g \n', std_chisq_all)
fprintf('sigma_chisq                            %g \n', sigma_chisq_all)
fprintf('chi^2_reduced                          %g \n', chisq_all_reduced)
fprintf('standard dev. of mean of chi^2_red     %g \n',std_chisq_red)
fprintf('sigma_chisq_red                        %g \n',sigma_chisq_red)
as = linspace(min(result_as),max(result_as),j);
bs = linspace(min(result_bs),max(result_bs),j);
cs = linspace(min(result_cs),max(result_cs),j);
chisquared2d= Chisq2(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
chisquared3d= Chisq3(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
cl_linspace = gammainc(chisquared2d/2,1);
cl_linspace_all = gammainc(chisquared3d/2,3/2);
chisquared2d_data= Chisq2(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
chisquared3d_data= Chisq3(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
cl_data = gammainc(chisquared2d_data/2,1);
cl_data_all = gammainc(chisquared3d_data/2,3/2);

% fprintf('a_data_lin     b_data_lin     c_data_lin     cl_lin(2d.o.f.)     cl_lin(3d.o.f)    cl_lin(2d.o.f.)-clin(3d.o.f.)\n')
% [as',bs',cs',cl_linspace',cl_linspace_all',(cl_linspace-cl_linspace_all)']
% fprintf('Mean of CL_linspace for 2 d.o.f. = %g \n',mean(cl_linspace))
% fprintf('Mean of CL_linspace for 3 d.o.f. =  %g \n',mean(cl_linspace_all))
% fprintf('Mean of CL_linspace(2 d.o.f.) -  Mean of CL_linspace(3 d.o.f.) =   %g \n',mean(cl_linspace)-mean(cl_linspace_all))
fprintf('                                         \n')
fprintf('                                         \n')
fprintf('                                         \n')
% [mean(cl_linspace),mean(cl_linspace_all)]
% fprintf('a_data     b_data     c_data     cl(2d.o.f.)     cl(3d.o.f)      cl_data(2d.o.f.)-cl_data(3d.o.f.)\n')
% [result_as',result_bs',result_cs',result_cl',result_cl_all',(result_cl - result_cl_all)']
fprintf('Mean of CL_data for A = %g \n',mean(result_cl))
fprintf('Mean of C.L_data for B =  %g \n',mean(result_cl_all))
fprintf('Mean of (CL_data(B) -  CL_data(A)) =  %g \n',mean(result_cl_all-result_cl))
fprintf('Standard deviation of (CL_data(B) - CL_data(A)) =  %g \n',std(result_cl_all - result_cl))
fprintf('Standard deviation of the mean of (CL_data(B) - CL_data(A)) =  %g \n',std(result_cl_all - result_cl)/sqrt(j))
[n,x]=hist(result_cl_all);
[n1,x1]=hist(result_cl);
[n2,x2]=hist((result_cl-result_cl_all));
% fprintf('Histogram frequency results for CL_data(2 d.o.f) = %g \n',[n,x])
% fprintf('Histogram frequency results for CL_data(3 d.o.f) = %g \n',[n1,x1])
% fprintf('Histogram frequency results for (CL_data(2 d.o.f)-CL_data(3 d.o.f)) = %g \n',[n2,x2])

% mean(cl_data)
% mean(cl_data_all)
% [mean(result_cl), mean(result_cl_all)]
% 
% 
% a = 1.1000;
% b = 1.5000;
% c = 5.8000;
% sigmaa = 0.001:0.001:0.1; % standard deviations of parameters chosen by me
% sigmab = 0.001:0.001:0.1;
% sigmac = 0.001:0.001:0.1;
% sigmax  = 0.0030;
% sigmay  = 0.3000;
% sigmaz  = 0.5000;
% as = 1.0975:0.00005:1.1025;

% as = linspace(min(result_as),max(result_as),j);
% bs = linspace(min(result_bs),max(result_bs),j);
% cs = linspace(min(result_cs),max(result_cs),j);
% as= reshape(result_as,length(result_as));
% bs= reshape(result_bs,length(result_bs));
% cs = reshape(result_cs,length(result_cs));
% % as =min(result_as):max(result_as);
% % bs =min(result_bs):max(result_bs);
% % cs =min(result_cs):max(result_cs);
% bs = 1.4975:0.00005:1.5025;
% cs = 5.7975:0.00005:5.8025;
% chisquared2d= Chisq2(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% chisquared2d=reshape(chisquare2d,length(result_as));
% chisquared3d= Chisq3(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% chisquared3d=reshape(chisquare3d,length(result_as));

% % [X1,X2] = meshgrid(as,bs);
% % chisquare3d = chisquared3d(X1(:),X2(:));
% % chisquare3d=reshape(chisquare3d,length(bs),length(as));
% figure(1)
% plot3(as,bs,chisquared2d)
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of two free parameters, a, b');
% grid on
% title(['\chi^2 vs. a_{data} and b_{data} for case of a =',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% print -deps plot2_as_bs_chi2_2dof;
% 
% figure(2)
% plot3(as,cs,chisquared2d)
% xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of two free parameters, a, c');
% grid on
% title(['\chi^2 vs. a_{data} and c_{data} for case of a =',num2str(a),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% print -deps plot3_as_cs_chi2_2dof;
% 
% 
% cl = gammainc(chisquared2d/2,1);
% figure(3)
% plot3(as,bs,cl)
% % xlim([0 max(as)]);
% % ylim([0 max(bs)]);
% % zlim([0 1.1]);
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Confidence level values for case of two free parameters, a, b')
% title(['Confidence level vs. a_{data} and b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% grid on
% print -deps plot2_as_bs_cl_2dof;
% 
figure(4)
plot(chisquared2d,cl)
xlim([-50.0 max(chisquared2d)]);
ylim([0 2]);
xlabel('\chi^{2} values for case of two free parameters, a, b'); ylabel('Confidence level values for case of two free parameters, a, b')
title(['Confidence level vs.\chi^2 for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ';  2 free parameters & 3 observables']);
grid on
print -deps plot_as_bs_cl_chi2_2dof;
% 
% figure(5)
% plot3(as,cs,cl)
% xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of two free parameters, a, c')
% title(['Confidence level vs. a_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% grid on
% print -deps plot3_as_cs_cl_2dof;
% 
% figure(6)
% plot3(bs,cs,chisquared2d)
% xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of two free parameters, b, c');
% grid on
% title(['\chi^2 vs. b_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% print -deps plot4_bs_cs_chi2_2dof;
% 
% figure(7)
% plot3(bs,cs,cl)
% xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of two free parameters, b, c')
% title(['Confidence level vs. b_{data} and c_{data} for case of a =',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% grid on
% print -deps plot4_bs_cs_cl_2dof;
% 
% chisquared3d= Chisq3(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% % % [X1,X2] = meshgrid(as,bs);
% % % chisquare3d = chisquared3d(X1(:),X2(:));
% % % chisquare3d=reshape(chisquare3d,length(bs),length(as));
% figure(8)
% plot3(as,bs,chisquared3d)
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of three free parameters, a, b, c');
% grid on
% title(['\chi^2 vs. a_{data} & b_{data} for case of a =',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% print -deps plot2_as_bs_chi2_3dof;
% 
% figure(9)
% plot3(as,cs,chisquared3d)
% xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of 3 free parameters, a, b, c');
% grid on
% title(['\chi^2 vs. a_{data} & c_{data} for case of a =',num2str(a),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% print -deps plot3_as_cs_chi2_3dof;
% 
% 
% cl_all = gammainc(chisquared3d/2,3/2);
% figure(10)
% plot3(as,bs,cl_all);
% % xlim([0 max(as)]);
% % ylim([0 max(bs)]);
% % zlim([0 1.5]);
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Confidence level values for case of three free parameters, a, b, c')
% title(['Confidence level vs. a_{data} and b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% grid on
% print -deps plot2_as_bs_cl_3dof;
% 
figure(11)
plot(chisquared3d,cl_all)
xlim([-100.0 max(chisquared3d)]);
ylim([0 2]);
xlabel('Chi^{2} values for case of 3 free parameters, a, b, c'); ylabel('Confidence level values for case of 3 free parameters, a, b, c')
title(['Confidence level vs.\chi^2 for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ';  3 free parameters & 3 observables']);
grid on
print -deps plot_cl_chi2_3dof;
% 
% figure(12)
% plot3(as,cs,cl_all)
% xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of 3 free parameters, a, b, c')
% title(['Confidence level vs. a_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% grid on
% print -deps plot3_as_cs_cl_3dof;
% 
% figure(13)
% plot3(bs,cs,chisquared3d)
% xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of 3 free parameters, a, b, c');
% grid on
% title(['\chi^2 vs. b_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% print -deps plot4_bs_cs_chi2_3dof;

% figure(14)
% plot3(bs,cs,cl_all)
% xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of 3 free parameters, a, b, c')
% title(['Confidence level vs. b_{data} & c_{data} for case of a =',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% grid on
% print -deps plot4_bs_cs_cl_3dof;

figure(14)
plot(chisquared2d,chisquared3d)
xlim([-1.0 max(chisquared2d)]);
ylim([0 max(chisquared3d)]);
xlabel('\chi^2(A) values'); ylabel('\chi^2(B) values')
title(['\chi^2 (A) vs.\chi^2 (B) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_chi2abcfree_chi2abfree;

figure(15)
plot(chisquared2d/2,chisquared3d/3)
xlim([-1.0 max(chisquared2d/2)]);
ylim([0 max(chisquared3d/3)]);
xlabel('Reduced \chi^2 (A) values'); ylabel('Reduced \chi^2 (B) values')
title(['Reduced \chi^2 (A) vs. Reduced \chi^2 (B) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_chi2redA_chi2redB;

figure(16)
plot(chisquared2d,chisquared3d-chisquared2d)
xlim([-1.0 max(chisquared2d)]);
ylim([-100.0 max(chisquared3d-chisquared2d)]);
xlabel('\chi^2(A)'); ylabel('\chi^2(B) - \chi^2(A)')
title(['\chi^2(B)-\chi^2(A) vs.\chi^2(A) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_chi2B-chi2A_vs_chi2A;

figure(17)
plot(chisquared2d/2,chisquared3d/3-chisquared2d/2)
xlim([-10.0 max(chisquared2d/2)]);
ylim([-10.0 max(chisquared3d/3-chisquared2d/2)]);
xlabel('\chi^2_reduced(A) '); ylabel('\chi^2_reduced(B) - \chi^2_reduced(A)')
title(['\chi^2_reduced(B)-\chi^2_reduced(A) vs.\chi^2(A) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_chi2Bred-chi2Ared_vs_chi2Ared;

% fprintf('chi-squared_A     chi-squared_B\n    chi-squared_B - chi-squared_A\n')
% [chisquared2d',chisquared3d',chisquared3d'-chisquared2d']
% fprintf('Mean of chi-squared for A = %g \n',mean(chisquared2d))
% fprintf('Mean of chi-squared for B =  %g \n',mean(chisquared3d))
% fprintf('Mean of (chisqr(B) -  chisqr(A)) =  %g \n',mean(chisquared3d-chisquared2d))
% fprintf('Standard deviation of ((chisqr(B)-chisqr(A)) =  %g \n',std(chisquared3d-chisquared2d))
% fprintf('Standard deviation of the mean of chisqr(B)-chisqr(A) =  %g \n',std(chisquared3d-chisquared2d)/sqrt(j))
% fprintf('                                         \n')
% fprintf('                                         \n')
% fprintf('chi-squared_reduced_A     chi-squared_reduced_B    chi-squared_reduced_B - chi-squared_reduced_A\n')
% [chisquared2d'/2,chisquared3d'/3,chisquared3d'/3 - chisquared2d'/2]
% fprintf('Mean of reduced chi-squared for A = %g \n',mean(chisquared2d/2))
% fprintf('Mean of reduced chi-squared for B =  %g \n',mean(chisquared3d/3))
% fprintf('Mean of (reduced chisqr(B) -  reduced chisqr(A)) =  %g \n',mean(chisquared3d/3-chisquared2d/2))
% fprintf('Standard deviation of (reduced chisqr(B) - reduced chisqr(A)) =  %g \n',std(chisquared3d/3-chisquared2d/2))
% fprintf('Standard deviation of the mean of (reduced chisqr(B) - reduced chisqr(A)) =  %g \n',std(chisquared3d/3-chisquared2d/2)/sqrt(j))

% figure(16)
% xspan = [0.0:0.001:1000];
% h=hist(chisqred2d,xspan);
% hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% bar(xspan,hn)
% hist(chisquared2d,xspan)
% title(['histogram of chi-squared(A) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
% ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% print -deps hist_chisqr_A
% 
% figure(17)
% xspan = [0.0:std(result_cl_all)/sqrt(j):1000];
% h=hist((result_cl),xspan);
% hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% bar(xspan,hn)
% hist(chisquared3d,xspan)
%  title(['histogram of chi-squared(B) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
% ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% print -deps hist_chisqr_B
% 
% figure(18)
% xspan = [0.0:std(result_cl - result_cl_all)/sqrt(j):1000];
% h=hist(chisquaerd2d-chisquared3d,xspan);
% hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% bar(xspan,hn)
% hist(chisquared2d-chisquared3d,xspan)
%  title(['histogram of chi-squared(A) - chi-squared(B) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
% ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% print -deps hist_chisqr_B



figure(18)
plot(cl_linspace,cl_linspace_all)
xlim([-1.1 1.1]);
ylim([0 1.1]);
xlabel('C.L.(A)'); ylabel('C.L.(B)')
title(['C.L. (A) vs. C.L.(B) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_CLabcfree_CLabfree;

figure(19)
plot(cl_linspace ,cl_linspace - cl_linspace_all)
xlim([-1.1 1.1]);
ylim([-20.0 1.1]);
xlabel('C.L.(A)'); ylabel('C.L.(A)-C.L.(B)')
title(['C.L.(A) - C.L.(B) vs. C.L.(A) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_CLA-ClBvsCLA;


% 
% figure(18)
% 
% % cl1 = gammainc((chisquared3d)/2,3/2);
% % cl2 = gammainc((chisquared2d)/2,2/2);
% plot3(as,bs,cl_all, ':*', as,bs,cl,'b-')
% title(['Confidence level vs. a_{data} and b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 & 3 free parameters & 3 observables']);
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('C.L.for case of two and three free parameters');
% grid on
% legend('3 free parameters (a,b,c)','2 free parameters (a,b)',0);
% print -deps plot_cl_2and3dof;
% 
% 
% 
% figure (19)
% plot3(as,bs,chisquared2d, ':*', as,bs,chisquared3d,'b-')
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of two and three free parameters');
% grid on
% legend('3 free parameters (a,b,c)','2 free parameters (a,b)',0);
% title(['\chi^2 (a,b,c, free parameters) and \chi^2 (a,b, free parameters for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% 
% print -deps plot_chi3dabcfree_chi2dabfree;
% 
% 
% chisquared2d= Chisq2(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
% chisquared3d= Chisq3(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
% 
% figure(20)
% plot3(result_as,result_bs,chisquared2d,'*')
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of 2 free parameters, a, b');
% grid on
% title(['\chi^2 vs. a_{data} and b_{data} for the case of a =',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% print -deps plot2_as_bs_chi2_2dof_data;
% 
% figure(21)
% plot3(result_as,result_cs,chisquared2d,'*')
% xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of 2 free parameters, a, c');
% grid on
% title(['\chi^2 vs. a_{data} & c_{data} for case of a =',num2str(a),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% print -deps plot3_as_cs_chi2_2dof_data;
% 
% 
% % cl = gammainc(chisquared2d/2,1);
% figure(22)
% plot3(result_as,result_bs,result_cl,'*')
% % xlim([0 max(as)]);
% % ylim([0 max(bs)]);
% % zlim([0 1.1]);
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Confidence level values for case of 2 free parameters, a, b')
% title(['Confidence level vs. a_{data} & b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% grid on
% print -deps plot2_as_bs_cl_2dof_data;

chisquared2d= Chisq2(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
chisquared3d= Chisq3(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);

figure(23)
plot(chisquared2d,result_cl,'*')
xlim([-50.0 max(chisquared2d)]);
ylim([0 2]);
xlabel('\chi^{2} values for case of 2 free parameters, a, b'); ylabel('Confidence level values for case of 2 free parameters, a, b')
title(['Graph of confidence level vs.\chi^2 for the case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ';  2 free parameters & 3 observables']);
grid on
print -deps plot_as_bs_cl_chi2_2dof_data;
% 
% figure(24)
% plot3(result_as,result_cs,result_cl,'*')
% xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of 2 free parameters, a, c')
% title(['Confidence level vs. a_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% grid on
% print -deps plot3_as_cs_cl_2dof_data;
% 
% figure(25)
% plot3(result_bs,result_cs,chisquared2d,'*')
% xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of 2 free parameters, b, c');
% grid on
% title(['\chi^2 vs. b_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% print -deps plot4_bs_cs_chi2_2dof_data;
% 
% figure(26)
% plot3(result_bs,result_cs,result_cl,'*')
% xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of two free parameters, b, c')
% title(['Confidence level vs. b_{data} & c_{data} for case of a =',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  2 free parameters & 3 observables']);
% grid on
% print -deps plot4_bs_cs_cl_2dof_data;
% 
% % chisquared3d= Chisq3(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% % % [X1,X2] = meshgrid(as,bs);
% % % chisquare3d = chisquared3d(X1(:),X2(:));
% % % chisquare3d=reshape(chisquare3d,length(bs),length(as));
% figure(27)
% plot3(result_as,result_bs,chisquared3d,'*')
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Chi^{2} values for case of 3 free parameters, a, b, c');
% grid on
% title(['\chi^2 vs. a_{data} & b_{data} for case of a =',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% print -deps plot2_as_bs_chi2_3dof_data;
% 
% figure(28)
% plot3(result_as,result_cs,chisquared3d,'*')
% xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of 3 free parameters, a, b, c');
% grid on
% title(['\chi^2 vs. a_{data} & c_{data} for case of a =',num2str(a),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% print -deps plot3_as_cs_chi2_3dof_data;
% 
% 
% % cl_all = gammainc(chisquared3d/2,3/2);
% figure(29)
% plot3(result_as,result_bs,result_cl_all,'*')
% % xlim([0 max(as)]);
% % ylim([0 max(bs)]);
% % zlim([0 1.5]);
% xlabel('a_{data}'); ylabel('b_{data}'); zlabel('Confidence level values for case of 3 free parameters, a, b, c')
% title(['Confidence level vs. a_{data} and b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% grid on
% print -deps plot2_as_bs_cl_3dof_data;
% 
% figure(30)
% plot(chisquared3d,result_cl_all,'*')
% xlim([-100.0 max(chisquared3d)]);
% ylim([0 2]);
% xlabel('\chi^{2} values for case of 3 free parameters, a, b, c'); ylabel('Confidence level values for case of 3 free parameters, a, b, c')
% title(['Confidence level vs.\chi^2 for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%      ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%      ';  3 free parameters & 3 observables']);
% grid on
% print -deps plot_cl_chi2_3dof_data;
% 
% figure(31)
% plot3(result_as,result_cs,result_cl_all,'*')
% xlabel('a_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of 3 free parameters, a, b, c')
% title(['Confidence level vs. a_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% grid on
% print -deps plot3_as_cs_cl_3dof_data;
% 
% figure(32)
% plot3(result_bs,result_cs,chisquared3d,'*')
% xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Chi^{2} values for case of 3 free parameters, a, b, c');
% grid on
% title(['\chi^2 vs. b_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% print -deps plot4_bs_cs_chi2_3dof_data;
% 
% figure(33)
% plot3(result_bs,result_cs,result_cl_all,'*')
% xlabel('b_{data}'); ylabel('c_{data}'); zlabel('Confidence level values for case of 3 free parameters, a, b, c')
% title(['Confidence level vs. b_{data} & c_{data} for case of a =',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
%     ';  3 free parameters & 3 observables']);
% grid on
% print -deps plot4_bs_cs_cl_3dof_data;
% 
chisquared2d= Chisq2(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
chisquared3d= Chisq3(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);

figure(34)
plot(chisquared2d,chisquared3d,'*')
xlim([-1.0 max(chisquared2d)]);
ylim([0 max(chisquared3d)]);
xlabel('\chi^2(A)(data)'); ylabel('\chi^2(B)(data)')
title(['\chi^2 (B)(data) vs.\chi^2 (A)(data) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_chi2abcfree_chi2abfree_data;

figure(35)
plot(result_cl,result_cl_all,'*')
xlim([-1.1 1.1]);
ylim([0 1.1]);
xlabel('C.L.(A)(data)'); ylabel('C.L.(B)(data)')
title(['C.L.(A)(data) vs. C.L.(B)(data) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_CLabcfree_CLabfree_data;

figure(36)
plot(result_cl, result_cl - result_cl_all,'*')
xlim([-1.1 1.1]);
ylim([-20.0 1.1]);
xlabel('C.L.(A)(data)'); ylabel('C.L.(A)(data)-C.L.(B)(data)')
title(['C.L. (A)(data)  vs. C.L.(B)(data) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_CL_A-CL_B_vs_CL_A_data;

figure(37)
xspan = [-1.0:std(result_cl-result_cl_all)/sqrt(j):1.1];
% h=hist((result_cl - result_cl_all),xspan);
% hn=h/((length(result_cl - result_cl_all)*(std(result_cl-result_cl_all)/sqrt(j))));
% bar(xspan,hn)
hist((result_cl - result_cl_all),xspan)
 title(['histogram of C.L.(A)(data) - C.L.(B)(data) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
print -deps hist_CLabfree-CLabcfree

figure(38)
xspan = [0.0:std(result_cl_all)/sqrt(j):1.1];
% h=hist((result_cl_all),xspan);
% hn=h/((length(result_cl_all)*(std(result_cl_all)/sqrt(j))));
% bar(xspan,hn)
hist(result_cl_all,xspan)
 title(['histogram of C.L.(B)(data) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
print -deps hist_CLabcfree

figure(39)
xspan = [0.0:std(result_cl)/sqrt(j):1.1];
% h=hist((result_cl),xspan);
% hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% bar(xspan,hn)
hist(result_cl,xspan)
 title(['histogram of C.L(A)(data) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
print -deps hist_CLabfree

figure(40)
p1 = wchi2cdf(chisquared2d_data,2); 
p2 = wchi2cdf(chisquared3d_data,3);
plot(p1,p2,'*')
grid on
xlabel('Confidence Level for theorist A (data)');
ylabel('Confidence Level for theorist B (data)');
title('Confidence Level for theorist B (data) vs. Confidence Level for theorist A (data)');
print -deps plot_chisq_A_B_data

% figure(41)
% plot(chisquared2d_data,chisquared3d_data,'*')
% xlim([0 max(chisquared2d_data)]);
% ylim([0 max(chisquared3d_data)]);
% xlabel('\chi^2 (A)(data)'); ylabel('\chi^2 (B)(data)')
% title(['\chi^2(B) vs.\chi^2(A) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
%     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% grid on
% print -deps plot_chi2B_data_chi2A_data;

figure(42)
plot(chisquared2d_data/2,chisquared3d_data/3,'*')
xlim([0 max(chisquared3d_data)]);
ylim([0 max(chisquared3d_data)]);
xlabel('Reduced \chi^2 (A)(data)'); ylabel('Reduced \chi^2 (B)(data)')
title(['Reduced \chi^2(B) vs. Reduced \chi^2(A) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_chi2B_data_red_chi2A_data_red;

figure(43)
plot(chisquared2d_data,chisquared3d_data-chisquared2d_data,'*')
xlim([-1.0 max(chisquared2d_data)]);
ylim([-100.0 max(chisquared3d_data-chisquared2d_data)]);
xlabel('\chi^2(A)(data)'); ylabel('\chi^2(B)(data) - \chi^2(A)(data)')
title(['\chi^2(B)(data)-\chi^2(A)(data) vs.\chi^2(A)(data) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_chi2B-chi2A_vs_chi2A_data;

figure(44)
plot(chisquared2d_data/2,chisquared3d_data/3-chisquared2d_data/2,'*')
xlim([-10.0 max(chisquared2d_data/2)]);
ylim([-10.0 max(chisquared3d_data/3-chisquared2d_data/2)]);
xlabel('\chi^2_reduced(A)(data) '); ylabel('\chi^2_reduced(B)(data) - \chi^2_reduced(A)(data)')
title(['\chi^2_reduced(B)(data)-\chi^2_reduced(A)(data) vs.\chi^2(A)(data) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_chi2Bred-chi2Ared_vs_chi2Ared_data;

figure(45)
sigmax = linspace(0.001,1.00,j);
plot(sigmax,result_cl,'*');
xlim([min(sigmax) max(sigmax)]);
xlabel('\sigma_x'); ylabel('Confidence Level (A)');
title(['Confidence Level (A) vs. \sigma_x for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
grid on
print -deps plot_CL_A-sigmax_data;


% fprintf('chi-squared_A_data     chi-squared_B_data    chi-squared_A_data - chi-squared_B_data\n')
% [chisquared2d_data',chisquared3d_data',chisquared3d_data' - chisquared2d_data']
fprintf('                                         \n')
fprintf('                                         \n')
fprintf('Mean of chi-squared for A (data) = %g \n',mean(chisquared2d_data))
fprintf('Mean of chi-squared for B (data)=  %g \n',mean(chisquared3d_data))
fprintf('Mean of (chisqr(B) -  chisqr(A)) =  %g \n',mean(chisquared3d_data - chisquared2d_data))
fprintf('Standard deviation of ((chisqr(B)-chisqr(A)) (data) =  %g \n',std(chisquared3d_data-chisquared2d_data))
fprintf('Standard deviation of the mean of chisqr(B)-chisqr(A) (data) =  %g \n',std(chisquared3d_data-chisquared2d_data)/sqrt(j))
fprintf('                                         \n')

% fprintf('chi-squared_A_data_red     chi-squared_B_data_red    chi-squared_A_data_red - chi-squared_B_data_red\n')
% [chisquared2d_data'/2,chisquared3d_data'/3,chisquared3d_data'/3-chisquared2d_data'/2]
fprintf('                                         \n')
fprintf('                                         \n')
fprintf('Mean of reduced chi-squared for A (data) = %g \n',mean(chisquared2d_data/2))
fprintf('Mean of reduced chi-squared for B (data) =  %g \n',mean(chisquared3d_data/3))
fprintf('Mean of (reduced chisqr(B) - reduced  chisqr(A)) (data) =  %g \n',mean(chisquared3d_data/3 - chisquared2d_data/2))
fprintf('Standard deviation of (reduced chisqr(B)- reduced chisqr(A)) (data) =  %g \n',std(chisquared3d_data/3-chisquared2d_data/2))
fprintf('Standard deviation of the mean of reduced chisqr(B) - reduced chisqr(A) (data) =  %g \n',std(chisquared3d_data/3-chisquared2d_data/2)/sqrt(j))



% figure(45)
% xspan = [0.0:std(result_cl)/sqrt(j):1000];
% h=hist(chisqred2d_data,xspan);
% hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% bar(xspan,hn)
% hist(chisquared2d_data,xspan)
 % title(['histogram of chi-squared(A) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
% ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% print -deps hist_chisqr_A_data

% figure(46)
% xspan = [0.0:std(result_cl_all)/sqrt(j):1000];
% h=hist((result_cl),xspan);
% hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% bar(xspan,hn)
% hist(chisquared3d_data,xspan)
%  title(['histogram of chi-squared(B) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
% ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% print -deps hist_chisqr_B_data

% figure(47)
% xspan = [0.0:std(result_cl - result_cl_all)/sqrt(j):1000];
% h=hist(chisquaerd2d-chisquared3d,xspan);
% hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% bar(xspan,hn)
% hist(chisquared2d_data - chisquared3d_data,xspan)
%  title(['histogram of chi-squared(A) - chi-squared(B) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
% ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% print -deps hist_chisqr_B_data
