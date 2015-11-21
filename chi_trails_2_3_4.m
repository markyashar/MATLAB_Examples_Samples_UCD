clear all
syms a b c d as bs cs ds as_draw bs_draw cs_draw ds_draw sigma_as sigma_bs sigma_cs sigma_ds sigmaa... 
    sigmab sigmac sigmad sigmax sigmay sigmaz sigmav real;
x = a*b*c*d; y = a^2 + b^2 + c^2 + d^2; z = a + b + c +d; v = a*d^2 + c*b^2;
Fdata = [1/(sigmax)^2 0 0 0; 0 1/(sigmay)^2 0 0; 0 0 1/(sigmaz)^2 0; 0 0 0 1/(sigmav)^2];
M = jacobian([x;y;z;v],[a,b,c,d]);
F = M'*Fdata*M;
% a = 2.9; b = 5.9; c = 0.8; sigmax = 0.005; sigmay = 0.009; sigmaz = 0.003;
% a = 2.9 
% b = 5.9 
% c = 0.8
% a = 21.96;
a = 1.16;
b = 1.63;
c = 1.18;
d = 1.55;
% sigmaa = 0.2;
% sigmaa = 0.00001;
sigmaa = 0.030;
sigmab = 0.002;
sigmac = 0.003;
sigmad = 0.004;
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
sigmav = 0.007;
% sigmax = sqrt((sigmaa^2)*(b*c*d)^2 +(sigmab^2)*(a*c*d)^2 + (sigmac^2)*(a*b*d)^2 + (sigmad^2)*(a*b*c)^2);
% sigmay = sqrt((sigmaa^2)*(2*a)^2 + (sigmab^2)*(2*b)^2 + (sigmac^2)*(2*c)^2 + (sigmad^2)*(2*d)^2 );
% sigmaz = sqrt(sigmaa^2 + sigmab^2 + sigmac^2 + sigmad^2);
% sigmav = sqrt((sigmaa^2)*(d^2)^2 + (sigmab^2)*(2*b*c)^2 + (sigmac^2)*(b^2)^2 + (sigmad^2)*(2*d*a)^2);
Fdata_subs = subs(Fdata);
Fsubs = subs(F);
eigsFsubs = eig(Fsubs);
eig1 = eigsFsubs(1);
eig2 = eigsFsubs(2);
eig3 = eigsFsubs(3);
eig4 = eigsFsubs(4);
sigma1 = 1/sqrt(eigsFsubs(1));
sigma2 = 1/sqrt(eigsFsubs(2));
sigma3 = 1/sqrt(eigsFsubs(3));
sigma4 = 1/sqrt(eigsFsubs(4));
Msubs = subs(M);
fom = sqrt(det(Fsubs));
xsubs = subs(x);
ysubs = subs(y);
zsubs = subs(z);
vsubs = subs(v);
precisionx = (sigmax/xsubs)*100;
precisiony = (sigmay/ysubs)*100;
precisionz = (sigmaz/zsubs)*100;
precisionv = (sigmav/vsubs)*100;
Fsubsrem4 = Fsubs(1:3,1:3);
fom_rem4 = sqrt(det(Fsubsrem4));
eigsFsubsrem4 = eig(Fsubsrem4);
eigrem4_1 = eigsFsubsrem4(1);
eigrem4_2 = eigsFsubsrem4(2);
% eigrem3_3 = eigsFsubsrem3(3)
eigrem4_3 = eigsFsubsrem4(3);
% eigrem4_4 = eigsFsubs(4);

sigmarem4_1 = 1/sqrt(eigsFsubsrem4(1));
sigmarem4_2 = 1/sqrt(eigsFsubsrem4(2));
sigmarem4_3 = 1/sqrt(eigsFsubsrem4(3));

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
j=10000;
for i = 1:j;
% as_draw=normrnd(a,sigmaa,[1]);
% bs_draw=normrnd(b,sigmab,[1]);
% cs_draw=normrnd(c,sigmac,[1]);
% as=as_draw;
as=normrnd(a,sigmaa,[1]);
bs=normrnd(b,sigmab,[1]);
cs=normrnd(c,sigmac,[1]);
ds=normrnd(d,sigmad,[1]);
result_as(i) = as;
% bs=bs_draw;
result_bs(i) = bs;
% cs=cs_draw;
result_cs(i) = cs;
result_ds(i) = ds;
% sigma_as =(std(as_draw))/1000
% sigma_bs = (std(bs_draw))/1000
% sigma_cs = (std(cs_draw))/1000
pt = [a;b];
ptall =[a;b;c];
ptall4 = [a;b;c;d];
% ps = [xs;ys;zs];
% ps = [as;bs];
ps = [as;bs];
% psall =  [as;bs;cs];
psall = [as;bs;cs];
psall4 = [as;bs;cs;ds];
% chisq_rem3 = (pt-ps)'*Fsubsrem3*(pt-ps);
chisq = (ptall-psall)'*Fsubsrem4*(ptall-psall);
chisq_all4 = (ptall4-psall4)'*Fsubs*(ptall4-psall4);
% result_chisq_rem3(i)=chisq_rem3;
result_chisq(i)=chisq;
result_chisq_all4(i)=chisq_all4;

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
% prob_dens = ((2*pi)^-(2/2))*sqrt(det(Fsubsrem3))*exp(-(chisq_rem4)/2);
% prob_dens_all = ((2*pi)^-(3/2))*sqrt(det(Fsubs))*exp(-chisq/2);
% result_prob_dens(i) = prob_dens;
% result_prob_dens_all(i) = prob_dens_all;
% n=2;
cl = gammainc(chisq/2,3/2);
cl_all4 = gammainc((chisq_all4)/2,4/2);
result_cl(i) =cl;
result_cl_all4(i) = cl_all4;
chisq_red(i) = chisq/3;
chisq_all4_red(i) = chisq_all4/4;
% result_chisq_red(i) = chisq_red;
result_chisq_red = chisq_red;
% result_chisq_all4_red(i) = chisq_all4_red;
result_chisq_all4_red = chisq_all4_red;

% alpha(i) = 1-cl;
% alpha_all(i) = 1 - cl_all;
% result_alpha(i) = 1-cl;
% result_alpha_all(i) = 1-cl_all;
%cl_percent(i) = 100*cl;
%cl_percent_all(i) = 100*cl_all;
% result_cl_percent(i) = 100*cl;
% result_cl_percent_all(i) = 100*cl_all;
% result_alpha_percent(i) = (1-cl)*100;
% result_alpha_percent_all(i) = (1-cl_all)*100;
% result_alpha_percent(i) = alpha_percent;
% result_alpha_percent_all(i) = alpha_percent_all;
% end
a;
b;
c;
d;
sigmaa;
sigmab;
sigmac;
sigmad;
xsubs;
ysubs;
zsubs;
vsubs;
sigmax;
sigmay;
sigmaz;
sigmav;
sigma_as = (std(result_as));
sigma_bs = (std(result_bs));
sigma_cs = (std(result_cs));
sigma_ds = (std(result_ds));
std_as = (std(result_as))/sqrt(j);
std_bs = (std(result_bs))/sqrt(j);
std_cs = (std(result_cs))/sqrt(j);
std_ds = (std(result_ds))/sqrt(j);
% pd = mean(result_prob_dens);
% pd_all = mean(result_prob_dens_all);
con_level = mean(result_cl);
con_level_all4 = mean(result_cl_all4);
chisq_full = mean(result_chisq);
chisq_full_all4 = mean(result_chisq_all4);
chisq_reduced = mean(result_chisq_red);
chisq_all4_reduced = mean(result_chisq_all4_red);
% alpha_sig = mean(result_alpha);
% alpha_sig_all = mean(result_alpha_all);
% conlevel_percent = mean(result_cl_percent);
% conlevel_percent_all = mean(result_cl_percent_all);
% alpha_sig_percent = mean(result_alpha_percent);
% alpha_sig_percent_all = mean(result_alpha_percent_all);
% std_pd = (std(result_prob_dens))/sqrt(j);
% std_pd_all = (std(result_prob_dens_all))/sqrt(j);
std_cl = (std(result_cl))/sqrt(j);
std_cl_all4 = (std(result_cl_all4))/sqrt(j);
sigma_cl = std(result_cl);
sigma_cl_all4 = std(result_cl_all4);
std_chisq_full = (std(result_chisq))/sqrt(j);
std_chisq_all4 = (std(result_chisq_all4))/sqrt(j);
sigma_chisq_full = std(result_chisq);
sigma_chisq_all4 = std(result_chisq_all4);
sigma_chisq_red = std(result_chisq_red);
sigma_chisq_all4_red = std(result_chisq_all4_red);
std_chisq_red = (std(result_chisq_red))/sqrt(j);
std_chisq_all4_red = (std(result_chisq_all4_red))/sqrt(j);
% std_alpha_sig = (std(result_alpha))/sqrt(j);
% std_alpha_sig_all = (std(result_alpha_all))/sqrt(j);
% sigma_alpha_sig = std(result_alpha);
% sigma_alpha_sig_all = std(result_alpha_all);
% std_cl_percent = (std(result_cl_percent))/sqrt(j);
std_cl = (std(result_cl))/sqrt(j);
% std_cl_percent_all = (std(result_cl_percent_all))/sqrt(j);
std_cl_all4 = (std(result_cl_all4))/sqrt(j);
% sigma_cl_percent = std(result_cl_percent);
sigma_cl = std(result_cl);
% sigma_cl_percent_all = std(result_cl_percent_all);
sigma_cl_all4 = std(result_cl_all4);
% std_alpha_percent = (std(result_alpha_percent))/sqrt(j);
% std_alpha_percent_all = (std(result_alpha_percent_all))/sqrt(j);
% sigma_alpha_percent = std(result_alpha_percent);
% sigma_alpha_percent_all = std(result_alpha_percent_all);

end
fprintf('                                         \n')
fprintf('                                         \n')
fprintf('                                         \n')
fprintf(' 4 observables, 3 free parameters and 1 parameter held fixed \n')
fprintf('a                                      %g \n',a)
fprintf('b                                      %g \n', b)
fprintf('c                                      %g \n', c)
fprintf('d                                      %g \n', d)
fprintf('sigmarem4_1                            %g \n', sigmarem4_1)
fprintf('sigmarem4_2                            %g \n', sigmarem4_2)
fprintf('sigmarem4_3                            %g \n', sigmarem4_3)
% fprintf('sigmarem4_4                            %g \n', sigmarem4_4)
fprintf('sigma_a                                %g \n', sigmaa)
fprintf('sigma_b                                %g \n', sigmab)
fprintf('sigma_c                                %g \n', sigmac)
fprintf('sigma_d                                %g \n', sigmad)
fprintf('x                                      %g \n', xsubs)
fprintf('y                                      %g \n', ysubs)
fprintf('z                                      %g \n', zsubs)
fprintf('v                                      %g \n', vsubs)
fprintf('sigma_x                                %g \n', sigmax)
fprintf('sigma_y                                %g \n', sigmay)
fprintf('sigma_z                                %g \n', sigmaz)
fprintf('sigma_v                                %g \n', sigmav)
fprintf('sigma_as                               %g \n', sigma_as)
fprintf('sigma_bs                               %g \n', sigma_bs)
fprintf('sigma_cs                               %g \n', sigma_cs)
fprintf('sigma_ds                               %g \n', sigma_ds)
fprintf('stddev_mean_as                         %g \n', std_as)
fprintf('stddev_mean_bs                         %g \n', std_bs)
fprintf('stddev_mean_cs                         %g \n', std_cs)
fprintf('stddev_mean_ds                         %g \n', std_ds)
% fprintf('probability                            %g \n', pd)
% fprintf('standard dev. of mean of probability   %g \n', std_pd)
% fprintf('CONFIDENCE LEVEL (percent)             %g \n', conlevel_percent)
fprintf('CONFIDENCE LEVEL                       %g \n', con_level)
% fprintf('standard deviation of mean of c.l(per) %g \n', std_cl_percent)
fprintf('standard deviation of mean of c.l.     %g \n', std_cl)
% fprintf('SIGMA_CL_PERCENT (percent)             %g \n', sigma_cl_percent)
fprintf('standard deviation of c.l.             %g \n', sigma_cl)
% fprintf('significance (alpha)                   %g \n',alpha_sig_percent)
% fprintf('standard deviation of mean of alpha    %g \n',std_alpha_percent)
% fprintf('sigma_alpha_percent                    %g \n', sigma_alpha_percent)
fprintf('chi^2                                  %g \n', chisq_full)
fprintf('standard deviation of mean of chi^2    %g \n', std_chisq_full)
fprintf('sigma_chisq_rem4_full                  %g \n', sigma_chisq_full)
fprintf('chi^2_rem4_reduced                          %g \n', chisq_reduced)
fprintf('standard dev. of mean of chi^2_red     %g \n',std_chisq_red)
fprintf('sigma_chisq_rem4_red                   %g \n',sigma_chisq_red)

fprintf('                                         \n')
fprintf('                                         \n')
fprintf('                                         \n')
fprintf(' 4 observables, 4 free parameters \n')
fprintf('a                                      %g \n',a)
fprintf('b                                      %g \n', b)
fprintf('c                                      %g \n', c)
fprintf('d                                      %g \n', d)
fprintf('sigma_1                                %g \n', sigma1)
fprintf('sigma_2                                %g \n', sigma2)
fprintf('sigma_3                                %g \n', sigma3)
fprintf('sigma_4                                %g \n', sigma4)
fprintf('sigma_a                                %g \n', sigmaa)
fprintf('sigma_b                                %g \n', sigmab)
fprintf('sigma_c                                %g \n', sigmac)
fprintf('sigma_d                                %g \n', sigmad)
fprintf('x                                      %g \n', xsubs)
fprintf('y                                      %g \n', ysubs)
fprintf('z                                      %g \n', zsubs)
fprintf('v                                      %g \n', vsubs)
fprintf('sigma_x                                %g \n', sigmax)
fprintf('sigma_y                                %g \n', sigmay)
fprintf('sigma_z                                %g \n', sigmaz)
fprintf('sigma_v                                %g \n', sigmav)
fprintf('sigma_as                               %g \n', sigma_as)
fprintf('sigma_bs                               %g \n', sigma_bs)
fprintf('sigma_cs                               %g \n', sigma_cs)
fprintf('sigma_ds                               %g \n', sigma_ds)
fprintf('stddev_mean_as                         %g \n', std_as)
fprintf('stddev_mean_bs                         %g \n', std_bs)
fprintf('stddev_mean_cs                         %g \n', std_cs)
fprintf('stddev_mean_ds                         %g \n', std_ds)
% fprintf('probability                            %g \n', pd_all)
% fprintf('standard dev. of mean of probability   %g \n', std_pd_all)
% fprintf('CONFIDENCE LEVEL (percent)             %g \n', conlevel_percent_all)
fprintf('CONFIDENCE LEVEL                       %g \n', con_level_all4)
% fprintf('standard deviation of mean of c.l(per) %g \n', std_cl_percent_all)
fprintf('standard deviation of mean of c.l.     %g \n', std_cl_all4)
% fprintf('SIGMA_CL_PERCENT                       %g \n', sigma_cl_percent_all)
fprintf('standard deviation of c.l.             %g \n', sigma_cl_all4)
% fprintf('significance (alpha)                   %g \n',alpha_sig_percent_all)
% fprintf('standard deviation of mean of alpha    %g \n',std_alpha_percent_all)
% fprintf('sigma_alpha_percent                    %g \n', sigma_alpha_percent_all)
fprintf('chi^2                                  %g \n', chisq_full_all4)
fprintf('standard deviation of mean of chi^2    %g \n', std_chisq_all4)
fprintf('sigma_chisq                            %g \n', sigma_chisq_all4)
fprintf('chi^2_reduced                          %g \n', chisq_all4_reduced)
fprintf('standard dev. of mean of chi^2_red     %g \n',std_chisq_all4_red)
fprintf('sigma_chisq_red                        %g \n',sigma_chisq_all4_red)
% as = linspace(min(result_as),max(result_as),j);
% bs = linspace(min(result_bs),max(result_bs),j);
% cs = linspace(min(result_cs),max(result_cs),j);
% ds = linspace(min(result_ds),max(result_ds),j);
% chisquared4d= Chisq2(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% chisquared3d= Chisq3(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% cl_linspace = gammainc(chisquared2d/2,1);
% cl_linspace_all = gammainc(chisquared3d/2,3/2);
% chisquared2d_data= Chisq2(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
% chisquared3d_data= Chisq3(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
% cl_data = gammainc(chisquared2d_data/2,1);
% cl_data_all = gammainc(chisquared3d_data/2,3/2);

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
fprintf('Mean of C.L_data for B =  %g \n',mean(result_cl_all4))
fprintf('Mean of (CL_data(B) -  CL_data(A)) =  %g \n',mean(result_cl_all4-result_cl))
fprintf('Standard deviation of (CL_data(B) - CL_data(A)) =  %g \n',std(result_cl_all4 - result_cl))
fprintf('Standard deviation of the mean of (CL_data(B) - CL_data(A)) =  %g \n',std(result_cl_all4 - result_cl)/sqrt(j))
fprintf('                                         \n')
fprintf('                                         \n')
fprintf('Mean of chi-squared for A (data) = %g \n',mean(result_chisq))
fprintf('Mean of chi-squared for B (data)=  %g \n',mean(result_chisq_all4))
fprintf('Mean of (chisqr(B) -  chisqr(A)) =  %g \n',mean(result_chisq_all4 - result_chisq))
fprintf('Standard deviation of ((chisqr(B)-chisqr(A)) (data) =  %g \n',std(result_chisq_all4 - result_chisq))
fprintf('Standard deviation of the mean of chisqr(B)-chisqr(A) (data) =  %g \n',std(result_chisq_all4 - result_chisq)/sqrt(j))
fprintf('                                         \n')

% fprintf('chi-squared_A_data_red     chi-squared_B_data_red    chi-squared_A_data_red - chi-squared_B_data_red\n')
% [chisquared2d_data'/2,chisquared3d_data'/3,chisquared3d_data'/3-chisquared2d_data'/2]
fprintf('                                         \n')
fprintf('                                         \n')
fprintf('Mean of reduced chi-squared for A (data) = %g \n',mean(result_chisq/3))
fprintf('Mean of reduced chi-squared for B (data) =  %g \n',mean(result_chisq_all4/4))
fprintf('Mean of (reduced chisqr(B) - reduced  chisqr(A)) (data) =  %g \n',mean(result_chisq_all4/4 - result_chisq/3))
fprintf('Standard deviation of (reduced chisqr(B)- reduced chisqr(A)) (data) =  %g \n',std(result_chisq_all4/4 - result_chisq/3))
fprintf('Standard deviation of the mean of reduced chisqr(B) - reduced chisqr(A) (data) =  %g \n',std(result_chisq_all4/4 - result_chisq/3)/sqrt(j))

% chisquared2d_data= Chisq2(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
% chisquared3d_data= Chisq3(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
cl_data = gammainc(result_chisq/2,3/2);
cl_data_all4 = gammainc(result_chisq_all4/2,4/2);


figure(1)
plot3(result_as,result_bs,cl_data_all4,'*')
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('C.L.(B)')
title(['C.L.(B) vs. a_{data} & b_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]);
grid on
print -deps plot3_CLall4_B_as_bs_data;


figure(2)
plot3(result_as,result_bs,cl_data,'*')
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('C.L.(A)')
title(['C.L.(A) vs. a_{data} & b_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
 ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
 ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]);   
grid on
print -deps plot3_CLall_A_as_bs_data;

figure(3)
plot3(result_as,result_cs,cl_data_all4,'*')
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('C.L.(B)')
title(['C.L.(A) vs. a_{data} & c_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
 ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
 ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]); 
grid on
print -deps plot3_CLall4_B_as_cs_data;


figure(4)
plot3(result_as,result_cs,cl_data,'*')
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('C.L.(A)')
title(['C.L.(A) vs. a_{data} & c_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
 ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
 ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]); 
grid on
print -deps plot3_CLall_A_as_cs_data;


figure(5)
plot3(result_bs,result_cs,cl_data_all4,'*')
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('C.L.(B)')
title(['C.L.(B) vs. b_{data} & c_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
 ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
 ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]); 
grid on
print -deps plot3_CLall4_B_bs_cs_data;


figure(6)
plot3(result_bs,result_cs,cl_data,'*')
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('C.L.(A)')
title(['C.L.(A) vs. b_{data} & c_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
 ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
 ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]); 
grid on
print -deps plot3_CLall_A_bs_cs_data;


figure(7)
plot3(result_as,result_bs,result_chisq_all4,'*')
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('\chi^2 (B)')
title(['\chi^2 (B) vs. a_{data} & b_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]);
grid on
print -deps plot3_chisq_all4_B_as_bs_data;


figure(8)
plot3(result_as,result_bs,result_chisq,'*')
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('\chi^2 (A)')
title(['\chi^2 (A) vs. a_{data} & b_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]);
grid on
print -deps plot3_chisq_A_as_bs_data;


figure(9)
plot3(result_as,result_cs, result_chisq_all4,'*')
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('\chi^2 (B)')
title(['\chi^2 (B) vs. a_{data} & c_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]);
grid on
print -deps plot3_chisq_all4_B_as_cs_data;


figure(10)
plot3(result_as,result_cs,result_chisq,'*')
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('\chi^2 (A)')
title(['\chi^2 (A) vs. a_{data} & c_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]);
grid on
print -deps plot3_chisq_all_A_as_cs_data;



figure(11)
plot3(result_bs,result_cs,result_chisq_all4,'*')
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('\chi^2 (B)')
title(['\chi^2 (B) vs. b_{data} & c_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]);
grid on
print -deps plot3_chisq_all4_B_bs_cs_data;



figure(12)
plot3(result_bs,result_cs,result_chisq,'*')
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('\chi^2 (A)')
title(['\chi^2 (A) vs. b_{data} & c_{data}, a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),',d =',num2str(d),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),', \sigma_{v} =',num2str(sigmav),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac),', \sigma_{d} =',num2str(sigmad)]);
grid on
print -deps plot3_chisq_all_A_bs_cs_data;




