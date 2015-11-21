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
j=5000;
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
% fprintf('cl(A) cl(B) cl(A)-cl(B) chi^2(A) chi^2(B) chi^2(B)-chi^2(A) chi^2(A)_red chi^2(B)_red chi^2(B)_red-chi^2(A)_red\n')
% [result_cl',result_cl_all',(result_cl_all-result_cl)',chisquared2d_data',chisquared3d_data',chisquared3d_data'-chisquared2d_data',chisquared2d_data'/2,chisquared3d_data'/3,chisquared3d_data'/3-chisquared2d_data'/2]
% as = linspace(min(result_as),max(result_as),j);
% bs = linspace(min(result_bs),max(result_bs),j);
% cs = linspace(min(result_cs),max(result_cs),j);
% chisquared2d= Chisq2(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% chisquared3d= Chisq3(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
% cl_linspace = gammainc(chisquared2d/2,1);
% cl_linspace_all = gammainc(chisquared3d/2,3/2);
chisquared2d_data= Chisq2(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
chisquared3d_data= Chisq3(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
cl_data = gammainc(chisquared2d_data/2,1);
cl_data_all = gammainc(chisquared3d_data/2,3/2);

% fprintf('cl(A) cl(B) cl(A)-cl(B) chi^2(A) chi^2(B) chi^2(B)-chi^2(A) chi^2(A)_red chi^2(B)_red chi^2(B)_red-chi^2(A)_red\n')
% [cl_data',cl_data_all',cl_data_all'-cl_data',chisquared2d_data',chisquared3d_data',chisquared3d_data'-chisquared2d_data',chisquared2d_data'/2,chisquared3d_data'/3,chisquared3d_data'/3-chisquared2d_data'/2]

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
fprintf('Mean of CL_data for A = %g \n',mean(cl_data))
fprintf('Mean of C.L_data for B =  %g \n',mean(cl_data_all))
fprintf('Mean of (CL_data(B) -  CL_data(A)) =  %g \n',mean(result_cl_all-result_cl))
fprintf('Standard deviation of (CL_data(B) - CL_data(A)) =  %g \n',std(result_cl_all - result_cl))
fprintf('Standard deviation of the mean of (CL_data(B) - CL_data(A)) =  %g \n',std(result_cl_all - result_cl)/sqrt(j))
fprintf('                                         \n')
fprintf('                                         \n')
fprintf('Mean of chi-squared for A (data) = %g \n',mean(chisquared2d_data))
fprintf('Mean of chi-squared for B (data)=  %g \n',mean(chisquared3d_data))
fprintf('Mean of (chisqr(B) -  chisqr(A)) =  %g \n',mean(abs(chisquared3d_data - chisquared2d_data)))
fprintf('Standard deviation of ((chisqr(B)-chisqr(A)) (data) =  %g \n',std(abs(chisquared3d_data-chisquared2d_data)))
fprintf('Standard deviation of the mean of chisqr(B)-chisqr(A) (data) =  %g \n',std(abs(chisquared3d_data-chisquared2d_data)/sqrt(j)))
fprintf('                                         \n')

% fprintf('chi-squared_A_data_red     chi-squared_B_data_red    chi-squared_A_data_red - chi-squared_B_data_red\n')
% [chisquared2d_data'/2,chisquared3d_data'/3,chisquared3d_data'/3-chisquared2d_data'/2]
fprintf('                                         \n')
fprintf('                                         \n')
fprintf('Mean of reduced chi-squared for A (data) = %g \n',mean(chisquared2d_data/2))
fprintf('Mean of reduced chi-squared for B (data) =  %g \n',mean(chisquared3d_data/3))
fprintf('Mean of (reduced chisqr(B) - reduced  chisqr(A)) (data) =  %g \n',mean(abs(chisquared3d_data/3 - chisquared2d_data/2)))
fprintf('Standard deviation of (reduced chisqr(B)- reduced chisqr(A)) (data) =  %g \n',std(abs(chisquared3d_data/3-chisquared2d_data/2)))
fprintf('Standard deviation of the mean of reduced chisqr(B) - reduced chisqr(A) (data) =  %g \n',std(abs(chisquared3d_data/3-chisquared2d_data/2)/sqrt(j)))

figure(1)
plot3(result_as,result_bs,chisquared3d_data,'O');
% zdir = [0 0 1];
% center = [10 10 0];
% rotate(h,zdir,45,center)
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('\chi^2 (B)')
title(['\chi^2 (B) vs. a_{data} & b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqB_as_bs_data;


figure(2)
plot3(result_as,result_bs,chisquared2d_data,'O');
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('\chi^2 (A)')
title(['\chi^2 (A) vs. a_{data} & b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqA_as_bs_data;


figure(3)
plot3(result_as,result_bs,abs(chisquared3d_data-chisquared2d_data),'O');
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('\chi^2 (A)')
title(['\chi^2(B)-\chi^2(A) vs. a_{data} & b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqB-chisqA_as_bs_data;


figure(4)
plot3(result_as,result_cs,chisquared3d_data,'O');
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('\chi^2 (B)')
title(['\chi^2 (B) vs. a_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqB_as_cs_data;


figure(5)
plot3(result_as,result_cs,chisquared2d_data,'O');
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('\chi^2 (A)')
title(['\chi^2 (A) vs. a_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqA_as_cs_data;

figure(6)
plot3(result_as,result_cs,abs(chisquared3d_data-chisquared2d_data),'O');
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('\chi^2 (B)')
title(['\chi^2 (B) vs. a_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqB-chisqA_as_cs_data;



figure(7)
plot3(result_bs,result_cs,chisquared3d_data,'O');
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('\chi^2 (B)')
title(['\chi^2 (B) vs. b_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqB_bs_cs_data;


figure(8)
plot3(result_bs,result_cs,chisquared2d_data,'O');
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('\chi^2 (A)')
title(['\chi^2 (A) vs. b_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqA_bs_cs_data;


figure(9)
plot3(result_bs,result_cs,abs(chisquared3d_data-chisquared2d_data),'O');
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('\chi^2(B)-\chi^2(A)')
title(['\chi^2 (B) vs. b_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqB-chisqA_bs_cs_data;


figure(10)
plot3(result_as,result_bs,cl_data_all,'O');
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('C.L.(B)')
title(['C.L.(B) vs. a_{data} & b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_CL_B_as_bs_data;


figure(11)
plot3(result_as,result_bs,cl_data,'O');
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('C.L.(A)')
title(['C.L.(A) vs. a_{data} & b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_CL_A_as_bs_data;

figure(12)
plot3(result_as,result_bs,cl_data_all-cl_data,'O');
xlabel('a_{data}'); ylabel('b_{data}'); zlabel('C.L.(B)-C.L.(A)')
title(['C.L.(B)-C.L.(A) vs. a_{data} & b_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_CLB-CLA_as_bs_data;

figure(13)
plot3(result_as,result_cs,cl_data_all,'O');
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('C.L.(O)')
title(['C.L.(B) vs. a_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_CL_B_as_cs_data;


figure(14)
plot3(result_as,result_cs,cl_data,'O');
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('C.L.(A)')
title(['C.L.(A) vs. a_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_CL_A_as_cs_data;

figure(15)
plot3(result_as,result_cs,cl_data_all-cl_data,'O');
xlabel('a_{data}'); ylabel('c_{data}'); zlabel('C.L.(B)-C.L.(A)')
title(['C.L.(B)-C.L.(A) vs. a_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_CLB-CLA_as_cs_data;


figure(16)
plot3(result_bs,result_cs,cl_data_all,'O');
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('C.L.(B)')
title(['C.L.(B) vs. b_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_CL_B_bs_cs_data;


figure(17)
plot3(result_bs,result_cs,cl_data,'O');
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('C.L.(A)')
title(['C.L.(A) vs. b_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_CL_A_bs_cs_data;

figure(18)
plot3(result_bs,result_cs,cl_data_all-cl_data,'O');
xlabel('b_{data}'); ylabel('c_{data}'); zlabel('C.L.(B)-C.L.(A)')
title(['C.L.(B)-C.L.(A) vs. b_{data} & c_{data} for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_CLB-CLA_bs_cs_data;

figure(19)
plot3(chisquared2d_data,chisquared3d_data,abs(chisquared3d_data-chisquared2d_data),'O');
xlabel('\chi^2 (A)'); ylabel('\chi^2 (B)'); zlabel('\chi^2(B)-\chi^2(A)')
title(['\chi^2 (B)-\chi^2(A) vs.  & \chi^2 (B) , \chi^2(A) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqB-chisqA_chisqA_chisqB_data;

figure(20)
plot3(chisquared2d_data/2,chisquared3d_data/3,abs(chisquared3d_data/3-chisquared2d_data/2),'O');
xlabel('reduced \chi^2 (A)'); ylabel('reduced \chi^2 (B)'); zlabel('reduced \chi^2(B)-reduced \chi^2(A)')
title(['\chi^2 (B)-\chi^2(A) vs.  & \chi^2 (B) , \chi^2(A) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_chisqBred-chisqAred_chisqAred_chisqBred_data;

figure(21)
plot3(cl_data,cl_data_all,cl_data_all-cl_data,'O')
xlabel('C.L.(A)'); ylabel('C.L.(B)'); zlabel('C.L.(B)-C.L.(A)')
title(['C.L.(B)-C.L.(A) vs. C.L.(A) & C.L.(B) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot3_CLB-CLA_CLA_CLB_data;


figure(22)
axis([-10.0 1.0e5 0.0 2.0]);
plot(chisquared2d_data,cl_data,'O')
axis([-10.0 1.0e5 0.0 2.0])
xlim([-10.0 max(chisquared2d_data)])
ylim([0 2])
xlabel('\chi^2 (A)'); ylabel('C.L.(A)')
title(['C.L.(A) vs.\chi^2 (A), a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_as_bs_CLA_chisqA_data;

figure(23)
plot(chisquared3d_data,cl_data_all,'O');
% xlim([-100.0 max(chisquared3d_data)]);
% ylim([0 2]);
xlabel('\chi^2 (B)'); ylabel('C.L. (B)')
title(['C.L.(B) vs.\chi^2(B), a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLB_chisqB_data;

figure(24)
xlim([-2.0 max(cl_data)]);
ylim([0 max(cl_data_all)])
plot(cl_data,cl_data_all,'O')
% xlim([-2.0 max(cl_data)])
% ylim([0 max(cl_data_all)])
xlabel('C.L.(A)'); ylabel('C.L.(B)')
title(['C.L.(B)) vs. C.L.(A) , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLB_CLA_data;

figure(25)
plot(cl_data,cl_data_all-cl_data,'O');
xlim([-2.0 max(cl_data)]);
ylim([0 max(cl_data_all-cl_data)]);
xlabel('C.L.(A)'); ylabel('C.L.(B)- C.L.(A)')
title(['C.L.(B)-C.L.(A) vs C.L.(A) , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLB-CLA_CLA_data;

figure(26)
plot(cl_data_all,cl_data_all-cl_data,'O');
xlim([-2.0 max(cl_data_all)])
ylim([0 max(cl_data_all-cl_data)])
xlabel('C.L.(B)'); ylabel('C.L.(B)-C.L(A)')
title(['C.L.(B)-C.L.(A) vs. C.L.(B) , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLB-CLA_CLB_data;


figure(27)
plot(chisquared2d_data,chisquared3d_data,'O');
xlim([-1.0 max(chisquared2d_data)]);
ylim([0 max(chisquared3d_data)]);
xlabel('\chi^2(A)'); ylabel('\chi^2(B)')
title(['\chi^2 (A) vs.\chi^2 (B) , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chisqB_chisqA_data;

figure(28)
plot(chisquared2d_data,chisquared3d_data-chisquared2d_data,'O');
xlim([-1.0 max(chisquared2d_data)]);
ylim([0 max(chisquared3d_data-chisquared2d_data)]);
xlabel('\chi^2(A)'); ylabel('\chi^2(B)-\chi^2(A)')
title(['\chi^2(B)-\chi^2(A) vs.\chi^2(A) , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chisqB-chisqA_chisqA_data;

figure(29)
plot(chisquared3d_data,chisquared3d_data-chisquared2d_data,'O');
xlim([-1.0 max(chisquared3d_data)]);
ylim([0 max(chisquared3d_data-chisquared2d_data)]);
xlabel('\chi^2(B)'); ylabel('\chi^2(B)-\chi^2(A)')
title(['\chi^2(B)-\chi^2(A) vs.\chi^2(B) , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chisqB-chisqA_chisqB_data;


figure(30)
plot(chisquared2d_data/2,chisquared3d_data/3,'O');
xlim([-1.0 max(chisquared2d_data/2)]);
ylim([0 max(chisquared3d_data/3)]);
xlabel('Reduced \chi^2(A)'); ylabel('\chi^2(B)')
title(['Reduced \chi^2 (A) vs. Reduced \chi^2 (B), = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2redA_chi2redB_data;

figure(31)
plot(chisquared2d_data/2,chisquared3d_data/3-chisquared2d_data/2,'O');
xlim([-1.0 max(chisquared2d_data/2)]);
ylim([-100.0 max(chisquared3d_data/3-chisquared2d_data/2)]);
xlabel('reduced \chi^2(A)'); ylabel('reduced \chi^2(B) - reduced \chi^2(A)')
title(['reduced \chi^2(B)- reduced \chi^2(A) vs.reduced \chi^2(A), a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chisqBred-chisqAred_vs_chisqAred_data;

figure(32)
plot(chisquared3d_data/2,chisquared3d_data/3-chisquared2d_data/2,'O');
xlim([-10.0 max(chisquared3d_data/2)]);
ylim([-10.0 max(chisquared3d_data/3-chisquared2d_data/2)]);
xlabel('\chi^2_{reduced}(B) '); ylabel('\chi^2_{reduced}(B) - \chi^2_{reduced}(A)')
title(['\chi^2_{reduced(B)}-\chi^2_{reduced}(A) vs.\chi^2_{reduced}(B), a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2Bred-chi2Ared_vs_chi2Bred_data;

figure(33)
sigmax=normrnd(0.004,0.003,[1 j])
plot(sigmax,cl_data,'O');
xlim([min(sigmax) max(sigmax)]);
xlabel('\sigma_x'); ylabel('Confidence Level (A)');
title(['C.L.(A) vs. \sigma_x , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CL_A-sigmax_data;


figure(34)
sigmay = normrnd(0.005,0.004,[1 j]);
plot(sigmay,cl_data,'O');
xlim([min(sigmay) max(sigmay)]);
xlabel('\sigma_y'); ylabel('Confidence Level (A)');
title(['C.L.(A) vs. \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CL_A-sigmay_data;

figure(35)
sigmaz = normrnd(0.006,0.005,[1 j]);
plot(sigmaz,cl_data,'O');
xlim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_z'); ylabel('Confidence Level (A)');
title(['C.L.(A) vs. \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CL_A-sigmaz_data;


figure(36)
plot(sigmax,cl_data_all,'O');
xlim([min(sigmax) max(sigmax)]);
xlabel('\sigma_x'); ylabel('Confidence Level (B)');
title(['C.L.(B) vs. \sigma_x , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CL_B-sigmax_data;


figure(37)
plot(sigmay,cl_data_all,'O');
xlim([min(sigmay) max(sigmay)]);
xlabel('\sigma_y'); ylabel('Confidence Level (B)');
title(['C.L.(B) vs. \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CL_B-sigmay_data;

figure(38)
plot(sigmaz,cl_data_all,'O');
xlim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_z'); ylabel('Confidence Level (B)');
title(['C.L.(B) vs. \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CL_B-sigmaz_data;


figure(39)
plot(sigmax,cl_data-cl_data_all,'O');
xlim([min(sigmax) max(sigmax)]);
xlabel('\sigma_x'); ylabel('C.L.(A) - C.L. (B)');
title(['C.L.(A) - C.L.(B) vs. \sigma_x , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLA-CLB-sigmax_data;


figure(40)
plot(sigmay,cl_data-cl_data_all,'O');
xlim([min(sigmay) max(sigmay)]);
xlabel('\sigma_y'); ylabel('C.L.(A) - C.L.(B)');
title(['C.L.(A)-C.L.(B) vs. \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLA-CLB-sigmay_data;

figure(41)
plot(sigmaz,cl_data-cl_data_all,'O');
xlim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_z'); ylabel('C.L.(A)-C.L.(B)');
title(['C.L.(A)-C.L.(B) vs. \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLA-CLB-sigmaz_data;

figure(42)
plot(sigmax,chisquared2d_data,'O');
xlim([min(sigmax) max(sigmax)]);
xlabel('\sigma_x'); ylabel('\chi^2 (A)');
title(['\chi^2 (A) vs. \sigma_x , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2_A-sigmax_data;

figure(43)
plot(sigmax,cl_data_all,'O');
xlim([min(sigmax) max(sigmax)]);
xlabel('\sigma_x'); ylabel('Confidence Level (B)');
title(['C.L.(B) vs. \sigma_x , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CL_B-sigmax_data;


figure(44)
plot(sigmay,cl_data_all,'O');
xlim([min(sigmay) max(sigmay)]);
xlabel('\sigma_y'); ylabel('Confidence Level (B)');
title(['C.L.(B) vs. \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CL_B-sigmay_data;

figure(45)
plot(sigmaz,cl_data_all,'O');
xlim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_z'); ylabel('Confidence Level (B)');
title(['C.L.(B) vs. \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CL_B-sigmaz_data;


figure(46)
plot(sigmax,cl_data-cl_data_all,'O');
xlim([min(sigmax) max(sigmax)]);
xlabel('\sigma_x'); ylabel('C.L.(A) - C.L. (B)');
title(['C.L.(A) - C.L.(B) vs. \sigma_x , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLA-CLB-sigmax_data;


figure(47)
plot(sigmay,cl_data-cl_data_all,'O');
xlim([min(sigmay) max(sigmay)]);
xlabel('\sigma_y'); ylabel('C.L.(A) - C.L.(B)');
title(['C.L.(A)-C.L.(B) vs. \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLA-CLB-sigmay_data;

figure(48)
plot(sigmaz,cl_data-cl_data_all,'O');
xlim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_z'); ylabel('C.L.(A)-C.L.(B)');
title(['C.L.(A)-C.L.(B) vs. \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLA-CLB-sigmaz_data;

figure(49)
plot(sigmax,chisquared2d_data,'O');
xlim([min(sigmax) max(sigmax)]);
xlabel('\sigma_x'); ylabel('\chi^2 (A)');
title(['\chi^2 (A) vs. \sigma_x , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2_A-sigmax_data;



figure(50)
plot(sigmay,chisquared2d_data,'O');
xlim([min(sigmay) max(sigmay)]);
xlabel('\sigma_y'); ylabel('\chi^2 (A)');
title(['\chi^2 (A) vs. \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2_A-sigmay_data;

figure(51)
plot(sigmaz,chisquared2d_data,'O');
xlim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_z'); ylabel('\chi^2 (A)');
title(['\chi^2 (A) vs. \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2_A-sigmaz_data;


figure(52)
plot(sigmax,chisquared3d_data,'O');
xlim([min(sigmax) max(sigmax)]);
xlabel('\sigma_x'); ylabel('\chi^2 (B)');
title(['\chi^2 (B) vs. \sigma_x , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2_B-sigmax_data;


figure(53)
plot(sigmay,chisquared3d_data,'O');
xlim([min(sigmay) max(sigmay)]);
xlabel('\sigma_y'); ylabel('\chi^2 (B)');
title(['\chi^2 (B) vs. \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2B-sigmay_data;

figure(54)
plot(sigmaz,chisquared3d_data,'O');
xlim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_z'); ylabel('\chi^2 (B)');
title(['\chi^2 (B) vs. \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2B-sigmaz_data;


figure(55)
plot(sigmax,chisquared2d_data-chisquared3d_data,'O');
xlim([min(sigmax) max(sigmax)]);
xlabel('\sigma_x'); ylabel('\chi^2(A) - \chi^2(B)');
title(['\chi^2(A) - \chi^2(B) vs. \sigma_x , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2A-chi2B-sigmax_data;


figure(56)
plot(sigmay,chisquared2d_data-chisquared3d_data,'O');
xlim([min(sigmay) max(sigmay)]);
xlabel('\sigma_y'); ylabel('\chi^2(A) - \chi^2(B)');
title(['\chi^2(A) - \chi^2(B) vs. \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2A-chi2B_sigmay_data;

figure(57)
plot(sigmaz,chisquared2d_data-chisquared3d_data,'O');
xlim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_z'); ylabel('\chi^2 (A) - \chi^2(B)');
title(['\chi^2 (A) - \chi^2(B) vs. \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2A-chi2B_sigmaz_data;

figure(58)
plot3(sigmax,sigmay,chisquared2d_data,'O');
xlim([min(sigmax) max(sigmax)]);
ylim([min(sigmay) max(sigmay)]);
xlabel('\sigma_x'); ylabel('\sigma_y'); zlabel('\chi^2 (A)');
title(['\chi^2 (A) vs. \sigma_x & \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2A_sigmax_sigmay_data;


figure(59)
plot3(sigmax,sigmay,chisquared3d_data,'O');
xlim([min(sigmax) max(sigmax)]);
ylim([min(sigmay) max(sigmay)]);
xlabel('\sigma_x'); ylabel('\sigma_y'); zlabel('\chi^2 (B)');
title(['\chi^2 (B) vs. \sigma_x , \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2B_sigmax_sigmay_data;

figure(60)
plot3(sigmax,sigmay,chisquared2d_data-chisquared3d_data,'O');
xlim([min(sigmax) max(sigmax)]);
ylim([min(sigmay) max(sigmay)]);
xlabel('\sigma_x'); ylabel('\sigma_y'); zlabel('\chi^2(A)-\chi^2(B)');
title(['\chi^2(A)-\chi^2(B) vs. \sigma_x & \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2A-chi2B-sigmax_sigmay_data;


figure(60)
plot3(sigmax,sigmaz,chisquared2d_data,'O');
xlim([min(sigmax) max(sigmax)]);
ylim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_x'); ylabel('\sigma_z'); zlabel('\chi^2 (A)');
title(['\chi^2 (A) vs. \sigma_x & \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2A_sigmax_sigmaz_data;


figure(61)
plot3(sigmax,sigmaz,chisquared3d_data,'O');
xlim([min(sigmax) max(sigmax)]);
ylim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_x'); ylabel('\sigma_z'); zlabel('\chi^2 (B)');
title(['\chi^2 (B) vs. \sigma_x , \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmax),', \sigma_{z} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2B_sigmax_sigmaz_data;

figure(62)
plot3(sigmay,sigmaz,chisquared2d_data-chisquared3d_data,'O');
xlim([min(sigmay) max(sigmay)]);
ylim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_y'); ylabel('\sigma_z'); zlabel('\chi^2(A)-\chi^2(B)');
title(['\chi^2(A)-\chi^2(B) vs. \sigma_y & \sigma_z , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_chi2A-chi2B-sigmay_sigmaz_data;


figure(63)
plot3(sigmax,sigmay,cl_data-cl_data_all,'O');
xlim([min(sigmay) max(sigmay)]);
ylim([min(sigmay) max(sigmay)]);
xlabel('\sigma_x'); ylabel('\sigma_y'); zlabel('C.L.(A)-C.L.(B)');
title(['C.L.(A)-C.L.(B) vs. \sigma_x & \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLA-CLB-sigmay_sigmaz_data;

figure(64)
plot3(sigmax,sigmaz,cl_data-cl_data_all,'O');
xlim([min(sigmax) max(sigmax)]);
ylim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_x'); ylabel('\sigma_y'); zlabel('C.L.(A)-C.L.(B)');
title(['C.L.(A)-C.L.(B) vs. \sigma_x & \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLA-CLB-sigmax_sigmaz_data;


figure(65)
plot3(sigmax,sigmaz,cl_data-cl_data_all,'O');
xlim([min(sigmax) max(sigmax)]);
ylim([min(sigmaz) max(sigmaz)]);
xlabel('\sigma_x'); ylabel('\sigma_y'); zlabel('C.L.(A)-C.L.(B)');
title(['C.L.(A)-C.L.(B) vs. \sigma_x & \sigma_y , a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
    ', \sigma_{x} =',num2str(sigmay),', \sigma_{y} =',num2str(sigmaz),...
    ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
grid on
print -deps plot_CLA-CLB-sigmax_sigmaz_data;

cl_data = cl_dataI;
c_data_all =  cl_data_allI;
chisquared2d_data = chisquared2d_dataI;
chisquared3d_data = chisquared3d_dataI;
result_as = result_asI;
result_bs = result_bsI;
result_cs = result_csI;



save cl_dataI;
save c_data_allI;
save chisquared2d_dataI;
save chisquared3d_dataI;
save result_asI;
save result_bsI;
save result_csI;


figure(67);
xspan = [0.0:0.01:100];
hist(chisquared3d_data-chisquared2d_data,xspan);
% % hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% % bar(xspan,hn)
title(['histogram of \chi^2 (B) - \chi^2 (A), a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
print -deps hist_chisqB-chisqA_data
 
figure(68);
xspan = [0.0:0.001:100];
hist(chisquared3d_data/3-chisquared2d_data/2,xspan);
ylim([0 100])
% % hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% % bar(xspan,hn)
title(['histogram of reduced \chi^2 (B) - reduced \chi^2 (A), a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
print -deps hist_chisqBred-chisqAred_data
 
 
figure(69);
xspan = [0.0:0.001:100];
hist(cl_data_all-cl_data,xspan);
% % hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% % bar(xspan,hn)
title(['histogram of C.L.(B) - C.L.(A), a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
     ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz),...
     ', \sigma_{a} =',num2str(sigmaa),', \sigma_{b} =',num2str(sigmab),', \sigma_{c} =',num2str(sigmac)]);
print -deps hist_CLA-CLB_data



% figure(28)
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





% figure(26)
% xspan = [0.0:std(cl_data)/sqrt(j):1000];
% h=hist(chisqred2d_data,xspan);
% hn=h/((length(cl_data)*(std(cl_data)/sqrt(j))));
% bar(xspan,hn)
% hist(chisquared2d_data,xspan)
% title(['histogram of \chi^2 (A), a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
% ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% print -deps hist_chisqr_A_data
% 
% figure(22)
% xspan = [0.0:std(result_cl_all)/sqrt(j):1000];
% h=hist((result_cl),xspan);
% hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% bar(xspan,hn)
% hist(chisquared3d_data,xspan)
% title(['histogram of chi-squared(B) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
% ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% print -deps hist_chisqr_B_data
% 
% figure(23)
% xspan = [0.0:std(result_cl - result_cl_all)/sqrt(j):1000];
% h=hist(chisquaerd2d-chisquared3d,xspan);
% hn=h/((length(result_cl)*(std(result_cl)/sqrt(j))));
% bar(xspan,hn)
% hist(chisquared2d_data - chisquared3d_data,xspan)
% title(['histogram of chi-squared(A) - chi-squared(B) for case of a = ',num2str(a),' ,b =',num2str(b),',c =',num2str(c),...
% ', \sigma_{x} =',num2str(sigmax),', \sigma_{y} =',num2str(sigmay),', \sigma_{z} =',num2str(sigmaz)]);
% print -deps hist_chisqr_B_data




