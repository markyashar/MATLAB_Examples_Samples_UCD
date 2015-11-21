clear all
j=100;
% for i = 1:j;
sigmax=linspace(0.005,.5,j);
sigmay=linspace(0.006,.6,j);
sigmaz=linspace(0.007,.7,j);
sigmaa=linspace(0.005,0.5,j);
sigmab=linspace(0.006,0.6,j);
sigmac=linspace(0.007,0.7,j);
a=linspace(1,5,j);
b=linspace(1,5,j);
c=linspace(1,5,j);
a=normrnd(a,sigmaa,[1 j]);
b=normrnd(b,sigmab,[1 j]);
c=normrnd(c,sigmac,[1 j]);
as=normrnd(a,sigmaa,[1 j]);
bs=normrnd(b,sigmab,[1 j]);
cs=normrnd(c,sigmac,[1 j]);
% result_as(i) = as;
% result_bs(i) = bs;
% result_cs(i) = cs;
% result_cl(i) =cl;
% result_cl_all(i) = cl_all;
% chisq_rem3_red = chisq_rem3/2;
% chisq_red(i) = chisq/n;
% result_chisq_rem3_red(i) = chisq_rem3_red;
% result_chisq_red = chisq_red;
% result_alpha(i) = 1-cl;
% result_alpha_all(i) = 1-cl_all;
% result_cl_percent(i) = 100*cl;
% result_cl_percent_all(i) = 100*cl_all;
% result_alpha_percent(i) = (1-cl)*100;
% result_alpha_percent_all(i) = (1-cl_all)*100;

chisquared2d= Chisq2(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
chisquared3d= Chisq3(a,b,c,as,bs,cs,sigmax,sigmay,sigmaz);
cl = gammainc(chisquared2d/2,1);
cl_all = gammainc(chisquared3d/2,3/2);
[sigmax',sigmay',sigmaz',sigmaa',sigmab',sigmac',a',b',c',as',bs',cs',cl',cl_all']
% chisquared2d= Chisq2(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);
% chisquared3d= Chisq3(a,b,c,result_as,result_bs,result_cs,sigmax,sigmay,sigmaz);


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
% bs = 1.4975:0.00005:1.5025;
% cs = 5.7975:0.00005:5.8025;