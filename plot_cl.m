
figure(1)

x = linspace(0.0,20.5,200);
p0 = wchi2cdf(x,1); p1 = wchi2cdf(x,2); p2 = wchi2cdf(x,3); p3 = wchi2cdf(x,4); p4 = wchi2cdf(x,5); p5 = wchi2cdf(x,6)
plot(x,p0,'*', x,p1,'b-',x,p2,'b--', x, p3,'b:+',x, p4,':ks')
grid on
legend('1 d.o.f.','2 d.o.f.','3 d.o.f.','4 d.o.f.','5 d.o.f.',0);
xlabel('\chi^2');
ylabel('Confidence Level');
title('Confidence Level vs. \chi^2 for different numbers of degrees of freedom');

figure(2)

plot(p0,p1,'*', p1,p2,'b-',p2,p3,'b--', p3, p4,'b:+',p4, p5,':ks')
grid on
legend('1 d.o.f. vs 2 d.o.f.','2 d.o.f. vs 3 d.o.f.','3 d.o.f. vs. 4 d.o.f.','4 d.o.f. vs 5 d.o.f.','5 d.o.f. vs. 6 d.o.f',0);
xlabel('Confidence Level, n degrees of freedom');
ylabel('Confidence Level, n+1 degrees of freedom');
title('Confidence Level for n+1 d.o.f. vs. Confidence Level for n d.o.f., where n = 1,2,3,4,5,6');
print -deps plot_cl_n_n+1

figure(3)

plot(p1,p5,'*',p1,p4,'b-',p1,p3,'b--',p1, p2,'b:+')
grid on
legend('2 d.o.f. vs 6 d.o.f.','2 d.o.f. vs 5 d.o.f.','2 d.o.f. vs. 4 d.o.f.','2 d.o.f. vs. 3 d.o.f.',0);
xlabel('Confidence Level, 2 degrees of freedom');
ylabel('Confidence Level, 3,4,5,6 degrees of freedom');
title('Confidence Level for 3,4,5,6 d.o.f. vs. Confidence Level for 2 d.o.f.');
print -deps plot_cl_2dof_3-6dof

figure(4)
y = linspace(0,10.5,200);
prob0 = wchi2pdf(y,1); prob1 = wchi2pdf(y,2); prob2 = wchi2pdf(y,3); prob3 = wchi2pdf(y,4); prob4 = wchi2pdf(y,9);
plot(y,prob0,'*', y,prob1,'b-',y,prob2,'b--', y, prob3,'b:+',y, prob4,':ks')
grid on
legend('1 d.o.f.','2 d.o.f.','3 d.o.f.','4 d.o.f.','9 d.o.f.',0);
xlabel('chi^2');
ylabel('Chi squared probability density function');
title('Chi squared prob. density function vs. Chi^2 for different numbers of degrees of freedom');
print -deps plot_density_chi_9dof
fprintf('C.L.(1 d.o.f)     C.L. (2 d.o.f.)     3 d.o.f     4 d.o.f.    5 d.o.f.     6 d.o.f.  %g \n')
          [p0'    , p1'    ,    p2'   ,    p3'   ,   p4'   ,  p5']
  fprintf('                                                                         \n')       
 fprintf('1 d.o.f-2 d.o.f.     2d.o.f-3d.o.f     3d.o.f.-4d.o.f.    4d.o.f.-5d.o.f.     5d.o.f.-6 d.o.f.  %g \n')
         [(p0-p1)'    , (p1-p2)'    ,    (p2-p3)'   ,    (p3-p4)'   ,   (p4-p5)']        
         
% figure(5)
% x = linspace(0.0,15.5,200);
% p1 = wchi2cdf(x,2); p2 = wchi2cdf(x,3); 
% plot3(p1,p2,x)
% grid on

figure(5)
y = linspace(0,20.5,200);
p1new = wchi2cdf(y,2); p2new = wchi2cdf(y+1,3); p3new = wchi2cdf(y+2,3); p4new = wchi2cdf(y+3,3); p5new = wchi2cdf(y+4,3); p6new = wchi2cdf(y+5,3);
plot(y,p1-p2,'b-',y,p1new-p2new,'*',y,p1new-p3new,'b--',y,p1new-p4new,'b:+',y,p1new-p5new,':ks')
grid on
legend('C.L. for 2 d.o.f - C.L. for 3 d.o.f for the same \chi^2','C.L. for 2 d.o.f and \chi^2 - C.L. for 3 d.o.f and \chi^2 + 1',...
    'C.L.(2 d.o.f , \Chi^2) - C.L.(3 d.o.f , \Chi^2 + 2)', 'C.L.(2 d.o.f , \chi^2) - C.L.(3 d.o.f , \chi^2 + 3)',...
    'C.L.(2 d.o.f , \chi^2) - C.L.(3 d.o.f , \chi^2 + 4)', 'C.L.(2 d.o.f , \chi^2) - C.L.(3 d.o.f , \chi^2 + 5)',0);
xlabel('\chi^2');
ylabel('C.L. for 2 d.o.f - C.L. for 3 d.o.f ');
title('Difference in C.L. vs. \chi^2 for 2 and 3 d.o.f.');
print -deps chi_CL_difference_2and3dof

figure(6)
plot(y,p1new,'*', y,p2new,'b-', y,p3new,'b--',y,p4new,'b:+',y,p5new,':ks')
grid on
legend('C.L.(2 d.o.f. and \chi^2)','C.L.(3 d.o.f. and \chi^2 + 1)','C.L.(3 d.o.f. and \chi^2 + 2)','C.L.(3 d.o.f. and \chi^2 + 3)',...
    'C.L.(3 d.o.f. and chi^2 + 4)',0);
plot(p1new,p2new,'b-')
grid on
legend('1 d.o.f. vs 2 d.o.f.','2 d.o.f. vs 3 d.o.f.','3 d.o.f. vs. 4 d.o.f.','4 d.o.f. vs 5 d.o.f.','5 d.o.f. vs. 6 d.o.f',0);
xlabel('\chi^2');
ylabel('Confidence Level with \chi^2 + n  and 3 d.o.f.');
title('Confidence Levels for different \chi^2 with 2 and 3 d.o.f.');
print -deps plot_cl_chi2_chi2+n

figure(7)
y = linspace(0,20.5,200);
p1new = wchi2cdf(y,2); p2new = wchi2cdf(y,3); p3new = wchi2cdf(y*2,3); p4new = wchi2cdf(y*3,3); p5new = wchi2cdf(y*4,3); p6new = wchi2cdf(y*5,3);
plot(y,p1-p2,'b-',y,p1new-p2new,'*',y,p1new-p3new,'b--',y,p1new-p4new,'b:+',y,p1new-p5new,':ks')
grid on
legend('C.L. for 2 d.o.f - C.L. for 3 d.o.f for the same \chi^2','C.L. for 2 d.o.f and \chi^2 - C.L. for 3 d.o.f and \chi^2 * 1',...
    'C.L.(2 d.o.f , \chi^2) - C.L.(3 d.o.f , \chi^2 * 2)', 'C.L.(2 d.o.f , \chi^2) - C.L.(3 d.o.f , \chi^2 * 3)',...
    'C.L.(2 d.o.f , \chi^2) - C.L.(3 d.o.f , \chi^2 * 4)', 'C.L.(2 d.o.f , \chi^2) - C.L.(3 d.o.f , \chi^2 * 5)',0);
xlabel('\chi^2');
ylabel('C.L. for 2 d.o.f - C.L. for 3 d.o.f ');
title('Difference in C.L. vs. \chi^2 for 2 and 3 d.o.f.');
print -deps chi_times_CL_difference_2and3dof


figure(8)
y = linspace(0,20.5,200);
plot(y,p1-p2,'b-',y,p1new-p2new,'*',y,p1new-p3new,'b--',y,p1new-p4new,'b:+',y,p1new-p5new,':ks')
grid on
legend('C.L. for 2 d.o.f - C.L. for 3 d.o.f for the same \chi^2','C.L. for 2 d.o.f and \chi^2 - C.L. for 3 d.o.f and \chi^2 * 1',...
    'C.L.(2 d.o.f , \Chi^2) - C.L.(3 d.o.f , \Chi^2 * 2)', 'C.L.(2 d.o.f , \chi^2) - C.L.(3 d.o.f , \chi^2 * 3)',...
    'C.L.(2 d.o.f , \chi^2) - C.L.(3 d.o.f , \chi^2 * 4)', 'C.L.(2 d.o.f , \chi^2) - C.L.(3 d.o.f , \chi^2 * 5)',0);
xlabel('\chi^2');
ylabel('C.L. for 2 d.o.f - C.L. for 3 d.o.f ');
title('Difference in C.L. vs. \chi^2 for 2 and 3 d.o.f.');
print -deps chi_times_CL_difference_2and3dof

figure(9)
plot(y,p1new,'*', y,p2new,'b-', y,p3new,'b--',y,p4new,'b:+',y,p5new,':ks')
grid on
legend('C.L.(2 d.o.f. and \chi^2)','C.L.(3 d.o.f. and \chi^2 * 1)','C.L.(3 d.o.f. and \chi^2 * 2)','C.L.(3 d.o.f. and \chi^2 * 3)',...
    'C.L.(3 d.o.f. and chi^2 * 4)',0);
plot(p1new,p2new,'b-')
grid on
legend('1 d.o.f. vs 2 d.o.f.','2 d.o.f. vs 3 d.o.f.','3 d.o.f. vs. 4 d.o.f.','4 d.o.f. vs 5 d.o.f.','5 d.o.f. vs. 6 d.o.f',0);
xlabel('\chi^2');
ylabel('Confidence Level with \chi^2 * n  and 3 d.o.f.');
title('Confidence Levels for different \chi^2 with 2 and 3 d.o.f.');
print -deps plot_cl_chi2_chi2_times_n
fprintf('chi-squared   C.L.(chi2 and 2 d.o.f.)     C.L.(chi2 + 1 and 3 d.o.f.)   C.L.(chi2 + 2 and 3 d.o.f.)   C.L.(chi2 + 3 and 3 d.o.f.)  %g \n')
         [y',  p1new'    , p2new'    ,    p3new'   ,    p4new']
  fprintf('                                                                         \n')       
 fprintf('chi-squared    (1) C.L.(chi2,2 d.o.f.) - C.L.(chi2,3 d.o.f.)   (2)C.L.(chi^2,2 d.o.f.) - (C.L.(chi^2 + 1, 3 d.o.f.)    %g \n')
         [y',      (p1-p2)'    ,   (p1new-p2new)']        
  fprintf('                                                                         \n')
  
   fprintf('mean(1)   stdv(1)    mean(2)    stdv(2)  %g \n')
   [mean(p1-p2),std(p1-p2),mean(p1new-p2new),std(p1new-p2new)]

   
figure(10)
plot(p1new,p2new,'*', p1new,p3new,'b-',p1new,p4new,'b--', p1new, p5new,'b:+')
grid on
legend('C.L.(2 d.o.f.,\chi^2) vs C.L. (3 d.o.f.,\chi^2)','C.L. (2 d.o.f., \chi^2) vs C.L.(3 d.o.f., \chi^2 * 2)','C.L. (2 d.o.f., \chi^2) vs C.L.(3 d.o.f., \chi^2 * 3)',...
    'C.L. (2 d.o.f., \chi^2) vs C.L.(3 d.o.f., \chi^2 * 3)',0);
xlabel('C.L.( 2 d.o.f.,\chi^2');
ylabel('C.L.(\chi^2 * n  and 3 d.o.f.)');
title('Confidence Levels for different \chi^2 with 2 and 3 d.o.f.');
print -deps plot_cl2_cl3_chi2_times_n

figure(11)
 semilogy(p1new,p2new,'*', p1new,p3new,'b-',p1new,p4new,'b--', p1new, p5new,'b:+')
grid on
legend('C.L.(2 d.o.f.,\chi^2) vs C.L. (3 d.o.f.,\chi^2)','C.L. (2 d.o.f., \chi^2) vs C.L.(3 d.o.f., \chi^2 * 2)','C.L. (2 d.o.f., \chi^2) vs C.L.(3 d.o.f., \chi^2 * 3)',...
    'C.L. (2 d.o.f., \chi^2) vs C.L.(3 d.o.f., \chi^2 * 3)',0);
title('semilog plot');
print -deps semilog_plot_cl2_cl3_chi2_times_n

figure(12)
y = linspace(0.0,20.5,200);
p1new = wchi2cdf(y,2); p2new = wchi2cdf(y*2,3); p3new = wchi2cdf(y*3,3); p4new = wchi2cdf(y*4,3); p5new = wchi2cdf(y*5,3); p6new = wchi2cdf(y*6,3);
plot(p1new,p2new-p1new,'*',p1new,p3new-p1new,'O',p1new,p4new-p1new,'b-',p1new,p5new-p1new,'b--',p1new,p6new-p1new,'b:+')
% xlim([5.0 20.5]);
grid on
legend('C.L.(3 d.o.f , \chi^2) - C.L.(2 d.o.f ,\chi^2)','C.L.(3 d.o.f , \chi^2 * 2) - C.L.(2 d.o.f ,\chi^2)', ...
    'C.L.(3 d.o.f , \chi^2 * 3) - C.L.(2 d.o.f ,\chi^2)','C.L.(3 d.o.f , \chi^2 * 4) - C.L.(2 d.o.f ,\chi^2)',...
    'C.L.(3 d.o.f , \chi^2 * 5) - C.L.(2 d.o.f ,\chi^2)',0);
xlabel('C.L for 2 d.o.f (\chi^2)');
ylabel('C.L. for 3 d.o.f (\chi^2 * 1,2,3,4,5,6)-C.L. for 2 d.o.f (\chi^2)');
title('Difference in C.L. vs. \chi^2 for 2 and 3 d.o.f.');
print -deps chi_chi_times_1_2_3_4_5_CL_diff_2_3dof

fprintf('C.L.(2 d.o.f, chi^2)  chi^2   C.L.(3 d.o.f., chi^2 * 2)  chi^2 * 2  C.L.(3 d.o.f., chi^2 * 3)  chi^2 * 3  C.L.(3 d.o.f., chi^2 * 4)  chi^2 * 4  C.L.(3 d.o.f., chi^2 * 5)  chi^2 * 5   %g \n')
         [p1new', y',p2new',2*y',p3new',3*y',p4new',4*y',p5new',5*y']
fprintf('                                                                         \n')       
fprintf('chi^2 chi^2*2 chi^2*3 chi^2*4 chi^2*5 C.L.(3 d.o.f., chi^2 * 2)-C.L.(2 d.o.f, chi^2)   C.L.(3 d.o.f., chi^2 * 3)-C.L.(2 d.o.f, chi^2)  C.L.(3 d.o.f., chi^2 * 4)-C.L.(2 d.o.f, chi^2)  C.L.(3 d.o.f., chi^2 * 5)-C.L.(2 d.o.f, chi^2) %g \n') 
[y',y'*2,y'*3,y'*4,y'*5,p2new'-p1new',p3new'-p1new',p4new'-p1new',p5new'-p1new']        

figure(13)
y = linspace(0.0,20.5,200);
% xlim([5 20.5]);
p1new = wchi2cdf(y,2); p2new = wchi2cdf(y,3); p3new = wchi2cdf(y*2,3); p4new = wchi2cdf(y*3,3); p5new = wchi2cdf(y*4,3); p6new = wchi2cdf(y*5,3);
plot(y,p2new-p1new,'*',y,p3new-p1new,'O',y,p4new-p1new,'b-',y,p5new-p1new,'b--')
xlim([5.0 20.0]);
grid on
legend('C.L.(3 d.o.f , \chi^2) - C.L.(2 d.o.f ,\chi^2)','C.L.(3 d.o.f , \chi^2 * 2) - C.L.(2 d.o.f ,\chi^2)', ...
    'C.L.(3 d.o.f , \chi^2 * 3) - C.L.(2 d.o.f ,\chi^2)','C.L.(3 d.o.f , \chi^2 * 4) - C.L.(2 d.o.f ,\chi^2)',0);
xlabel('\chi^2');
ylabel('C.L. for 3 d.o.f (\chi^2 * 1,2,3,4)-C.L. for 2 d.o.f (\chi^2)');
title('Difference in C.L. for different \chi^2 and 2 and 3 d.o.f.');
print -deps chi_chi_times_1_10_15_20_CL_diff_2_3dof_2

figure(14)
y = linspace(0.0,20.5,200);
% xlim([0.6 1.0]);
p1new = wchi2cdf(y,2); p2new = wchi2cdf(y,3); p3new = wchi2cdf(y*10,3); p4new = wchi2cdf(y*15,3); p5new = wchi2cdf(y*20,3); p6new = wchi2cdf(y*25,3);
plot(y,p2new-p1new,'*',y,p3new-p1new,'O',y,p4new-p1new,'b-',y,p5new-p1new,'b--')
xlim([5.0 20.0]);
grid on
legend('C.L.(3 d.o.f , \chi^2) - C.L.(2 d.o.f ,\chi^2)','C.L.(3 d.o.f , \chi^2 * 10) - C.L.(2 d.o.f ,\chi^2)', ...
    'C.L.(3 d.o.f , \chi^2 * 15) - C.L.(2 d.o.f ,\chi^2)','C.L.(3 d.o.f , \chi^2 * 20) - C.L.(2 d.o.f ,\chi^2)',0);
xlabel('\chi^2');
ylabel('C.L. for 3 d.o.f (\chi^2 * 1,10,15,20)-C.L. for 2 d.o.f (\chi^2)');
title('Difference in C.L. for different \chi^2 and 2 and 3 d.o.f.');
print -deps chi_chi_times_1_10_15_20_CL_diff_2_3dof_2


fprintf('C.L.(2 d.o.f, chi^2)  chi^2   C.L.(3 d.o.f., chi^2 * 10)  chi^2 * 10  C.L.(3 d.o.f., chi^2 * 15)  chi^2 * 15  C.L.(3 d.o.f., chi^2 * 20)  chi^2 * 20  C.L.(3 d.o.f., chi^2 * 25)  chi^2 * 25   %g \n')
         [p1new', y',p2new',2*y',p3new',3*y',p4new',4*y',p5new',5*y']
fprintf('                                                                         \n')       
fprintf('chi^2 chi^2*2 chi^2*3 chi^2*4 chi^2*5 C.L.(3 d.o.f., chi^2 * 10)-C.L.(2 d.o.f, chi^2)   C.L.(3 d.o.f., chi^2 * 15)-C.L.(2 d.o.f, chi^2)  C.L.(3 d.o.f., chi^2 * 20)-C.L.(2 d.o.f, chi^2)  C.L.(3 d.o.f., chi^2 * 25)-C.L.(2 d.o.f, chi^2) %g \n') 
[y',y'*2,y'*3,y'*4,y'*5,p2new'-p1new',p3new'-p1new',p4new'-p1new',p5new'-p1new']        


% xlabel('Confidence Level, n degrees of freedom');
% ylabel('Confidence Level, n+1 degrees of freedom');
% title('Confidence Level for n+1 d.o.f. vs. Confidence Level for n d.o.f., where n = 1,2,3,4,5,6');
% print -deps plot_cl_n_n+1

   
   