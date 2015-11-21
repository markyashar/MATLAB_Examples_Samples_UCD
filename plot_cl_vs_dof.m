figure(10)

n = 1:50;
x=2.30;
y=6.17;
z=9.21;
p0 = wchi2cdf(x,n);
p1 = wchi2cdf(y,n);
p2 = wchi2cdf(z,n);
% p1 = wchi2cdf(x,2); p2 = wchi2cdf(x,3); p3 = wchi2cdf(x,4); p4 = wchi2cdf(x,5); p5 = wchi2cdf(x,6)
% plot(x,p0,'*', x,p1,'b-',x,p2,'b--', x, p3,'b:+',x, p4,':ks')
plot(n,p0,'*',n,p1,'b-',n,p2,'b:+');
grid on
legend('\chi^2 = 2.30',',chi^2 = 6.17','chi^2 = 9.21',0);
xlabel('n = number of degrees of freedom');
ylabel('Confidence Level');
title('Confidence Level vs. number of d.o.f. for \chi^2 = 2.30, 6.17, 9.21');