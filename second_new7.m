clc;clear;close all;
C=[1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10];
B=[0.001 0.071 0.131 0.166 0.156 0.106 0.081 0.051 0.041 0.026 0.021 0.011 0.010 0.009 0.008 0.007 0.006 0.005 0.004];
r=asinh(sqrt(B.*C));
a=sqrt(C-(sinh(r)).^2);
subplot(2,1,1)
plot(C,r,'--*k','Linewidth',2)
axis([1 6.5 0 0.7])
xlabel({'Channel energy C';'(a)'})
ylabel('Optimal squeezing parameter r')
set(gca,'fontsize',15)
subplot(2,1,2)
plot(C,B,'--*k','Linewidth',2)
axis([1 6.5 0 0.18])
xlabel({'Channel energy C';'(b)'})
ylabel('Optimal squeezing fraction β')
set(gca,'fontsize',15)