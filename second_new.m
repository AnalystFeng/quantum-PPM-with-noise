clear;
clc;
close all;
%% 挤压真空态向量
r = 2;
theta = pi;
% aimuda = tanh(r)*exp(1i*theta);
% wei = 20;
% ket = zeros(wei,1);
% for n = 0:2:wei/2
%     ket(n+1) = aimuda^n*sqrt(sech(r)) * sqrt(factorial(2*n))/(2^n*factorial(n));
% end
%% 位移挤压态
Wei = 100;
alpha = 4;

ket_sq = zeros(Wei,1);

for n = 0:1:Wei-1
    ket_sq(n+1) = exp(-0.5*alpha^2+alpha^2*exp(1i*theta)*tanh(r)/2)*...
        (1i)^n*sqrt(exp(1i*n*theta/2)/(factorial(n)*cosh(r)))*...
        (tanh(r)/2)^(n/2)*...
        hermiteH(n,-0.5*1i*exp(-1i*theta/2)*sqrt(2/tanh(r))*(alpha-exp(1i*theta)*tanh(r)*alpha));
end % hermiteH(i-1,mu*v/(B^2))
plot(0:1:n,abs(ket_sq))
title('挤压态包含n个光子数的概率')









%% H 函数
function H = H(n,x)
H = 0;
for j = 0:1:floor(n/2)
    H = H+x^j/( factorial(n-2*j)*factorial(j) );
end
end