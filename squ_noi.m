clc;close all;
clear;
Wei=20;
N=0.7; % 噪声光子数
I = eye(Wei);
P = N/7;
%%
rou_alpha0 = P*(I(:,1)*I(:,1)');
%%
alpha1 = 1;
ket_alpha_1 = zeros(Wei,1);
for n = 1:1:Wei
    ket_alpha_1(n) = exp(-1/2)*alpha1^(n-1)/sqrt(factorial(n-1));
end
rou_alpha1 = ket_alpha_1*ket_alpha_1';
%%
alpha2 = 1/2+1i*sqrt(3)/2;
ket_alpha_2 = zeros(Wei,1);
for n = 1:1:Wei
    ket_alpha_2(n) = exp(-1/2)*alpha2^(n-1)/sqrt(factorial(n-1));
end
rou_alpha2 = ket_alpha_2*ket_alpha_2';
%%
alpha3 = -1/2+1i*sqrt(3)/2;
ket_alpha_3 = zeros(Wei,1);
for n = 1:1:Wei
    ket_alpha_3(n) = exp(-1/2)*alpha3^(n-1)/sqrt(factorial(n-1));
end
rou_alpha3 = ket_alpha_3*ket_alpha_3';
%%
alpha4 = -1;
ket_alpha_4 = zeros(Wei,1);
for n = 1:1:Wei
    ket_alpha_4(n) = exp(-1/2)*alpha4^(n-1)/sqrt(factorial(n-1));
end
rou_alpha4 = ket_alpha_4*ket_alpha_4';
%%
alpha5 = -1/2-1i*sqrt(3)/2;
ket_alpha_5 = zeros(Wei,1);
for n = 1:1:Wei
    ket_alpha_5(n) = exp(-1/2)*alpha5^(n-1)/sqrt(factorial(n-1));
end
rou_alpha5 = ket_alpha_5*ket_alpha_5';
%%
alpha6 = 1/2-1i*sqrt(3)/2;
ket_alpha_6 = zeros(Wei,1);
for n = 1:1:Wei
    ket_alpha_6(n) = exp(-1/2)*alpha6^(n-1)/sqrt(factorial(n-1));
end
rou_alpha6 = ket_alpha_6*ket_alpha_6';
%%近似噪声密度算子
rou_noi = P*(rou_alpha0+rou_alpha1+rou_alpha2+rou_alpha3+rou_alpha4+rou_alpha5+rou_alpha6);
rou_noi = rou_noi/trace(rou_noi);   % 迹归一化
%%
X = 1 : 1 : Wei;
Y = 1 : 1 : Wei;
figure(1)
mesh(X, Y, rou_noi)
title('近似噪声密度算子')
%%
rou_noi_real = zeros(Wei);
for i=1:1:Wei
    rou_noi_real(i,i) = N^i/(1+N)^(i+1);
end
rou_noi_real = rou_noi_real/trace(rou_noi_real);
figure(2)
mesh(X, Y, rou_noi_real)
title('噪声密度算子')
%%
r = 0.5;
alpha = 8;
theta = pi;
mu = cosh (r);
v = sinh(r)*exp(1*i*theta);
B = mu*alpha-v*alpha';
ket_squeezed = zeros(Wei,1);
for i = 1:1:Wei
    ket_squeezed(i) = (sqrt(factorial(i))/mu) * (B/mu)^i * hermiteH(i,mu*v/(B^2)) * exp(  -0.5*(abs(B))^2 - B^2*(v'/2*mu));
end
rou_squeezed_pure = ket_squeezed*ket_squeezed';
rou_squeezed_pure = rou_squeezed_pure/trace(rou_squeezed_pure);  % 纯挤压态的密度算子  ； 
figure(3)
mesh(X, Y,rou_squeezed_pure)
title('纯挤压态的密度算子')
%%
