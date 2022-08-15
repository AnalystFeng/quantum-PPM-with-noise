%% 噪声光子分布
clc;close all;
clear;
%%  不可更改
%%  不可更改
%%  不可更改
%%

Wei = 20;
% noise = eye(Wei);
theta = pi;
% N = 0.1;
ket_0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
% ket_noi = zeros(Wei,1);
% % 定义 noise
% for i =1:1:Wei
%     noise(i,i) = N^i/(1+N)^(1+i);
%     ket_noi(i) =  sqrt(noise(i,i));
% end
% ket_noi = ket_noi/norm(ket_noi);
% subplot(1,3,1)
% bar(0:1:19,abs(ket_noi).^2,'k')
% axis([-0.6 20 0 1])
% ylabel('Probability')
% xlabel({'Photon number N';'(a)'})
% set(gca,'fontsize',15)
% %%
% noise = eye(Wei);
% theta = pi;
% N = 0.2;
% ket_0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
% ket_noi = zeros(Wei,1);
% % 定义 noise
% for i =1:1:Wei
%     noise(i,i) = N^i/(1+N)^(1+i);
%     ket_noi(i) =  sqrt(noise(i,i));
% end
% ket_noi = ket_noi/norm(ket_noi);
% subplot(1,3,2)
% bar(0:1:19,abs(ket_noi).^2,'k')
% axis([-0.6 20 0 1])
% 
% xlabel({'Photon number N';'(b)'})
% set(gca,'fontsize',15)
%%
% noise = eye(Wei);
% theta = pi;
% N = 0.3;
% ket_0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
% ket_noi = zeros(Wei,1);
% % 定义 noise
% for i =1:1:Wei
%     noise(i,i) = N^i/(1+N)^(1+i);
%     ket_noi(i) =  sqrt(noise(i,i));
% end
% ket_noi = ket_noi/norm(ket_noi);
% subplot(1,3,3)
% bar(0:1:19,abs(ket_noi).^2,'k')
% axis([-0.6 20 0 1])
% 
% xlabel({'Photon number N';'(c)'})
% set(gca,'fontsize',15)
% ket_0_noi = (ket_0+ket_noi)/norm(ket_0+ket_noi);
%% 噪声光子分布
%%
%%
%%
%% 相干态 挤压态
% alpha = 2;
% [ket_co,rou_coh] = rou_co(alpha);
% figure(1)
% subplot(1,2,1)
% bar(0:1:19,(abs(ket_0)).^2,'k')
% axis([-0.6 20 0 1])
% xlabel({'Photon number N';'(a)'})
% ylabel('Probability')
% set(gca,'fontsize',15)
% 
% subplot(1,2,2)
% bar(0:1:19,(abs(ket_co)).^2,'k')
% axis([-0.6 20 0 1])
% xlabel({'Photon number N';'(b)'})
% set(gca,'fontsize',15)
% 
% r = 0.1;
% theta = pi;
% [ket_squ,rou_squ] = rou_sq(r,alpha,theta);
% 
% figure(2)
% subplot(1,2,1)
% bar(0:1:19,(abs(ket_co)).^2,'k')
% axis([-0.6 20 0 1])
% xlabel({'Photon number N';'(a)'})
% ylabel('Probability')
% set(gca,'fontsize',15)
% 
% subplot(1,2,2)
% bar(0:1:19,(abs(ket_squ)).^2,'k')
% axis([-0.6 20 0 1])
% xlabel({'Photon number N';'(b)'})
% set(gca,'fontsize',15)
%%
%%
%%
%%
Wei = 20;
noise = eye(Wei);
theta = pi;
N = 0.2;
ket_0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
ket_noi = zeros(Wei,1);
% 定义 noise
for i =1:1:Wei
    noise(i,i) = N^i/(1+N)^(1+i);
    ket_noi(i) =  sqrt(noise(i,i));
end
ket_noi = ket_noi/norm(ket_noi);
alpha = 2;
[ket_co,rou_coh] = rou_co(alpha);
subplot(2,3,1)
bar(0:1:19,(4*(ket_co.^2)+0.2*(ket_noi.^2))/4.2,'m')
axis([-0.6 20 0 0.25])
xlabel({'Photon number N';'(a)'})
ylabel('Probability')
legend('coherent state')
set(gca,'fontsize',15)


N = 0.4;
ket_0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
ket_noi = zeros(Wei,1);
% 定义 noise
for i =1:1:Wei
    noise(i,i) = N^i/(1+N)^(1+i);
    ket_noi(i) =  sqrt(noise(i,i));
end
ket_noi = ket_noi/norm(ket_noi);
alpha = 2;
[ket_co,rou_coh] = rou_co(alpha);
subplot(2,3,2)
bar(0:1:19,(4*(ket_co.^2)+0.4*(ket_noi.^2))/4.4,'m')
axis([-0.6 20 0 0.25])
xlabel({'Photon number N';'(b)'})

set(gca,'fontsize',15)

N = 0.6;
ket_0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
ket_noi = zeros(Wei,1);
% 定义 noise
for i =1:1:Wei
    noise(i,i) = N^i/(1+N)^(1+i);
    ket_noi(i) =  sqrt(noise(i,i));
end
ket_noi = ket_noi/norm(ket_noi);
alpha = 2;
[ket_co,rou_coh] = rou_co(alpha);
subplot(2,3,3)
bar(0:1:19,(4*(ket_co.^2)+0.6*(ket_noi.^2))/4.6,'m')
axis([-0.6 20 0 0.25])
xlabel({'Photon number N';'(c)'})

set(gca,'fontsize',15)

%%
%%
%%
%% 挤压态加噪声的光子数分布
r = 0.3;
Wei = 20;
noise = eye(Wei);
theta = pi;
N = 0.2;
ket_0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
ket_noi = zeros(Wei,1);
% 定义 noise
for i =1:1:Wei
    noise(i,i) = N^i/(1+N)^(1+i);
    ket_noi(i) =  sqrt(noise(i,i));
end
ket_noi = ket_noi/norm(ket_noi);
alpha = 2;

[ket_squ,rou_squ] = rou_sq(r,alpha,theta);
subplot(2,3,4)
bar(0:1:19,(4*(abs(ket_squ).^2)+0.2*(ket_noi.^2))/4.2,'g')
axis([-0.6 20 0 0.25])
xlabel({'Photon number N';'(d)'})
ylabel('Probability')
legend('squeezed state')
set(gca,'fontsize',15)


N = 0.4;
ket_0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
ket_noi = zeros(Wei,1);
% 定义 noise
for i =1:1:Wei
    noise(i,i) = N^i/(1+N)^(1+i);
    ket_noi(i) =  sqrt(noise(i,i));
end
ket_noi = ket_noi/norm(ket_noi);
alpha = 2;
[ket_squ,rou_squ] = rou_sq(r,alpha,theta);
subplot(2,3,5)
bar(0:1:19,(4*(abs(ket_squ).^2)+0.4*(ket_noi.^2))/4.4,'g')
axis([-0.6 20 0 0.25])
xlabel({'Photon number N';'(e)'})

set(gca,'fontsize',15)

N = 0.6;
ket_0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
ket_noi = zeros(Wei,1);
% 定义 noise
for i =1:1:Wei
    noise(i,i) = N^i/(1+N)^(1+i);
    ket_noi(i) =  sqrt(noise(i,i));
end
ket_noi = ket_noi/norm(ket_noi);
alpha = 2;
[ket_squ,rou_squ] = rou_sq(r,alpha,theta);
subplot(2,3,6)
bar(0:1:19,(4*(abs(ket_squ).^2)+0.6*(ket_noi.^2))/4.6,'g')
axis([-0.6 20 0 0.25])
xlabel({'Photon number N';'(f)'})

set(gca,'fontsize',15)




%% 生成相干态，参数为 α
function [ket_co,rou_co] = rou_co(alpha)
Wei = 20;
ket_co = zeros(Wei,1);
for n = 1:1:Wei
    ket_co(n) = exp((-1/2)*alpha^2)   *alpha^(n-1)/(factorial(n-1))^0.5;
end
% ket_co = ket_co/norm(ket_co);
rou_co = ket_co*ket_co';

end
%% 生成一个纯挤压态的密度算子，参数为 r,α，θ
function [ket_sq,rou_sq] = rou_sq(r,alpha,theta)
Wei = 20;
ket_sq = zeros(Wei,1);

for n = 0:1:Wei-1
    ket_sq(n+1) = exp(-0.5*alpha^2+alpha^2*exp(1i*theta)*tanh(r)/2)*...
        (1i)^n*sqrt(exp(1i*n*theta/2)/(factorial(n)*cosh(r)))*...
        (tanh(r)/2)^(n/2)*...
        hermiteH(n,-0.5*1i*exp(-1i*theta/2)*sqrt(2/tanh(r))*(alpha-exp(1i*theta)*tanh(r)*alpha));
end % hermiteH(i-1,mu*v/(B^2))
% ket_sq = ket_sq/norm(ket_sq);

rou_sq = ket_sq*ket_sq';
rou_sq = rou_sq/trace(rou_sq); % 使挤压态的迹归一化，可能有问题
end