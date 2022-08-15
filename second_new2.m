clc;close all;
clear;
%%  不可更改
%%  不可更改
%%  不可更改
%%

Wei = 20;
noise = eye(Wei);
theta = pi;
N = 0.1;
ket_0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
ket_noi = zeros(Wei,1);
% 定义 noise
for i =1:1:Wei
    noise(i,i) = N^i/(1+N)^(1+i);
    ket_noi(i) =  sqrt(noise(i,i));
end
ket_0_noi = (ket_0+ket_noi)/norm(ket_0+ket_noi);
%% coherent state
alpha_list = 0.001:0.2:2.001;
Pc_coh = zeros(1,length(alpha_list));
Pc_coh_noi = zeros(1,length(alpha_list));
for i = 1:1:length(alpha_list)
    alpha = alpha_list(i); % α = ∆
    [ket_coh,rou_coh] = rou_co(alpha) ;  % r=0,,迹为1，随着r的增加，迹增大
    gamma_0 = kron(ket_0,kron(ket_0,ket_coh));
    gamma_1 = kron(ket_0,kron(ket_coh,ket_0));
    gamma_2 = kron(ket_coh,kron(ket_0,ket_0));
    tao_co = [gamma_0 gamma_1 gamma_2];
    % T = tao_co*tao_co';
    G = tao_co'*tao_co;
    X = G^(1/2);
    diagg = diag(X);
    Pc(i) = 0.3333*sum(diagg.^2)
    %%
    ket_coh = (ket_coh+ket_noi)/norm(ket_coh+ket_noi);
    gamma_0 = kron(ket_0_noi,kron(ket_0_noi,ket_coh));
    gamma_1 = kron(ket_0_noi,kron(ket_coh,ket_0_noi));
    gamma_2 = kron(ket_coh,kron(ket_0_noi,ket_0_noi));
    tao_co = [gamma_0 gamma_1 gamma_2];
    % T = tao_co*tao_co';
    G = tao_co'*tao_co;
    X = G^(1/2);
    diagg = diag(X);
    Pc_n(i) = 0.3333*sum(diagg.^2)
end
semilogy(alpha_list,1-Pc,'--*b','Linewidth',2)
hold on
semilogy(alpha_list,1-Pc_n,'-db','Linewidth',2)
% legend('pure coherent state','mixed coherent state')
hold on

%% squeezed state
r = 0.2;
Pc_sq = zeros(length(r),length(alpha_list));
Pc_sq_noi = zeros(length(r),length(alpha_list));
for j = 1:1:length(r)
    for i = 1:1:length(alpha_list)
        alpha = alpha_list(i); % α = ∆
        [ket_sq,rou_squ] = rou_sq(r(j),alpha,theta) ;  % r=0,,迹为1，随着r的增加，迹增大
        gamma_sq_0 = kron(ket_0,kron(ket_0,ket_sq));
        gamma_sq_1 = kron(ket_0,kron(ket_sq,ket_0));
        gamma_sq_2 = kron(ket_sq,kron(ket_0,ket_0));
        tao_sq = [gamma_sq_0 gamma_sq_1 gamma_sq_2];
        % T_sq = tao_sq*tao_sq';
        G = tao_sq'*tao_sq;
        X = G^(1/2);
        diagg = diag(X);
        Pc(j,i) = 0.3333*sum(diagg.^2)
%%
        gamma_sq_0 = kron(ket_0_noi,kron(ket_0_noi,ket_sq));
        gamma_sq_1 = kron(ket_0_noi,kron(ket_sq,ket_0_noi));
        gamma_sq_2 = kron(ket_sq,kron(ket_0_noi,ket_0_noi));
        tao_sq = [gamma_sq_0 gamma_sq_1 gamma_sq_2];
        G = tao_sq'*tao_sq;
        X = G^(1/2);
        diagg = diag(X);
        Pc_n(j,i) = real(0.3333*sum(diagg.^2))


    end
end

semilogy(alpha_list,1-Pc,'--pr','Linewidth',2)
hold on
semilogy(alpha_list,1-Pc_n,'-+r','Linewidth',2)
legend('pure coherent state','coherent state with noise','pure squeezed state','squeezed state with noise')
axis([0 2.01 0.0001 1])
ylabel('Error probability Pe')
xlabel('|α|')
grid ;
%% 


%% 生成一个纯相干态的密度算子，参数为α
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
%% 密度算子加噪
function rou_noi = rou_noi(rou,np) % np : noise parameter:σ
Wei = 20;
for n = 1:1:Wei
    for m = 1:1:Wei
        rou(n,m) = rou(n,m)*exp(-(n-m)^2*np^4/4);
    end
end
% rou_noi = rou/trace(rou);
rou_noi = rou;
end

%% H 函数
function H = H(n,x)
H = 0;
for j = 0:1:floor(n/2)
    H = H+x^j/( factorial(n-2*j)*factorial(j) );
end
end