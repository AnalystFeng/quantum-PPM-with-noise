clc;close all;
clear;
r = 0.9;
theta = pi;
alpha = 5; 
Wei = 10;
ket_squeezed = zeros(Wei,1);
%%  纯挤压态密度算子
mu = cosh (r);
v = sinh(r)*exp(1*i*theta);
B = mu*alpha-v*alpha';
for i = 1:1:Wei
    ket_squeezed(i) = (sqrt(factorial(i))/mu) * (B/mu)^i * hermiteH(i,mu*v/(B^2)) * exp(  -0.5*(abs(B))^2 - B^2*(v'/2*mu));
end
rou_squeezed_pure = ket_squeezed*ket_squeezed';
rou_squeezed_pure = rou_squeezed_pure/trace(rou_squeezed_pure);  % rou 挤压态的密度算子  ；
X = 1 : 1 : Wei;
Y = 1 : 1 : Wei;
figure(1)
mesh(X, Y, abs(rou_squeezed_pure))
title('纯挤压态密度算子')

sita = 0.2;
deta = sita^2/2;
rou_squeezed_noi = zeros(Wei,Wei);
for n = 1:Wei
    for m = 1:Wei
        rou_squeezed_noi(n,m) = rou_squeezed_pure(n,m)*exp(  -(n-m)^2*deta^2     );
    end
end
figure(2)
mesh(X, Y, abs(rou_squeezed_noi))
title('混合挤压态密度算子')
rou0 = zeros(Wei,Wei);
rou0(1,1) = 1;
P = 0.5*(   1-0.5*norm(rou_squeezed_noi-rou0))



%% 混合态挤压态密度算子
% load rou_noi.mat
% load rou_squeezed_pure
% rou_squ_noi = rou_squeezed_pure+rou_noi;
% rou_squ_noi = rou_squ_noi/trace(rou_squ_noi);
% figure(2)
% mesh(X, Y, abs(rou_squ_noi))
% title('混合态挤压态密度算子')
% save rou_squ_noi
% %%
% [A,E]=eig(rou_squ_noi-rou_noi);
% q = 0;
% for i = 1:1:20
%     if E(i,i)>0
%         q=q+E(i,i);
%     end
% end
% Pe = 0.5*(1-q)
%%
% h=7;
% R_deta = rou_squ_noi;
% R_0 = rou_noi;
% [A,E] = eig(R_deta);
% [C,D] = eig(R_0);
% rou_up_1 = A(:,end-h+1:end)*E(end-h+1:end,end-h+1:end)*A(:,end-h+1:end)';
% rou_up_0 = C(:,end-h+1:end)*D(end-h+1:end,end-h+1:end)*C(:,end-h+1:end)';
% rou_down_0 = kron(rou_up_0,kron(rou_up_0,rou_up_1));
% 
% gamma_up_1 = A(:,end-h+1:end)*sqrt(E(end-h+1:end,end-h+1:end));
% gamma_up_0 = C(:,end-h+1:end)*sqrt(D(end-h+1:end,end-h+1:end));
% gamma_down_0 = kron(gamma_up_0,kron(gamma_up_0,gamma_up_1));
% gamma_down_1 = kron(gamma_up_0,kron(gamma_up_1,gamma_up_0));
% gamma_down_2 = kron(gamma_up_1,kron(gamma_up_0,gamma_up_1));
% States_matrix = [gamma_down_0 gamma_down_1 gamma_down_2];
% T = States_matrix*States_matrix';
% % G = States_matrix'*States_matrix;
% % S = (rou_down_0*T^(-1/2))^2;
% Pc = trace(  (rou_down_0*     pinv(T^(0.5))    )^2           )