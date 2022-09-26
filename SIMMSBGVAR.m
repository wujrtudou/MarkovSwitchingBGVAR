
%符佳丽毕业论文模拟数据的代码，也就是生成模拟数据的代码 -- final 


function [X, B_p, B_0,S] = SIMMSBGVAR(N, T, lag)

rng(123);

T0 = 5e2;  n = 5;    k = round(N/n);
T = T + T0;


B01 = zeros(5);  
B01(1,3) = -0.8;  
B01(4,2) = 0.5;   
% B01(4,5) =-0.5;


B02 = zeros(5);  
B02(1,2) = 0.8;  
B02(1,5) =-0.9;  
% B02(2,3) =0.8;




B11 = zeros(5);  
B11(1,1) =-0.8;  
B11(2,1) =0.6;
B11(3,1)=0.7; 
B11(5,2) = -0.6; 
B11(2,3) = 0.5;
B11(3,3) = -0.5;
B11(4,3) = 0.5;
B11(4,4) = 0.7;
B11(5,5) = 0.6;



B12 = zeros(5);  
B12(1,1) =0.6;  
B12(1,4) =0.8;
B12(2,1)=-0.9; 
% B12(3,3) = 0.6; 
B12(3,3) = 0.7;
% B12(4,1) = 0.6;
B12(4,2) = 0.6;
% B12(5,1) =0.6;
B12(5,3) = -0.7;
 



B21 = zeros(5);  
% B21(2,1) = 0.5;  
% B21(4,1) = -0.6;  
B21(3,5) = -0.4;
B21(5,3) = -0.4;

B22 = zeros(5);  
B22(4,3) = 0.2; 
% B22(1,5) = -0.3;
% 
B31 = zeros(5);  
% B31(3,4) = -0.4; 
B31(5,2) = 0.4; 
% 
B32 = zeros(5);  
% B32(2,3) = 0.1; 

%U是N列T行的，每行数据都是来自标准正态分布的随机数
U  = gennormrnd(0,1,N,T);
B_01 = kron(eye(k),B01);
B_02 = kron(eye(k),B02);

B_0(:,:,1)=B_01;
B_0(:,:,2)=B_02;


B_11 = kron(eye(k),B11);
B_12 = kron(eye(k),B12);

B_21 = kron(eye(k),B21);
B_22 = kron(eye(k),B22);
% 
% 
B_31 = kron(eye(k),B31);
B_32 = kron(eye(k),B32);

X  = gennormrnd(0,1,N,T);
%U是的T行N列，每行数据都是来自标准正态分布的随机数，均值为0，标准差为1

a = .05;
b = .05;

S = zeros(T,1); S(1) = (rand < .5);
for t=2:T
    if S(t-1) == 0
        S(t) = (rand < a);
    else
        S(t) = 1 - (rand < b);
    end
end


if lag == 0
    for t=lag+1:T  %lag+1:T相当于2：T
        if S(t)==0
            
            X(:,t) = B_01*X(:,t) + U(:,t);
        else
            X(:,t) = B_02*X(:,t) + U(:,t);
        end
    end    
    X = X(:,T0+1:end)';
    B_p = 0.*B_1;    

elseif lag == 1
    for t=lag+1:T 
        if S(t)==0
            X(:,t) = B_01*X(:,t) + B_11*X(:,t-1) +U(:,t);
        else
            X(:,t) = B_02*X(:,t) + B_12 *X(:,t-1) + 10*U(:,t);
            %X(:,t) = B_02*X(:,t) + B_12 *X(:,t-1) + U(:,t); %这种设置也是可以运行的，只是识别效果没有10好
        end
    end    
    X = X(:,T0+1:end)';%我终于懂了，因为这里有一个转置的符号
    B_p(:,:,1) = B_11;
    B_p(:,:,2) = B_12;

elseif lag == 2
    for t=lag+1:T
        if S(t)==0
            X(:,t) = B_01*X(:,t) + B_11*X(:,t-1) + B_21*X(:,t-2) + U(:,t);
        else
            %X(:,t) = B_02*X(:,t) + B_12*X(:,t-1) + B_22*X(:,t-2) + 10*U(:,t);
            X(:,t) = B_02*X(:,t) + B_12*X(:,t-1) + B_22*X(:,t-2) + U(:,t);
        end
    end  
    X = X(:,T0+1:end)';
    B_p(:,:,1) = [B_11,B_21];
    B_p(:,:,2) = [B_12,B_22];
    
elseif lag == 3
    for t=lag+1:T
        if S(t)==0
            X(:,t) = 5*B_01*X(:,t) + 3*B_11*X(:,t-1) +2* B_21*X(:,t-2) + ...
            3*B_31*X(:,t-3) + U(:,t);
            
        else
            X(:,t) = B_02*X(:,t) + B_12*X(:,t-1) + B_22*X(:,t-2) + ...
            B_32*X(:,t-3) + U(:,t);
        end
    end  
    
    X = X(:,T0+1:end)';
    B_p(:,:,1) = [B_11,B_21,B_31];
    B_p(:,:,2) = [B_12,B_22,B_32];
    
end
S = S(T0+1:end,1);

function r = gennormrnd(mu,sigma,m,n)
r = randn(m,n) .* sigma + mu;