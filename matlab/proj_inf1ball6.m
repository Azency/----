function [X, FX] = proj_inf1ball6(B, tau, max_iter, delt)
% project the matrix B onto the l-1-inf ball of radius tau, i.e.
% X = argmin |X-B|_F
% s.t. sum(max(abs(B), [] ,2)) <= tau
% max_iter=100, delt=1e-8 by default
if nargin==2
    max_iter = 100;
    delt = 1e-8;
end
if norm_inf_1(B)-tau<eps
    X = B;
    FX = 0;
    return;
end
[m,n] = size(B); %获取矩阵B的维度m,n
x_inf = zeros(m, 1); % 解矩阵X的每行的绝对值最大值组成的列向量
kk = zeros(m, 1); 
r = n; 
mu = 0;
z = zeros(m,1);
FX = zeros(max_iter, 1); 
G = ones(m,1); 
sB = sign(B); %矩阵B的符号矩阵
aB = abs(B); %矩阵B的绝对值矩阵，下文全用绝对值矩阵计算，最后再乘以符号矩阵，和矩阵B符号保持一致
saB = sum(aB, 2); % 列向量，对aB矩阵每行求和生成
B2 = zeros(m, n);
B3 = zeros(m, n);
aBt = aB'; % 矩阵aB的转置，将行向量的计算转成列向量
B2t = B2'; 
for ii=1:m
    B2t(:,ii) = sort(aBt(:,ii), 'ascend');
    % 对aBt每列（原aB矩阵的每行）做升序排序
end
B2 = B2t'; % 矩阵aB每行升序排序后的矩阵

%% 计算引理1的b_k^* 序列，具体公式见论文

B3(:,1) = - saB; % B3矩阵第一列，aB矩阵每行和的相反数
saB = saB - B2(:,1); % 列向量，将矩阵每行和减去矩阵每行最小值
for k=2:n
    B3(:,k) = - saB + (n+1-k+r)*B2(:,k-1);
    saB = saB - B2(:,k);
end
B3 = B3/r; % 计算结束，B3的每行为B2的每行的b*向量


B2t = B2'; % 转置，将升序矩阵B2的行向量的计算转成列向量
B3t = B3'; % 转置，将b*矩阵B3的行向量的计算转成列向量

for iter=1:max_iter
    % update X by scattering Algorithm 1
    for ii=1:m
        if G(ii)
            [x_inf(ii), kk(ii)] = prox_inf5(B2t(:,ii), B3t(:,ii), z(ii)-mu/r, r, kk(ii));
            % 函数prox_inf5 为计算Lemma1的结论的子程序            
            % x_inf(ii)为Lemma1的x_n，kk(ii)为Lemma1中b_k^*< c <b_k+1^*里k的值，
            % 对行ii=1:m,得到向量x_inf(解向量最大值)和向量kk
        end
    end
    
    %% 对于解矩阵x_inf， 删去其中全为0的矩阵，使其不参与接下来的计算
    f = x_inf;
%     f = max(abs(X), [], 2);
    G = f>eps;
    m_prime = sum(G);
    
    %% Algorithm 2 计算mu和z，且只计算上文矩阵中的非零行
    % update mu
    mu = mu + r/m_prime*(sum(f) - tau);
    % update z   
    z = f - 1/m_prime*(sum(f) - tau);
    FX(iter) = abs(sum(f) - tau); %% 目标函数误差值
      
%% 停止条件，误差值小于delt，或者误差序列前后两个值不变
  if iter>1
    eFX = 1e8.*FX(iter-1) - 1e8.*FX(iter);
    if FX(iter)<delt || eFX<delt 
        break;
    end
    
  end
      
end



FX = FX(1:iter);
X = aB; % 绝对值矩阵aB
for ii=1:m
    % 按行计算
    xx = X(ii,:);
    xx(xx>x_inf(ii)) = x_inf(ii); % 将aB中大于x_inf的值全替换成x_inf
    X(ii,:) = xx;
end
X = sB.*X; % 将解矩阵与原B矩阵符号一致