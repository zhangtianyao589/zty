function [u, u_hat, omega] = VMD(signal, alpha, tau, K, DC, init, tol)
% Variational Mode Decomposition
% Input and Parameters:
% ---------------------
% signal  - the time domain signal (1D) to be decomposed
% alpha   - the balancing parameter of the data-fidelity constraint 惩罚因子
% tau     - time-step of the dual ascent ( pick 0 for noise-slack )
% K       - the number of modes to be recovered 模态分量
% DC      - true if the first mode is put and kept at DC (0-freq)
% init    - 0 = all omegas start at 0
%                    1 = all omegas start uniformly distributed
%                    2 = all omegas initialized randomly
% tol     - tolerance of convergence criterion; typically around 1e-6
%
% Output:
% -------
% u       - the collection of decomposed modes
% u_hat   - spectra of the modes
% omega   - estimated mode center-frequencies
 
 
 
%---------- Preparations
% 输入信号的采样周期和频率
save_T = length(signal);
fs = 1/save_T;
 
% extend the signal by mirroring  通过镜像扩展信号
T = save_T;
f_mirror(1:T/2) = signal(T/2:-1:1);
f_mirror(T/2+1:3*T/2) = signal;
f_mirror(3*T/2+1:2*T) = signal(T:-1:T/2+1);
f = f_mirror;
 
% Time Domain 0 to T (of mirrored signal) 时域 0-T
T = length(f);
t = (1:T)/T;
 
% Spectral Domain discretization 谱域离散化
freqs = t-0.5-1/T;
 
% Maximum number of iterations最大迭代次数 (if not converged yet, then it won't anyway)
N = 500;

% For future generalizations: individual alpha for each mode  每个模式的单个alpha
Alpha = alpha*ones(1,K);

% Construct and center f_hat  建造并居中
f_hat = fftshift((fft(f)));
f_hat_plus = f_hat;
f_hat_plus(1:T/2) = 0;

% matrix keeping track of every iterant 每个迭代的矩阵保持跟踪// could be discarded for mem
u_hat_plus = zeros(N, length(freqs), K);

% Initialization of omega_k 初始化omega_k
omega_plus = zeros(N, K);
switch init
    case 1
        for i = 1:K
            omega_plus(1,i) = (0.5/K)*(i-1);
        end
    case 2
        omega_plus(1,:) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,K)));
    otherwise
        omega_plus(1,:) = 0;
end
% if DC mode imposed, set its omega to 0，直流
if DC
    omega_plus(1,1) = 0;
end
% start with empty dual variables  从空的对偶变量开始
lambda_hat = zeros(N, length(freqs));
% other inits
uDiff = tol+eps; % update step
n = 1; % loop counter
sum_uk = 0; % accumulator 
% ----------- Main loop for iterative updates  迭代更新的主循环
while ( uDiff > tol &&  n < N ) % not converged and below iterations limit 未收敛且低于迭代限制
    
    % update first mode accumulator，累加器
    k = 1;
    sum_uk = u_hat_plus(n,:,K) + sum_uk - u_hat_plus(n,:,1);
    
    % update spectrum of first mode through Wiener filter of residuals，利用残差维纳滤波器更新第一模谱
    u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
    
    % update first omega if not held at 0
    if ~DC
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
    end
    
    % update of any other mode
    for k=2:K
        
        % accumulator
        sum_uk = u_hat_plus(n+1,:,k-1) + sum_uk - u_hat_plus(n,:,k);
        
        % mode spectrum，谱
        u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
        
        % center frequencies
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
        
    end
    
    % Dual ascent   双重上升是一种优化算法，用于解决凸优化问题中的对偶问题。
    lambda_hat(n+1,:) = lambda_hat(n,:) + tau*(sum(u_hat_plus(n+1,:,:),3) - f_hat_plus);
    
    % loop counter 循环计数器
    n = n+1;
    
    % converged yet?未收敛
    uDiff = eps; 
    for i=1:K
        uDiff = uDiff + 1/T*(u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i))*conj((u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i)))';
    end
    uDiff = abs(uDiff);
    
end
 
 
%------ Postprocessing and cleanup  后处理和清理
 
 
% discard empty space if converged early  如果提前收敛，则丢弃空空间
N = min(N,n);
omega = omega_plus(1:N,:);
 
% Signal reconstruction  信号重构
u_hat = zeros(T, K);
u_hat((T/2+1):T,:) = squeeze(u_hat_plus(N,(T/2+1):T,:));
u_hat((T/2+1):-1:2,:) = squeeze(conj(u_hat_plus(N,(T/2+1):T,:)));
u_hat(1,:) = conj(u_hat(end,:));
 
u = zeros(K,length(t));
 
for k = 1:K
    u(k,:)=real(ifft(ifftshift(u_hat(:,k))));
end
 
% remove mirror part  移除镜像部分
u = u(:,T/4+1:3*T/4);
 
% recompute spectrum，谱，重新计算频谱
clear u_hat;
for k = 1:K
    u_hat(:,k)=fftshift(fft(u(k,:)))';
end
end