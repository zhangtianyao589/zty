%% 周期性脉冲
close all
clear
clc
L = 2560;
fs = 25600;                  % 采样频率
% fn = 3000;                   % 固有频率
y0 = 5;                      % 位移常数
g = 0.1;                     % 阻尼系数
T = 0.01;                   % 重复周期
N = 25600;                  % 采样点数
NT = round(fs*T);      % 单周期采样点数
t = 0:1/fs:(N-1)/fs;      % 采样时刻
t0 = 0:1/fs:(NT-1)/fs;  % 单周期采样时刻
K = 100;       % 重复次数
x1 = [];
for i = 1:K
%     y = [y,exp(-g*2*pi*fn*t0).*sin(2*pi*fn*sqrt(1-g^2)*t0)];
    x1 = [x1,exp(-1000*(t0)).*sin(2*pi*1700*(t0-i/100-0.01))];
end
x1 = x1(1:N);
Yf = fft(x1);                % 频谱
figure(1)
plot(t,x1);
axis([0,inf,-4,4])
title('')
xlabel('Time [s]');
ylabel('Amplitude')
xlim([0, 1]); 
    xticks(0:0.2:1); % 设置X轴刻度的位置
    xticklabels(0:0.2:1); % 设置X轴刻度的标签
ylim([-1, 1]);
   yticks(-1:1:1);
   fontSize = 12;     
        set(gca,'FontSize', fontSize,'color','w');

        %%  随机脉冲

L2 = 25600;
fs2 = 25600;                  % 采样频率
y02 = 5;                      % 位移常数
g2 = 0.1;                     % 阻尼系数
T2_min = 0.1;                 % 最小重复周期
T2_max = 0.3;                 % 最大重复周期
N2 = 25600;                   % 采样点数
t2 = 0:1/fs2:(N2-1)/fs2;      % 采样时刻
K2 = 7;                       % 重复次数
x2 = [];

for i = 1:K2
    random_amplitude = 0.5 + (1 - 0.5) * rand();     % 生成0.5到1之间的随机数作为幅值
    random_period = (T2_max - T2_min) * rand() + T2_min;   % 生成随机的周期值
    random_start = rand() * (1 - random_period);  % 生成0到(1-周期)之间的随机数作为起始时刻
    t02 = 0:1/fs2:random_period - 1/fs2;     % 重新生成单周期采样时刻
    x2 = [x2, random_amplitude*exp(-800*t02).*sin(2*pi*2000*t02)];
    
    % 添加时间间隔
    if i < K2
        interval_time = rand() * (1 - random_period);  % 生成0到(1-周期)之间的随机数作为时间间隔
        interval_samples = round(interval_time * fs2);   % 时间间隔对应的采样点数
        x2 = [x2, zeros(1, interval_samples)];
    end
end


x2 = x2(1:L2);  % 截取指定长度的信号



Yf = fft(x2);                % 频谱
figure(2)
plot(t2,x2);
axis([0,inf,-4,4])
title('')
xlim([0, 1]); 
    xticks(0:0.2:1); % 设置X轴刻度的位置
    xticklabels(0:0.2:1); % 设置X轴刻度的标签
ylim([-1, 1]);
   yticks(-1:1:1);
ylabel('Amplitude')
xlabel('Time [s]');
   fontSize = 12;     
   set(gca,'FontSize', fontSize,'color','w');

%%  离散谐波干扰


L3 = 25600;
fs3 = 25600;                  % 采样频率
N3 = 25600;                  % 采样点数
t3 = 0:1/fs3:(N3-1)/fs3;
x3 = 0.025*sin(2*pi*320*t3+pi/6)+0.025*sin(2*pi*530*t3-pi/3)+0.05*((0.5+0.3*sin(2*pi*7.5*t3)).*cos(2*pi*260*t3+1.2*sin(2*pi*14*t3)));
x3 = x3(1:25600);
% 绘制图形
figure(3)
plot(t3,x3);
title('')
ylim([-0.1, 0.1]);
   yticks(-0.1:0.1:0.1);
ylabel('Amplitude')
xlabel('Time [s]');
   fontSize = 12;     
   set(gca,'FontSize', fontSize,'color','w');

%%
% 参数设置
fs4 = 25600; % 采样率为1000Hz
n4 = fs4; % 信号长度为一秒钟的采样数
mean_val4 = 0; % 均值
std_dev4 = 0.5; % 标准差

% 生成高斯白噪声
x4 = mean_val4 + std_dev4 * randn(1, n4);

% 时间向量
t4 = (0:n4-1) / fs4; % 生成时间向量，单位为秒

% 绘制图形
figure(4)
plot(t4, x4)
title('')
ylabel('Amplitude')
xlabel('Time [s]');
ylim([-5, 5]);
   yticks(-5:5:5);
   fontSize = 12;     
   set(gca,'FontSize', fontSize,'color','w');

%%
X =x1+x2+ x3+x4;