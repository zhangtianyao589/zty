
clc
clear all
fs=25600;%采样频率
Ts=1/fs;%采样周期
L=25600;%采样点数
t=(0:L-1)*Ts;%时间序列
STA=1; %采样起始位置


%----------------导入数据-----------------
load 100fz.mat
X = X(1:L)';
% A=xlsread('20151124_08_15Bin粗卡阀.xlsx');
% x=A(1:3072,1)';
%ITD分解，分解成6个pr分量
iterated_max=10;
[H,L]=Itd4(X,6); %ITD分解
%% ------------分量信号显示-------------
figure(1)
n=size(H,1);
subplot(n+1,1,1);
plot(t,X);%原始信号
ylabel('原始信号');
for n1=1:n
    subplot(n+1,1,n1+1);
    plot(t,H(n1,:));%分量信号
  
    ylabel(['IMF' int2str(n1)]);
    fontSize = 12;     
        set(gca,'FontSize', fontSize,'color','w');
end
xlabel('Time [s]');
fontSize = 12;     
        set(gca,'FontSize', fontSize,'color','w');


 %% 求各imf分量与原始信号的相关系数，并筛选分量
xiangguanxishu=zeros(1,n);
for i=1:n
    a=H(i,:);
    b=corrcoef(a,X');
    xiangguanxishu(1,i)=b(1,2);
end

%提取每一个IMF分量的信息熵特征
xingxishang=zeros(1,n);
for i=1:n
    a=H(i,:);
    x1=mapminmax(a, 0, 1); %将数据归一化处理。
    %求数据的概率分布
   for m=1:3072
       x=x1(m)/sum(x1);
   end
   shannon=-sum(x.*log2(x));%计算信息熵
   xingxishang(1,i)=shannon;%把每一个分量的信息熵构造成一个特征向量
end
