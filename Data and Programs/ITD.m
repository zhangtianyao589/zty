
clc
clear all
fs=25600;%����Ƶ��
Ts=1/fs;%��������
L=25600;%��������
t=(0:L-1)*Ts;%ʱ������
STA=1; %������ʼλ��


%----------------��������-----------------
load 100fz.mat
X = X(1:L)';
% A=xlsread('20151124_08_15Bin�ֿ���.xlsx');
% x=A(1:3072,1)';
%ITD�ֽ⣬�ֽ��6��pr����
iterated_max=10;
[H,L]=Itd4(X,6); %ITD�ֽ�
%% ------------�����ź���ʾ-------------
figure(1)
n=size(H,1);
subplot(n+1,1,1);
plot(t,X);%ԭʼ�ź�
ylabel('ԭʼ�ź�');
for n1=1:n
    subplot(n+1,1,n1+1);
    plot(t,H(n1,:));%�����ź�
  
    ylabel(['IMF' int2str(n1)]);
    fontSize = 12;     
        set(gca,'FontSize', fontSize,'color','w');
end
xlabel('Time [s]');
fontSize = 12;     
        set(gca,'FontSize', fontSize,'color','w');


 %% ���imf������ԭʼ�źŵ����ϵ������ɸѡ����
xiangguanxishu=zeros(1,n);
for i=1:n
    a=H(i,:);
    b=corrcoef(a,X');
    xiangguanxishu(1,i)=b(1,2);
end

%��ȡÿһ��IMF��������Ϣ������
xingxishang=zeros(1,n);
for i=1:n
    a=H(i,:);
    x1=mapminmax(a, 0, 1); %�����ݹ�һ������
    %�����ݵĸ��ʷֲ�
   for m=1:3072
       x=x1(m)/sum(x1);
   end
   shannon=-sum(x.*log2(x));%������Ϣ��
   xingxishang(1,i)=shannon;%��ÿһ����������Ϣ�ع����һ����������
end
