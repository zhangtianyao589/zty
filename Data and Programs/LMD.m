clear all;
clc;
% 
fs=25600;
N=25600;
n=0:N-1;
t=n/fs;
L=25600;
PF=4;
load  96fz.mat
X = X(1:L)';

%% 先运行eMD,再运行这部分画图
figure(1);
imfn=PF;
subplot(PF+1,1,1); 
plot(t,X); %故障信号
ylabel('s','fontsize',12,'fontname','宋体');
 
for n1=1:PF
    subplot(PF+1,1,n1+1);
    plot(t,PF(n1,:));%输出IMF分量，a(:,n)则表示矩阵a的第n列元素，u(n1,:)表示矩阵u的n1行元素
    
    ylabel(['IMF' int2str(n1)]);%int2str(i)是将数值i四舍五入后转变成字符，y轴命名
    fontSize = 12;     
        set(gca,'FontSize', fontSize,'color','w');
        
end
 xlabel('Time [s]');
 fontSize = 12;     
        set(gca,'FontSize', fontSize,'color','w');
%% 绘制仿真信号和其频谱图
figure(1)
subplot(211)
plot(t,X)
subplot(212)
y2=X;
L=length(y2);
NFFT = 2^nextpow2(L);
Y = fft(y2,NFFT)/L;
f = fs/2*linspace(0,1,NFFT/2);
plot(f,2*abs(Y(1:NFFT/2)))
% 局域均值分析
X=X';
c = X';
N = length(X);
 
A = ones(1,N);
PF = [];
aii = 2*A;
 
while(1)
 
  si = c;
  a = 1;
  
   while(1)
    h = si;
    
      maxVec = [];
      minVec = [];
      
   % 寻找极大值点和极小值点
      for i = 2: N - 1
         if h (i - 1) < h (i) & h (i) > h (i + 1)
            maxVec = [maxVec i]; 		
         end
         if h (i - 1) > h (i) & h (i) < h (i + 1)
            minVec = [minVec i]; 		
         end         
      end
      
   % 检查是否有残余
      if (length (maxVec) + length (minVec)) < 2
         break;
      end
           
  % 原始信号中的两边两个点的判断 
      lenmax=length(maxVec);
      lenmin=length(minVec);
      %先是左边这个点
      if h(1)>0
          if(maxVec(1)<minVec(1))
              yleft_max=h(maxVec(1));
              yleft_min=-h(1);
          else
              yleft_max=h(1);
              yleft_min=h(minVec(1));
          end
      else
          if (maxVec(1)<minVec(1))
              yleft_max=h(maxVec(1));
              yleft_min=h(1);
          else
              yleft_max=-h(1);
              yleft_min=h(minVec(1));
          end
      end
      %然后判断右边这个点
      if h(N)>0
          if(maxVec(lenmax)<minVec(lenmin))
             yright_max=h(N);
             yright_min=h(minVec(lenmin));
          else
              yright_max=h(maxVec(lenmax));
              yright_min=-h(N);
          end
      else
          if(maxVec(lenmax)<minVec(lenmin))
              yright_max=-h(N);
              yright_min=h(minVec(lenmin));
          else
              yright_max=h(maxVec(lenmax));
              yright_min=h(N);
          end
      end
      %使用三次样条插值方法，对极大值向量和极小值向量进行插值
      %spline interpolate
      maxEnv=spline([1 maxVec N],[yleft_max h(maxVec) yright_max],1:N);
      minEnv=spline([1 minVec N],[yleft_min h(minVec) yright_min],1:N);
      
    mm = (maxEnv + minEnv)/2;%得到局部均值函数
    aa = abs(maxEnv - minEnv)/2;%得到包络函数
    
    mmm = mm;
    aaa = aa;
 
    preh = h;
    h = h-mmm;%从原始信号中分离处局部均值函数
    si = h./aaa;%对分离出的信号进行解调
    a = a.*aaa;    
    
aii = aaa;
 
    B = length(aii);
    C = ones(1,B);
    bb = norm(aii-C);%返回aii-C的最大奇异值，aii就是那个包络函数
    if(bb < 1000)%如果bb<1000，就得到了纯调频函数
        break;
    end     
    
   end %分解1个Pf分量在这结束
   
  pf = a.*si;%包络函数和纯调频函数相乘，得到PF分量
  
  PF = [PF; pf];
  
  bbb = length (maxVec) + length (minVec);
 % 简单的一个结束分解的条件
      if (length (maxVec) + length (minVec)) < 20
         break;
      end
           
  c = c-pf;
 
end
m=X'-PF(1,:)-PF(2,:);%如果分解出2个，从原始信号中减去，得到残余分量
figure(2);
subplot(4,1,1),plot(n,X),ylabel('X(t)');
subplot(4,1,2),plot(n,PF(1,:)),ylabel('PF_1(t)');
subplot(4,1,3),plot(n,PF(2,:)),ylabel('PF_2(t)');
subplot(4,1,4),plot(n,m),ylabel('u(t)');
xlabel('Time / S'); 