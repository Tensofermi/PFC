%预设
clc;
clear;
format long;
%tic;
%加载Hon-------------------(-0.10,0.10)
load('.\basic\Hon\1_Hon\fM.mat');
Hon_=fM;
clear fM;

%加载Dim-------------------(0.05,0.10)
load('.\basic\Dim\1_Dim\fM.mat');
Dim_=fM;
clear fM;

%加载Rec-------------------(0.00,-0.05)
load('.\basic\Rec\1_Rec\fM.mat');
Rec_=fM;
clear fM;

%加载Kag-------------------(0.10,0.00)
load('.\basic\Kag\1_Kag\fM.mat');
Kag_=fM;
clear fM;

%加载Tri-------------------(-0.40,-0.05)
load('.\basic\Tri\1_Tri\fM.mat');
Tri_=fM;
clear fM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_____mainfunction_____%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%迭代参数
Nx=128;   Ny=128;    %空间大小
dx=pi/8;  dy=pi/8;   %空间步长
nsteps=5000;   %时间步数
dtime=0.1;  %时间步长
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=dx:dx:Nx*dx;  y=dy:dy:Ny*dy;%设置空间标度矩阵
nprint=500;   %输出图像步长间隔
nstart=1;      %起始迭代步长

%自由能参数
den0=-0.2; %初始密度——重要可调参量
r=-0.15; %温度相关量——重要可调参量
lamda=0.02;  %可调参数，通常为定值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
km1=1; km2=sqrt(3); km3=2;  %矢量模
b1=0;  b2=0.05; b3=.4;     %权重
co=Tri_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0;     %立方项系数
V0=0;    %薄膜-基底势能
size=50; %放置晶体大小,注意要小于128,且必须是偶数




%添加势能(三角晶格)
for i=1:Nx
    for j=1:Ny
     %   V(i,j)=2*V0;
V(i,j)=2*V0*(cos(-(sqrt(3)/2)*i-(1/2)*j)+cos(j)+cos((sqrt(3)/2)*i-(1/2)*j));
    end
end

%添加势能(Tri)
V=V0.*Tri_;

%添加势能(Hon)
%V=V0.*Hon_;

%添加势能(Dim)
%V=V0.*Dim_;

%添加势能(Rec)
%V=V0.*Rec_;

%添加势能(Kag)
%V=V0.*Kag_;


%初始化倒空间
[kx,ky,k2,k4] = prepare_fft2(Nx,Ny,dx,dy);



px=128/2;   py=128/2;%放置中心
p1x=Nx-128;   p1y=128;%放置中心
p2x=128;   p2y=Nx-128;%放置中心
p3x=Nx-128;   p3y=Nx-128;%放置中心
ppp1=rot90(Hon_); %转动90度
ppp2=rot90(rot90(Hon_));%转动180度
ppp3=rot90(rot90(rot90(Hon_)));%转动270度
%初始化原子密度
    for i=1:Nx
      for j=1:Ny
          if((i>(px-size/2))&&(i<(px+size/2))&&(j<(py+size/2))&&(j>(py-size/2)))   den(i,j)= co(i-(px-size/2),j-(py-size/2));
     %     elseif ((i>(p1x-size/2))&&(i<(p1x+size/2))&&(j<(p1y+size/2))&&(j>(p1y-size/2)))   den(i,j)=ppp1(i-(p1x-size/2),j-(p1y-size/2));
     %         elseif ((i>(p2x-size/2))&&(i<(p2x+size/2))&&(j<(p2y+size/2))&&(j>(p2y-size/2)))   den(i,j)=ppp2(i-(p2x-size/2),j-(p2y-size/2));
     %                  elseif ((i>(p3x-size/2))&&(i<(p3x+size/2))&&(j<(p3y+size/2))&&(j>(p3y-size/2)))   den(i,j)=ppp3(i-(p3x-size/2),j-(p3y-size/2));
          else den(i,j)=den0+den0*0.01*(0.5-rand);  
            end
      %  den(i,j)=den0+den0*0.01*(0.5-rand);
      end
    end

%保存初态矩阵    
iM=den;
save('iM.mat','iM');

%den=fM;   %加载之前计算的结果


%生成倒空间直相关函数
  for i=1:Nx
   for j=1:Ny
    C(i,j)=lamda*(k4(i,j)-2*km1^2*k2(i,j)+km1^4+b1).*(k4(i,j)-2*km2^2*k2(i,j)+km2^4+b2).*(k4(i,j)-2*km3^2*k2(i,j)+km3^4+b3);   %根据C(2)的具体形式赋值
   end
  end
  
%在倒空间中对原子密度进行演化
    for istep=nstart:nsteps
        f_den=fftshift(fft2(den));              %对den二维傅里叶变换
        
        den2=den.^2;
        f_den2=fftshift(fft2(den2));              %对den^2二维傅里叶变换
        
        den3=den.^3;
        f_den3=fftshift(fft2(den3));              %对den^3二维傅里叶变换
        
        f_V=fftshift(fft2(V));              %对V二维傅里叶变换
        
      f_den=(f_den+dtime*k2.*(f_den2*a-f_den3-f_V))./(1+dtime.*k2.*(r+C));       %迭代
        
        den=real(ifft2(ifftshift(f_den)));         %逆傅里叶变换回实空间，注意取real值
        
        energy=get_energy(den,r,a,Nx,Ny,k2,k4,lamda,km1,km2,km3,b1,b2,b3,V);  %获得当前体系的自由能
        energylist(istep)=energy;                  %并加入energylist当中
        
        
  if((mod(istep,nprint)==0) ||(istep==1) ) %默认nprint步更新一次图像
%输出密度图像
            figure(1);  
            pcolor(x,y,den); 
             axis square;
           shading interp;
           colorbar;
%保存密度图像
      saveas(1,['istep',num2str(istep),'dt=',num2str(dtime),'.png']);
     
%绘制自由能演化图
             figure(2);
            plot(energylist);   
            end    
        pause(eps);%输出时间间隔为matlab最小值eps
    end
    
%保存终态矩阵以及能量曲线    
fM=den;
save('fM.mat','fM');

out1 = fopen('den(i,j).out','w');%%

fprintf(out1,'system parameters (dx, dy, dt, Nx, Ny, T) = ');
fprintf(out1,'(%d, %d, %d, %d, %d, %d) \n',dx,dy,dtime,Nx,Ny,nsteps);
fprintf(out1,'free function parameters (den0, r, lamda, {km1, km2, km3}, {b1, b2, b3}, a, V0) = ');
fprintf(out1,'(%d, %d, %d, {%d, %d, %d}, {%d, %d, %d}, %d, %d) \n',den0,r,lamda,km1,km2,km3,b1,b2,b3,a,V0);
fprintf(out1,'den(i,j):\n');
for i=1:Nx
for j=1:Ny
fprintf(out1,'%5d    ',den(i,j));
end
fprintf(out1,'\n ');
end


out2 = fopen('energylist.out','w');

for i=1:nsteps
fprintf(out2,' %d    %d\n',energylist(i),i);
end

fclose('all');

%输出计算所花时间
%toc;

min(energylist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_____subfunction_____%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare_fft2
%将原点定在矩阵的中心
function [kx,ky,k2,k4]=prepare_fft2(Nx,Ny,dx,dy)
for i=1:Nx+1                                   %构建k空间的kx轴
kx(i)=(i-(Nx/2+1))*(2*pi)/(Nx*dx);
end

for i=1:Ny+1                                   %构建k空间的ky轴
ky(i)=(i-(Ny/2+1))*(2*pi)/(Ny*dy);
end

for i=1:Nx                                         %构建kx+ky的模方，即对应k^2和nabla^2
for j=1:Ny
k2(i,j)=kx(i)^2+ky(j)^2;
end
end

k4=k2.^2;                                         %构建kx+ky的模方的平方，即对应k^4和nabla^4

k6=k2.^3;                                         %构建kx+ky的模方的平方，即对应k^6和nabla^6
end

%get_energy
%不同自由能表达式能量不同，需要及时修改！！
function energy=get_energy(den,r,a,Nx,Ny,k2,k4,lamda,km1,km2,km3,b1,b2,b3,V)

ss2=den.^2;  ss3=den.^3;  ss4=ss2.^2;     %获得密度不同幂次

f_den=fftshift(fft2(den));     %k空间的密度场
modes=(k4-2*km1^2*k2+km1^4+b1).*(k4-2*km2^2*k2+km2^4+b2).*(k4-2*km3^2*k2+km3^4+b3);  %三模
f_ff=0.5*lamda*f_den.*modes; %k空间计算拉普拉斯算符项
ff=real(ifft2(ifftshift(f_ff))).*den+0.5*r*ss2-(a/3)*ss3+0.25*ss4+V.*den; %根据自由能表达式输出能量矩阵

dx=pi/8;  dy=pi/8;
energy=sum(ff(:))/(Nx*Ny*dx*dy);      %对每个能量矩阵元素求和

end
