%Ԥ��
clc;
clear;
format long;
%tic;
%����Hon-------------------(-0.10,0.10)
load('.\basic\Hon\1_Hon\fM.mat');
Hon_=fM;
clear fM;

%����Dim-------------------(0.05,0.10)
load('.\basic\Dim\1_Dim\fM.mat');
Dim_=fM;
clear fM;

%����Rec-------------------(0.00,-0.05)
load('.\basic\Rec\1_Rec\fM.mat');
Rec_=fM;
clear fM;

%����Kag-------------------(0.10,0.00)
load('.\basic\Kag\1_Kag\fM.mat');
Kag_=fM;
clear fM;

%����Tri-------------------(-0.40,-0.05)
load('.\basic\Tri\1_Tri\fM.mat');
Tri_=fM;
clear fM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_____mainfunction_____%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%��������
Nx=128;   Ny=128;    %�ռ��С
dx=pi/8;  dy=pi/8;   %�ռ䲽��
nsteps=5000;   %ʱ�䲽��
dtime=0.1;  %ʱ�䲽��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=dx:dx:Nx*dx;  y=dy:dy:Ny*dy;%���ÿռ��Ⱦ���
nprint=500;   %���ͼ�񲽳����
nstart=1;      %��ʼ��������

%�����ܲ���
den0=-0.2; %��ʼ�ܶȡ�����Ҫ�ɵ�����
r=-0.15; %�¶������������Ҫ�ɵ�����
lamda=0.02;  %�ɵ�������ͨ��Ϊ��ֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
km1=1; km2=sqrt(3); km3=2;  %ʸ��ģ
b1=0;  b2=0.05; b3=.4;     %Ȩ��
co=Tri_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=0;     %������ϵ��
V0=0;    %��Ĥ-��������
size=50; %���þ����С,ע��ҪС��128,�ұ�����ż��




%�������(���Ǿ���)
for i=1:Nx
    for j=1:Ny
     %   V(i,j)=2*V0;
V(i,j)=2*V0*(cos(-(sqrt(3)/2)*i-(1/2)*j)+cos(j)+cos((sqrt(3)/2)*i-(1/2)*j));
    end
end

%�������(Tri)
V=V0.*Tri_;

%�������(Hon)
%V=V0.*Hon_;

%�������(Dim)
%V=V0.*Dim_;

%�������(Rec)
%V=V0.*Rec_;

%�������(Kag)
%V=V0.*Kag_;


%��ʼ�����ռ�
[kx,ky,k2,k4] = prepare_fft2(Nx,Ny,dx,dy);



px=128/2;   py=128/2;%��������
p1x=Nx-128;   p1y=128;%��������
p2x=128;   p2y=Nx-128;%��������
p3x=Nx-128;   p3y=Nx-128;%��������
ppp1=rot90(Hon_); %ת��90��
ppp2=rot90(rot90(Hon_));%ת��180��
ppp3=rot90(rot90(rot90(Hon_)));%ת��270��
%��ʼ��ԭ���ܶ�
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

%�����̬����    
iM=den;
save('iM.mat','iM');

%den=fM;   %����֮ǰ����Ľ��


%���ɵ��ռ�ֱ��غ���
  for i=1:Nx
   for j=1:Ny
    C(i,j)=lamda*(k4(i,j)-2*km1^2*k2(i,j)+km1^4+b1).*(k4(i,j)-2*km2^2*k2(i,j)+km2^4+b2).*(k4(i,j)-2*km3^2*k2(i,j)+km3^4+b3);   %����C(2)�ľ�����ʽ��ֵ
   end
  end
  
%�ڵ��ռ��ж�ԭ���ܶȽ����ݻ�
    for istep=nstart:nsteps
        f_den=fftshift(fft2(den));              %��den��ά����Ҷ�任
        
        den2=den.^2;
        f_den2=fftshift(fft2(den2));              %��den^2��ά����Ҷ�任
        
        den3=den.^3;
        f_den3=fftshift(fft2(den3));              %��den^3��ά����Ҷ�任
        
        f_V=fftshift(fft2(V));              %��V��ά����Ҷ�任
        
      f_den=(f_den+dtime*k2.*(f_den2*a-f_den3-f_V))./(1+dtime.*k2.*(r+C));       %����
        
        den=real(ifft2(ifftshift(f_den)));         %�渵��Ҷ�任��ʵ�ռ䣬ע��ȡrealֵ
        
        energy=get_energy(den,r,a,Nx,Ny,k2,k4,lamda,km1,km2,km3,b1,b2,b3,V);  %��õ�ǰ��ϵ��������
        energylist(istep)=energy;                  %������energylist����
        
        
  if((mod(istep,nprint)==0) ||(istep==1) ) %Ĭ��nprint������һ��ͼ��
%����ܶ�ͼ��
            figure(1);  
            pcolor(x,y,den); 
             axis square;
           shading interp;
           colorbar;
%�����ܶ�ͼ��
      saveas(1,['istep',num2str(istep),'dt=',num2str(dtime),'.png']);
     
%�����������ݻ�ͼ
             figure(2);
            plot(energylist);   
            end    
        pause(eps);%���ʱ����Ϊmatlab��Сֵeps
    end
    
%������̬�����Լ���������    
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

%�����������ʱ��
%toc;

min(energylist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_____subfunction_____%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare_fft2
%��ԭ�㶨�ھ��������
function [kx,ky,k2,k4]=prepare_fft2(Nx,Ny,dx,dy)
for i=1:Nx+1                                   %����k�ռ��kx��
kx(i)=(i-(Nx/2+1))*(2*pi)/(Nx*dx);
end

for i=1:Ny+1                                   %����k�ռ��ky��
ky(i)=(i-(Ny/2+1))*(2*pi)/(Ny*dy);
end

for i=1:Nx                                         %����kx+ky��ģ��������Ӧk^2��nabla^2
for j=1:Ny
k2(i,j)=kx(i)^2+ky(j)^2;
end
end

k4=k2.^2;                                         %����kx+ky��ģ����ƽ��������Ӧk^4��nabla^4

k6=k2.^3;                                         %����kx+ky��ģ����ƽ��������Ӧk^6��nabla^6
end

%get_energy
%��ͬ�����ܱ��ʽ������ͬ����Ҫ��ʱ�޸ģ���
function energy=get_energy(den,r,a,Nx,Ny,k2,k4,lamda,km1,km2,km3,b1,b2,b3,V)

ss2=den.^2;  ss3=den.^3;  ss4=ss2.^2;     %����ܶȲ�ͬ�ݴ�

f_den=fftshift(fft2(den));     %k�ռ���ܶȳ�
modes=(k4-2*km1^2*k2+km1^4+b1).*(k4-2*km2^2*k2+km2^4+b2).*(k4-2*km3^2*k2+km3^4+b3);  %��ģ
f_ff=0.5*lamda*f_den.*modes; %k�ռ����������˹�����
ff=real(ifft2(ifftshift(f_ff))).*den+0.5*r*ss2-(a/3)*ss3+0.25*ss4+V.*den; %���������ܱ��ʽ�����������

dx=pi/8;  dy=pi/8;
energy=sum(ff(:))/(Nx*Ny*dx*dy);      %��ÿ����������Ԫ�����

end
