%输出衍射图样
figure(1)
S=(f_den).*conj(f_den);
pcolor(x*(2*pi)/(Nx*dx),y*(2*pi)/(Nx*dx),S);
colormap(flipud(gray));
axis square;
shading interp;
saveas(1,['diffraction_pattern','.png']);


%衍射中心坐标(通常为(Nx/2+1,Ny/2+1))
x0=Nx/2+1;
y0=Ny/2+1;

%离散化的距离矩阵
dist2 = [];
index=0;
for i=-200:200
    for j=-200:200
    index=index+1;
    dist2(index)=i*i+j*j;
    end
end
dist2=unique(sort(dist2));
dist2(1)=[];

%开始遍历
pp=[];
for m=1:400
   %初始化_计数count和衍射权重p
count=0;
p=0;
for i=1:Nx
    for j=1:Ny
        %中心点到遍历点的距离平方满足离散矩阵的距离，即开始操作
   if((i-x0)^2+(j-y0)^2==dist2(m))
       p=p+S(i,j);
       count=count+1;
   end
    end
end
pp(m)=p/count;
end
%绘制衍射因子
figure(2)

for i=1:400
q_index(i)=(2*pi)/(Nx*dx)*sqrt(dist2(i));
end
plot(q_index,pp)


out3 = fopen('structure_factor.out','w');

for i=1:400
fprintf(out3,' %d    %d,\n',sqrt(dist2(i)),pp(i));
end
fclose('all');