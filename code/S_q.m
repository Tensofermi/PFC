%�������ͼ��
figure(1)
S=(f_den).*conj(f_den);
pcolor(x*(2*pi)/(Nx*dx),y*(2*pi)/(Nx*dx),S);
colormap(flipud(gray));
axis square;
shading interp;
saveas(1,['diffraction_pattern','.png']);


%������������(ͨ��Ϊ(Nx/2+1,Ny/2+1))
x0=Nx/2+1;
y0=Ny/2+1;

%��ɢ���ľ������
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

%��ʼ����
pp=[];
for m=1:400
   %��ʼ��_����count������Ȩ��p
count=0;
p=0;
for i=1:Nx
    for j=1:Ny
        %���ĵ㵽������ľ���ƽ��������ɢ����ľ��룬����ʼ����
   if((i-x0)^2+(j-y0)^2==dist2(m))
       p=p+S(i,j);
       count=count+1;
   end
    end
end
pp(m)=p/count;
end
%������������
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