VMAX=max(max(V));
VMIN=min(min(V));
DENMAX=max(max(den));
DENMIN=min(min(den));

den=den+abs(DENMIN);
V=V+abs(VMIN);


t1=[VMIN,VMAX,VMAX-VMIN];
t2=[DENMIN,DENMAX,DENMAX-DENMIN];

t=sum(sum(abs(den)))/sum(sum(abs(V)));


temp=[];

for i=1:Nx
    for j=1:Ny
    temp(i,j)=0;
    end
end


for k=1:Ny
for i=1:Nx
    for j=1:k
    temp(i,j)=den(i,Ny-k+j);
    end
end

hhh=temp+t.*V;

for i=1:Nx
    for j=1:Ny
 if(hhh(i,j)>max(max(den)))hhh(i,j)=1;
 end
    end
end
     pcolor(x,y,hhh); 
             axis square;
           shading interp;
           colorbar;
            saveas(1,['moire_k=',num2str(k),'.png']);
           pause(eps);

end



