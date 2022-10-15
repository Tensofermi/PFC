VMAX=max(max(V));
VMIN=min(min(V));
DENMAX=max(max(den));
DENMIN=min(min(den));

den=den+abs(DENMIN);
V=V+abs(VMIN);


t1=[VMIN,VMAX,VMAX-VMIN];
t2=[DENMIN,DENMAX,DENMAX-DENMIN];

t=sum(sum(abs(den)))/sum(sum(abs(V)));
%t=max(max(den))/max(max(V));

figure(1)
    pcolor(x,y,den); 
             axis square;
           shading interp;
           colorbar;
saveas(1,['film_pattern','.png']);           
figure(2)
    pcolor(x,y,V); 
             axis square;
           shading interp;
           colorbar;
saveas(2,['substrate_pattern','.png']);
figure(3)
    pcolor(x,y,den+t.*V); 
   
             axis square;
           shading interp;
           colorbar;
    %±£¥Ê√‹∂»ÕºœÒ
saveas(3,['moire_pattern','.png']);

fclose('all');
            