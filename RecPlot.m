function RecPlot(x0,y0,Lx,Ly)
x=[x0-Lx/2, x0+Lx/2, x0+Lx/2, x0-Lx/2];
y=[y0-Ly/2, y0-Ly/2, y0+Ly/2, y0+Ly/2];
index=zeros(4,2);
index(1,:)=[1 2];
index(2,:)=[2 3];
index(3,:)=[3 4];
index(4,:)=[4 1];
for k=1:size(index,1)
    plot(x(index(k,:)),y(index(k,:)),'k')
    hold on
end