x0=0;y0=100;%光源
Ax=0;Ay=10;
Bx=50;By=20;
km=(By-Ay)/(Bx-Ax);
jm=Ay-km*Ax;
z=0:0.0001:50;
M=km*z+jm;

Cx=90;Cy=60;
Dx=30;Dy=100;
kn=(Cy-Dy)/(Cx-Dx);
jn=Cy-kn*Cx;
N=kn*z+jn;

Nn=5;
Px=z(ranperm(numel(z),Nn));
Py=km*Px+jm;
plot([Ax,Bx],[Ay,By],'k','Linewidth',2);
plot([Cx,Dx],[Cy,Dy],'r','Linewidth',2)

for i=1:Nn
    plot([x0,Px(i)],[y0,Py(i)],'Linewidth',2);
end
hold on;

daspect([1,1,1]);

alpha=zeros(1,Nn);
inan=zeros(1,Nn);
for i=1:Nn
    L=norm([])
    alpha(i)=acos((m^2+n^2-i^2)/(2*m*n));
    inan(i)=(pi/2-alpha(i))/pi*180;
end
for i=1:Nn
    kr(i)=tan(alpha+atan(km));
    jr(i)=
    if kr(i)>=k_AC & kr(i)<=k-AD
        [X(i),Y(i)]=linecross()
        plot(Y ...)
    else
        Xx=3*Px(i);
        Yy=kd*xx
    end
end





