%%
x0=0;y0=100;
Ax=25;Ay=10;
Bx=60;By=12;
Km=(By-Ay)/(Bx-Ax); %mirror斜率
Jm=Ay-Km*Ax;
Z=25:1:60;
M=Km*Z+Jm;

Cx=80;Cy=70;
Dx=50;Dy=80;
Kn=(Dy-Cy)/(Dx-Cx); %cd斜率
Jn=Cy-Kn*Cx;
N=Kn*Z+Jn;

nn=5;
Px=Z(randperm(numel(Z),nn));
Py=Km.*Px+Jm;
for i=1:nn
plot([x0,Px(i)],[y0,Py(i)],'Linewidth',1.5);
hold on;
end
daspect([1 1 1]);
plot([Ax,Bx],[Ay,By],'k','Linewidth',2);
plot([Cx,Dx],[Cy,Dy],'r','Linewidth',2);
alpha=zeros(1,nn);
inan=zeros(1,nn);
for i=1:nn
l=norm([x0-Ax,y0-Ay]);
m=norm([x0-Px(i),y0-Py(i)]);
n=norm([Ax-Px(i),Ay-Py(i)]);
alpha(i)=acos((m^2+n^2-l^2)/(2*m*n));
inan(i)=(pi/2-alpha(i))/pi*180;
end
inan;

Kr=zeros(1,nn);
Jr=zeros(1,nn);
x=zeros(1,nn);
y=zeros(1,nn);
for i=1:nn
Kr(i)=tan(alpha(i)+atan(Km));
Jr(i)=Py(i)-Kr(i)*Px(i);
[x(i),y(i)]=linecross(Kr(i),Jr(i),Kn,Jn);
if x(i)>=Dx && x(i)<=Cx
plot([Px(i),x(i)],[Py(i),y(i)],'Linewidth',1.5);
else
xx=x(i)*1.05;
yy=Kr(i)*xx+Jr(i);
plot([Px(i),xx],[Py(i),yy],'Linewidth',1.5);
end
end
%%
function [x,y]=linecross(k1,b1,k2,b2)
x=[];
y=[];
if k1==k2&&b1==b2
disp('重合');
elseif k1==k2&&b1~=b2
disp('无交点');
else
x=(b2-b1)/(k1-k2);
y=k1*x+b1;
end
end 