format compact; clc
[A,B]=cursor_info.Position
dis=((A(1)-B(1))^2+(A(2)-B(2))^2)^0.5
kAB=((A(2)-B(2))/(A(1)-B(1)))
kccd=-1/kAB;
xf=60;yf=0;%焦点参数
%ccd
delM=5;delm=0.8;
hold on;
plot([xf-delM,xf-delm],[yf-delM*kccd,yf-delm*kccd],'r','LineWidth',2);
plot([xf+delM,xf+delm],[yf+delM*kccd,yf+delm*kccd],'r','LineWidth',2);
axis equal
hold off;