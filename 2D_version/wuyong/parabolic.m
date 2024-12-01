x0=4;y0=10;
focus=[x0,y0];
k=2;
b=5;
X=-4:0.1:0;


% f0=@(x) k.*x+b;
% f1=@(x) -1/k.*x+y0+x0/k;
hold on;
grid on;
daspect([1 1 1]);

% fplot(f0,'LineWidth',2);%准线
% plot(X-focus(1),M-focus(2),'b-','LineWidth',1);%准线why是加
% plot(X,M,'r--','LineWidth',1);%准线why是加
% fplot(f1,'--','LineWidth',1);%中轴线
% scatter(x0,y0,'*r','LineWidth',1);%光源

%red 平移变换
M=k.*(X)+b+focus(2);
plot(X+focus(1),M,'g.','LineWidth',1.5); %在同一条直线上但是点取得不好

N=k.*X+b;
% plot(X,N,'g-','LineWidth',1);%准线why是加
plot(X+focus(1),N+focus(2),'y-','LineWidth',1);%让已经画好的每一个点像焦点一样移动

f5=@(x) k.*(x-focus(1))+b+focus(2);
fplot(f5,'r--','LineWidth',1);
scatter(x0,y0,'r*','LineWidth',1);%光源
f2=@(x,y) (x-focus(1)+k.*(y-focus(2))).^2-2.*k.*b.*(x-focus(1))+2.*b.*(y-focus(2))-b^2;
fimplicit(f2,'r','LineWidth',2);%抛物面镜
f1=@(x) -1/k.*x+y0+x0/k;
fplot(f1,'r--','LineWidth',1);%中轴线

%black 原
% N=k.*X+b;
% plot(X,N,'g.','LineWidth',1);
f0=@(x) k.*x+b;
fplot(f0,'k--','LineWidth',1);
scatter(0,0,'k*','LineWidth',1);
f3=@(x,y) (x+k.*(y)).^2-2.*k.*b.*(x)+2.*b.*(y)-b^2;
fimplicit(f3,'k','LineWidth',2);
f4=@(x)-1/k.*x;
fplot(f4,'k--','LineWidth',1);%中轴线


hold off;