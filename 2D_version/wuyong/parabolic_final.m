clear;
%比例尺
daspect([1 1 1]);
hold on;
%焦点
x0=0;y0=0;
focus0=[x0,y0];
scatter(x0,y0,'k*');
%准线斜率截距(普通方程需要的量)
k0=0.5;b0=100;
%换算为 参数方程 需要的量
p0=b0/sqrt(k0^2+1); %(0,0)代入距离公式
phi0=atand(k0)-90; %旋转角度 
t0=10/180*pi:0.005:15/180*pi;
%
xp0=2*p0.*t0.^2-p0/2;
yp0=2*p0.*t0;
%对每个点,绕原点旋转变换
for i=1:numel(xp0)
    temp=num2cell([cosd(phi0),-sind(phi0);sind(phi0),cosd(phi0)]*[xp0(i);yp0(i)]);%旋转矩阵
    [xp0(i),yp0(i)]=deal(temp{:});%好像必须这样赋值(先变为cell数组，然后deal),有好方法告诉我
end
plot(xp0+x0,yp0+y0,'r-','LineWidth',2);%抛物线
Nn=3;%光线数
Pt=t0(randperm(numel(t0),Nn));
Px=2.*p0.*Pt.^2-p0/2+x0;
Py=2.*p0.*Pt+y0;
for i=1:numel(Px)
    temp=num2cell([cosd(phi0),-sind(phi0);sind(phi0),cosd(phi0)]*[Px(i);Py(i)]);%旋转矩阵
    [Px(i),Py(i)]=deal(temp{:});%好像必须这样赋值(先变为cell数组，然后deal),有好方法告诉我
end
scatter(Px,Py,'k*');
%平面镜1
Ax=24;Ay=20; %平面镜左端
Bx=43;By=20; %平面镜右端
km1=(By-Ay)/(Bx-Ax); %镜斜率
jm1=Ay-km1*Ax; %镜截距
plot([Ax,Bx],[Ay,By],'k','Linewidth',2);
%抛物面镜2
%焦点
x2=60;y2=0;scatter(x2,y2,'k*')
k2=-1/2; b2=100;
p2=b2/sqrt(k2^2+1); %(0,0)代入距离公式
phi2=atand(k2)-90; %旋转角度
t2=-8/180*pi:-0.005:-14/180*pi;
xp2=2*p2.*t2.^2-p2/2;
yp2=2*p2.*t2;
for i=1:numel(xp2)
    temp=num2cell([cosd(phi2),-sind(phi2);sind(phi2),cosd(phi2)]*[xp2(i);yp2(i)]);%旋转矩阵
    [xp2(i),yp2(i)]=deal(temp{:});%好像必须这样赋值(先变为cell数组，然后deal),有好方法告诉我
end
plot(xp2+x2,yp2+y2,'r-','LineWidth',2);%抛物线
KJ0=[];
KJ1=[];
KJ2=[];
syms x y
for i=1:Nn
    coArray=['#0072BD';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#4DBEEE';'#A2142F'];%颜色库
    color=coArray(mod(i,7)+1,:); %根据i选定一种颜色
    kr0=pslope(Px(i),Py(i),x0,y0,k0,b0);
    [xr1,yr1,flag]=reflect(x0,y0,Px(i),Py(i),kr0,km1,jm1,Ax,Bx,i,1);
    if flag==0
        continue;
    end
    [xr2,yr2,flag]=reflect(Px(i),Py(i),xr1,yr1,0,0,40,10,20,i,0);
    %求交点
    s=solve(((x-x2)+k2.*(y-y2)).^2-2.*k2.*b2.*(x-x2)+2.*b2.*(y-y2)-b2.^2==0, ...
        (xr2-xr1).*(y-yr2)==(yr2-yr1).*(x-xr2),x,y);
    x_c0=double(s.x); 
    y_c0=double(s.y);
    plot([xr1,x_c0],[yr1,y_c0],'Color',color,'Linewidth',1.5);
    kr02=pslope(x_c0,y_c0,x2,y2,k2,b2);
    [xr03,yr03,flag]=reflect(xr1,yr1,x_c0,y_c0,kr02,0,0,55,65,i,1);
    kr03=(yr03-y_c0)/(xr03-x_c0);jr03=y_c0-kr03*x_c0;
    KJ0(1,i)=kr03;KJ0(2,i)=jr03;
    
    %grating初试1
    s=solve(((x-x2)+k2.*(y-y2)).^2-2.*k2.*b2.*(x-x2)+2.*b2.*(y-y2)-b2.^2==0, ...
        (xr2-xr1).*(y-yr2)==2*(yr2-yr1).*(x-xr2),x,y);
    x_c1=double(s.x);  x_c1= x_c1(2);
    y_c1=double(s.y);  y_c1=y_c1(2);
    plot([xr1,x_c1],[yr1,y_c1],'Color',color,'Linewidth',1.5);
    kr12=pslope(x_c1,y_c1,x2,y2,k2,b2);
    [xr13,yr13,flag]=reflect(xr1,yr1,x_c1,y_c1,kr12,0,-10,55,65,i,1);
    kr13=(yr13-y_c1)/(xr13-x_c1);jr13=y_c1-kr13*x_c1;
    KJ1(1,i)=kr13;KJ1(2,i)=jr13;
    %grating初试2
    s=solve(((x-x2)+k2.*(y-y2)).^2-2.*k2.*b2.*(x-x2)+2.*b2.*(y-y2)-b2.^2==0, ...
        (xr2-xr1).*(y-yr2)==5*(yr2-yr1).*(x-xr2),x,y);
    x_c2=double(s.x);  x_c2= x_c2(2);
    y_c2=double(s.y); y_c2=y_c2(2);
    plot([xr1,x_c2],[yr1,y_c2],'Color',color,'Linewidth',1.5);
    kr22=pslope(x_c2,y_c2,x2,y2,k2,b2);
    [xr23,yr23,flag]=reflect(xr1,yr1,x_c2,y_c2,kr22,0,-10,55,65,i,1);
    kr23=(yr23-y_c2)/(xr23-x_c2);jr23=y_c2-kr23*x_c2;
    KJ2(1,i)=kr23;KJ2(2,i)=jr23;
end
[xfocus0,yfocus0]=linecross(KJ0(1,1),KJ0(2,1),KJ0(1,2),KJ0(2,2));
scatter(xfocus0,yfocus0,'b*');
[xfocus1,yfocus1]=linecross(KJ1(1,1),KJ1(2,1),KJ1(1,2),KJ1(2,2));
scatter(xfocus1,yfocus1,'b*');
[xfocus2,yfocus2]=linecross(KJ2(1,1),KJ2(2,1),KJ2(1,2),KJ2(2,2));
scatter(xfocus2,yfocus2,'b*');
d=((xfocus2-xfocus0)^2+(yfocus2-yfocus0)^2)^0.5
%% 主函数中用到的函数
%input:入射光路上一点x1y1 反射点x2y2 反射镜斜率km
% 目标线段(cd)的斜率ko和截距jo,端点横坐标CxDx,i用于选颜色,hua 1画 0不画
%output:画出从入射光路上一点出发,经一次反射,打到cd上的完整光路
% x3 y3反射光与cd所在直线的交点,flag表示能否打中
function [x3,y3,flag]=reflect(x1,y1,x2,y2,km,ko,jo,Cx,Dx,i,hua)
    coArray=['#0072BD';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#4DBEEE';'#A2142F'];%颜色库
    color=coArray(mod(i,7)+1,:); %根据i选定一种颜色
    if hua==1
    plot([x1,x2],[y1,y2],'Color',color,'Linewidth',1.5);%画入射光线
    end
    %km=0,法线斜率为Inf,故分类讨论
    if km==0
        x3=2*x2-x1;
        y3=y1;
    else
        %过x1y1作与镜子平行的直线 斜率km    截距y1-km.*x1
        %过x2y2作与镜子垂直的直线 斜率-1/km 截距y2+x2./km
        [x_c,y_c]=linecross(km,y1-km.*x1, -1/km,y2+x2./km);%求两直线交点
        x3=2.*x_c-x1;
        y3=2.*y_c-y1;
    end
    if abs(x2-x3)<=1e-10   %matlab算力有限
        y3=ko.*x2+jo;
    else
        kr=(y3-y2)/(x3-x2); %反射光斜率
        jr=y2-kr*x2;
        [x3,y3]=linecross(kr,jr,ko,jo);%反射光与cd的交点
    end
    if hua==1
        plot([x2,x3],[y2,y3],'Color',color,'Linewidth',1.5);%画反射光线
    end
    if x3>=Cx && x3<=Dx
        flag=1;%打中
    else
        flag=0;%没打中
    end
end
%% 求直线与抛物线的交点
%准线斜率k 截距b,入射光斜率m,焦点x0,y0,flag取哪一边的点1or2
%flag选1还是2，要试验，我们用的抛物面镜要么是抛物线对称轴左半部分的一小部分，
% 要么是右半部分的一小部分
function [x_c,y_c]=lpcross(k,b,m,x0,y0,flag)
    if flag==1
        x_c=(b.*(k-m)+b.*(k.^2+m.^2+k.^2.*m.^2+1).^0.5)/(k.*m+1)^2;
    elseif flag==2
        x_c=(b.*(k-m)-b.*(k.^2+m.^2+k.^2.*m.^2+1).^0.5)/(k.*m+1)^2;
    end
    y_c=m.*x_c;
    x_c=x_c+x0;
    y_c=y_c+y0;
end
%% 抛物线上一点x,y的斜率
%抛物线焦点x0,y0移到原点后的准线的斜率k与截距b
function k=pslope(x,y,x0,y0,k,b)
    x=x-x0;y=y-y0; %先把点移到原点抛物线上
    k=(k.*b-x-k.*y)./(k^2.*y+k.*x+b);%原点抛物线上一点的切线斜率
end
%% 求已知点关于已知直线的对称点
%x1,y1 已知点
% k b 对称轴斜率、截距
function [x2,y2]=symmetry(x1,y1,k,b)
    if k==-inf || inf
        
    end
    [x_c,y_c]=linecross(k,b,-1/k,y1+x1/k); %已知点和对称点的连线与对称轴的交点
    x2=2*x_c-x1; %中点在对称轴上
    y2=2*y_c-y1;
end
%推导: 1.中点在对称轴上 2. 连线与对称轴垂直
%% 已知两条直线的斜率和截距，求交点坐标
function [x,y]=linecross(k1,b1,k2,b2)
  x=[];
  y=[];
  if k1==k2 && b1==b2
      disp('重合');
  elseif k1==k2 && b1~=b2
      disp('无交点');
  else
     x=(b2-b1)/(k1-k2);
     y=k1*x+b1;
  end
end
%联立解方程