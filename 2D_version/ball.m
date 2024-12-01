x0=0;y0=0;
axis equal
R=304.8;
F=0.5*R;
t=88.6/180*pi:0.005:92.4*pi/180;
center=[13 -R/2];
X=R*cos(t)+center(1);
Y=R*sin(t)+center(2);
hold on;
% scatter(center(1),center(2),'r*');
scatter(x0,y0,'k*');
plot(X,Y,'r','Linewidth',2);
Nn=3;%光线数
Pt=t([1,round(length(t)/2),end]);
Px=R.*cos(Pt)+center(1);
Py=R.*sin(Pt)+center(2);
kr=0;
%平面镜1
Ax=10;Ay=0; %平面镜左端
Bx=40;By=0; %平面镜右端
km1=(By-Ay)/(Bx-Ax); %镜斜率
jm1=Ay-km1*Ax; %镜截距
plot([Ax,Bx],[Ay,By],'k','Linewidth',2);
%球面镜2
R2=304.8;
center2=[39,-R2/2];
t2=t;
X2=R2*cos(t2)+center2(1);
Y2=R2*sin(t2)+center2(2);
plot(X2,Y2,'r','Linewidth',2);
%探测器
Cx=15.5;Cy=0; %mirror2左
Dx=16.5;Dy=0; %mirror2右
km2=(Dy-Cy)/(Dx-Cx); %镜斜率
jm2=Cy-km2*Cx; %镜截距
plot([Cx,Dx],[Cy,Dy],'k','Linewidth',2);
% x3=nan;y3=nan;
syms x y;
for i=1:Nn
    coArray=['#0072BD';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#4DBEEE';'#A2142F'];%颜色库
    color=coArray(mod(i,7)+1,:); %根据i选定一种颜色
%     x3=nan;y3=nan;
    kr=(Py(i)-center(2))/(Px(i)-center(1));
    [x1,y1,flag]=reflect(x0,y0,Px(i),Py(i),-1/kr,km1,jm1,Ax,Bx,1);
    if flag==0
        continue;
    end
    [x2,y2,flag]=reflect(Px(i),Py(i),x1,y1,km1,0,50,40,50,0);
    s=solve((x-center2(1)).^2+(y-center2(2)).^2==R2^2, ...
        (x1-x2).*(y-y2)==(y1-y2).*(x-x2),x,y);
    x3=double(s.x); x3=x3(2);
    y3=double(s.y); y3=y3(2);
    plot([x1,x3],[y1,y3],'Color',color,'Linewidth',1.5);
    [x4,y4,flag]=reflect(x1,y1,x3,y3,-(x3-center2(1))/(y3-center2(2)),km2,jm2,Cx,Dx,1);
end
axis equal
set(gca,'color','none');
axis off;
hold off;
axis equal
%% 主函数中用到的函数
%input:入射光路上一点x1y1 反射点x2y2 反射镜斜率km
% 目标线段(cd)的斜率ko和截距jo,端点横坐标CxDx,i用于选颜色,hua 1画(画图可在后面加i参数) 0不画
%output:画出从入射光路上一点出发,经一次反射,打到cd上的完整光路
% x3 y3反射光与cd所在直线的交点,flag表示能否打中
function [x3,y3,flag]=reflect(x1,y1,x2,y2,km,ko,jo,Cx,Dx,hua,i)
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
    %反射光线斜率有可能Inf
    if abs(x2-x3)<=1e-10   %matlab算力有限(可以理解为x2==x3)
        y3=ko.*x2+jo;
    else
        kr=(y3-y2)/(x3-x2); %反射光斜率
        jr=y2-kr*x2;
        [x3,y3]=linecross(kr,jr,ko,jo);%反射光与cd的交点
    end
    %判断是否打中
    if x3>=Cx && x3<=Dx
        flag=1;
    else
        flag=0;
    end    
    %是否画图
    if hua==1
        if(~exist('i','var'))
            color='k';
        else
            coArray=["#7E2F8E";"#0072BD";"#77AC30"];%颜色库
            color=coArray(mod(i-1,3)+1); %根据i选定一种颜色
        end
        plot([x1,x2],[y1,y2],'Color',color,'Linewidth',1.5);%画入射光线
        plot([x2,x3],[y2,y3],'Color',color,'Linewidth',1.5);%画反射光线
    end
end
%% 抛物线上一点x,y的斜率
%抛物线焦点x0,y0移到原点后的准线的斜率k与截距b
function k=pslope(x,y,x0,y0,k,b)
    x=x-x0;y=y-y0; %先把点移到原点抛物线上
    k=(k.*b-x-k.*y)./(k^2.*y+k.*x+b);%原点抛物线上一点的切线斜率
end
%% 已知两条直线的斜率和截距，求交点坐标
function [x,y]=linecross(k1,b1,k2,b2)
  x=[];
  y=[];
  if k1==k2 && b1==b2
      disp('重合');
  elseif k1==k2 && b1~=b2
      disp('无交点');
  else
     x=(b2-b1)./(k1-k2);
     y=k1.*x+b1;
  end
end
%联立解方程