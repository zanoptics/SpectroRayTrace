%% 主函数
tic %开始计时
x0=0;y0=90; %光源
Ax=25;Ay=10; %平面镜左端
Bx=60;By=12; %平面镜右端
km1=(By-Ay)/(Bx-Ax); %镜斜率
jm1=Ay-km1*Ax; %镜截距

Cx=50;Cy=70; %mirror2左
Dx=90;Dy=80; %mirror2右
km2=(Cy-Dy)/(Cx-Dx);
jm2=Cy-km2*Cx; 

Ex=90;Ey=40;%mirror3左
Fx=120;Fy=44.5;
km3=(Fy-Ey)/(Fx-Ex);
jm3=Fy-km3*Fx;

Gx=118;Gy=80;%detector左
Hx=133;Hy=70;
kd=(Hy-Gy)/(Hx-Gx);
jd=Gy-kd*Gx;

%画
hold on; %画在同一张
daspect([1,1,1]); %设置比例尺

%两点式
M=plot([Ax,Bx],[Ay,By],'k','Linewidth',2);%平面镜1
plot([Cx,Dx],[Cy,Dy],'k','Linewidth',2);%平面镜2
% plot([Ex,Fx],[Ey,Fy],'k','Linewidth',2);%平面镜3
D=plot([Gx,Hx],[Gy,Hy],'r','Linewidth',2);%探测器
S=scatter(x0,y0,'r','*'); %画出光源
Nn=4; %入射光线数

z=Ax:1:Bx; %设置离散点
%在平面镜上随机取Nn个点
Px=z(randperm(numel(z),Nn)); 
Py=km1*Px+jm1;

for i=1:Nn
    [x3,y3,flag]=reflect(x0,y0,Px(i),Py(i),km1,km2,jm2,Cx,Dx,i);%第一次反射
    if flag==0
        continue; %flag=0 表示没打中,不会执行接下来的反射
    end
    
end
title('反射定律可视化','FontWeight','bold')%标题加粗
legend([M D S],'mirror','detector','source','Location','southeast');%图例locte in东南
hold off;
toc %停止计时
%% 主函数中用到的函数
%input:入射光路上一点x1y1 反射点x2y2 反射镜斜率km
% 目标线段(cd)的斜率ko和截距jo,端点横坐标CxDx,i用于选颜色
%output:画出从入射光路上一点出发,经一次反射,打到cd上的完整光路
% x3 y3反射光与cd所在直线的交点,flag表示能否打中
function [x3,y3,flag]=reflect(x1,y1,x2,y2,km,ko,jo,Cx,Dx,i)
    coArray=['#0072BD';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#4DBEEE';'#A2142F'];%颜色库
    color=coArray(mod(i,7)+1,:); %根据i选定一种颜色
    if Cx>Dx %保证Cx较小
        t=Cx;
        Cx=Dx;
        Dx=t;
    end
    plot([x1,x2],[y1,y2],'Color',color,'Linewidth',1.5);%画入射光线
    [x3,y3]=symmetry(x1,y1,-1/km,y2+x2/km); %关于法线的对称点 %法线的斜率和截距
    kr=(y3-y2)/(x3-x2); %反射光斜率
    jr=y2-kr*x2;
    [x3,y3]=linecross(kr,jr,ko,jo);%反射光与cd的交点
    if x3>=Cx && x3<=Dx
        plot([x2,x3],[y2,y3],'Color',color,'Linewidth',1.5);%画反射光线
        flag=1;%打中
    else
        xx=1.03*x3;  %取一个稍远一点的点,表示没打中
        yy=kr*xx+jr; 
        plot([x2,xx],[y2,yy],'Color',color,'Linewidth',1.5);
        flag=0;%没打中
    end
end
%% 求已知点关于已知直线的对称点
%x1,y1 已知点
% k b 对称轴斜率、截距
function [x2,y2]=symmetry(x1,y1,k,b)
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