clear;
hold on;
%焦点
x0=0;y0=0;
scatter(x0,y0,'k*');
%准线斜率截距(普通方程需要的量)
k0=0;b0=50.8;
%换算为 参数方程 需要的量
p0=b0/sqrt(k0^2+1); %(0,0)代入距离公式
phi0=atand(k0)-90; %旋转角度 
t0=0.35:0.00001:0.6;
xp0=2*p0.*t0.^2-p0/2;
yp0=2*p0.*t0;
%对每个点,绕原点旋转变换

%向量化
temp=zeros(2,numel(xp0));
temp=[cosd(phi0),-sind(phi0);sind(phi0),cosd(phi0)]*[xp0;yp0];%旋转矩阵
xp0=temp(1,:);yp0=temp(2,:);
%向量化
plot(xp0+x0,yp0+y0,'r-','LineWidth',2);%抛物线

%参数调节区
psi=31.5039; %平面镜绕A旋转角度ψ
Nn=3;%光线数
Pt=t0([1,round(length(t0)/2),end]); %1 7 end
Px=2.*p0.*Pt.^2-p0/2+x0;
Py=2.*p0.*Pt+y0;

%对每个点,绕原点旋转变换
%向量化
temp=zeros(2,numel(Px));
temp=[cosd(phi0),-sind(phi0);sind(phi0),cosd(phi0)]*[Px;Py];%旋转矩阵
Px=temp(1,:);Py=temp(2,:);
%向量化
scatter(Px,Py,'k*');
%平面镜1
Ax=35.55;Ay=-36.227; %平面镜左端
% plot([Ax,33.8913],[Ay,49.9434],'r--');辅助红虚线
lAB=32;
Bx=Ax+lAB*cosd(psi);
By=Ay+lAB*sind(psi);
km1=(By-Ay)/(Bx-Ax); %镜斜率
jm1=Ay-km1*Ax; %镜截距
plot([Ax,Bx],[Ay,By],'r','Linewidth',2);
%抛物面镜2
%焦点
x2=16;y2=-(19.92+50.8+5);scatter(x2,y2,'k*')
p2=p0;
t2=t0;
xp2=2*p2.*t2.^2-p2/2;
yp2=2*p2.*t2;
plot(xp2+x2,yp2+y2,'r-','LineWidth',2);%抛物线

%光栅方程
phi=psi;
d=0.001/600; %1.67μm 即600条/mm
%波长    紫     蓝         绿
lambda=[400e-9 475e-9 550e-9 625e-9];
% theta_z1=asind(sind(phi)-lambda/d);
theta_f1=asind(sind(phi)+lambda/d);
theta0=9.5; %闪耀角 衍射后的黑线就是单槽衍射的中央主极大方向(但颜色用复色光的黑色不恰当)
%有颜色的线表示不同波长的第一级干涉主极大
%原来的单槽衍射主极大就是符合反射定律的反射光线,现在要加个2*theta0
%其能量最大(可能不一定),且其附近的光能量都较强(一定)
%具体见第一个solve函数上面一行
syms x y
for i=1:Nn
    kr0=pslope(Px(i),Py(i),x0,y0,k0,b0);
    [xr1,yr1,flag]=reflect(x0,y0,Px(i),Py(i),kr0,km1,jm1,Ax,Bx,1);
    if flag==0
        continue;
    end
    [xr2,yr2,flag]=reflect(Px(i),Py(i),xr1,yr1,km1,km1,jm1+5,-inf,inf,0); %注意检查
    %求交点
    kg0=(yr2-yr1)/(xr2-xr1);  %衍射光斜率
    kg0=tand(atand(kg0)+2*theta0);
    s=solve((y-y2).^2-2.*p2.*(x-x2+p2./2)==0, ...
        (y-yr1)==kg0.*(x-xr1),x,y);
    x_c0=double(s.x); 
    y_c0=double(s.y);
    if(kg0 ~=0)
        if abs(x_c0(1))<abs(x_c0(2))
            x_c0=x_c0(1);y_c0=y_c0(1);
        else
            x_c0=x_c0(2);y_c0=y_c0(2);
        end
    end
    kr02=p2/(y_c0-y2); %抛物线斜率
    [xr03,yr03,flag]=reflect(xr1,yr1,x_c0,y_c0,kr02,0,-90,-inf,inf,1);

    for j=1:numel(lambda)
        %grating -1
        kg_f1=tand(90+theta_f1(j)+psi); %衍射光斜率
        s=solve((y-y2).^2-2.*p2.*(x-x2+p2./2)==0, ...
            (y-yr1)==kg_f1.*(x-xr1),x,y);
        x_c1=double(s.x);  
        y_c1=double(s.y);  
        if(kg_f1~=0)
            if abs(x_c1(1))<abs(x_c1(2))
                x_c1=x_c1(1);y_c1=y_c1(1);
            else
                x_c1=x_c1(2);y_c1=y_c1(2);
            end
        end
        kr12=p2/(y_c1-y2); %抛物线斜率
        [xr13,yr13,flag]=reflect(xr1,yr1,x_c1,y_c1,kr12,0,-90,55,65,1,j);
    end
end
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
    elseif km^2==Inf %平面镜斜率为Inf
        x3=x1;
        y3=2*y2-y1;
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
        %接受屏斜率±Inf
        if ko^2 == Inf
            x3=Cx;
            y3=kr.*x3+jr;
        else
            [x3,y3]=linecross(kr,jr,ko,jo);%反射光与cd的交点
        end
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
            coArray=["#7E2F8E";"#0072BD";"#77AC30";"#D95319"];%颜色库
            color=coArray(mod(i-1,4)+1); %根据i选定一种颜色
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