clear;clc
x0=0;y0=0;
R=304.8;
F=0.5*R;
t=linspace(-2.39/180*pi,2.39*pi/180,21); %直径25.4mm
center=[-R/2,-26];
X=R*cos(t)+center(1);
Y=R*sin(t)+center(2);
hold on;
% scatter(center(1),center(2),'r*');
scatter(x0,y0,'k*');
plot(X,Y,'r','Linewidth',2);
Pt=t([round(length(t)/2)-2,round(length(t)/2)-1,round(length(t)/2),round(length(t)/2)+1,round(length(t)/2)+2]);
Nn=numel(Pt);%光线数
Px=R.*cos(Pt)+center(1);
Py=R.*sin(Pt)+center(2);
kr=0;
%平面镜1
psi=-85.2;  % 可见光范围83.3-86.8  85.6802
Ax=-80;Ay=-53.27; %平面镜左端
length=25;
Bx=Ax+length*cosd(psi);By=Ay+length*sind(psi); %平面镜右端
km1=(By-Ay)/(Bx-Ax); %镜斜率
jm1=Ay-km1*Ax; %镜截距
plot([Ax,Bx],[Ay,By],'r','Linewidth',2);
%球面镜2
R2=R;
center2=[-R2/2,-105.297];
t2=t;
X2=R2*cos(t2)+center2(1);
Y2=R2*sin(t2)+center2(2);
plot(X2,Y2,'r','Linewidth',2);
%探测器
Cx=9.4;Cy=-131.77; %mirror2左
Dx=9.4;Dy=-127.77; %mirror2右
km2=(Dy-Cy)/(Dx-Cx); %镜斜率
jm2=Cy-km2*Cx; %镜截距
plot([Cx,Dx],[Cy,Dy],'r','Linewidth',2);

% %防杂散光1--光路
% kz=-1.15;jz=25;%y=-1.15x+15
% xslit=24.7578;yslit=21.52853;
% xzuo=22.95806;yzuo=xzuo*kz+jz;
% xyou=27.11955;yyou=xyou*kz+jz;
% Fx=13;Gx=33;
% scatter(xslit,yslit,'*k');
% plot([xzuo,xslit,xyou],[yzuo,yslit,yyou],'k','Linewidth',1.5)
% plot([Fx,Gx],[Fx*kz+jz,Gx*kz+jz],'r','Linewidth',2)
% 放杂散光2--挡板
da=0;db=0;
plot([0,40+da,40+da,nan,40+da,40+da],[-24+db,-24+db,-12+db,nan,0+db,21],'r-','Linewidth',2)
plot([90+da,230+da],[-65.65+db,-65.65+db],'r:','Linewidth',2)
plot([-10+da,25+da],[-110+db,-110+db],'r:','Linewidth',2)
syms x y;
for i=1:Nn
    kr=-(Px(i)-center(1))/(Py(i)-center(2));
    [x1,y1,flag]=reflect(x0,y0,Px(i),Py(i),kr,km1,jm1,Ax,Bx,1); %光源0--球一(i)--光栅1
    if flag==0
        continue;
    end

    %光栅
    ki=(y1-Py(i))/(x1-Px(i));
    phi=atand(ki)-psi-90;
    d=0.001/300; %3.33μm 即300条/mm
    %波长    紫     绿        红
    lambda=[400e-9 520e-9 720e-9];%j 遍历
%     lambda=[400e-9 475e-9 550e-9 620e-9 740e-9]; %j 遍历
    theta_z1=asind(sind(phi)-lambda/d); %干涉+1级
    theta_f1=asind(sind(phi)+lambda/d) %干涉-1级
    theta_shan=4.3; %闪耀角
    lambda_shan=2.*d.*cosd(phi+theta_shan).*sind(theta_shan); %闪耀波长497nm左右
    %闪耀光线方向
    theta_0=phi+2*theta_shan;
    
    %0级
    [x2,y2,flag]=reflect(Px(i),Py(i),x1,y1,km1,Inf,Inf,20,20,0); %球一(i)--光栅1--2(算反射斜率)
    kg_0=tand(psi+90-theta_0); %0级衍射斜率 %闪耀角psi+90-theta_0 %正常psi+90-phi
    if kg_0^2==Inf
        s=solve((x-center2(1)).^2+(y-center2(2)).^2==R2^2,x==x1,x,y);
    else
        s=solve((x-center2(1)).^2+(y-center2(2)).^2==R2^2, ...
            (y-y1)==kg_0.*(x-x1),x,y);
    end
    x_c0=double(s.x);
    y_c0=double(s.y);
    if numel(x_c0)>1
        if x_c0(1)>0
            x_c0=x_c0(1);y_c0=y_c0(1);
        else
            x_c0=x_c0(2);y_c0=y_c0(2);
        end
    end
    if y_c0<=-91 && y_c0>=-120
        [x4,y4,flag]=reflect(x1,y1,x_c0,y_c0,-(x_c0-center2(1))/(y_c0-center2(2)),km2,jm2,9.2,9.2,1,lambda_shan);
    end
    
    for j=1:numel(lambda)
        %grating -1
        kg_f1=tand(psi+90-theta_f1(j)); %衍射光斜率
        if kg_f1^2==Inf
            s=solve((x-center2(1)).^2+(y-center2(2)).^2==R2^2,x==x1,x,y);
        else
            s=solve((x-center2(1)).^2+(y-center2(2)).^2==R2^2, ...
                (y-y1)==kg_f1.*(x-x1),x,y);
        end
        x_cf1=double(s.x);
        y_cf1=double(s.y);
        if numel(x_cf1)>1
            if x_cf1(1)>0
                x_cf1=x_cf1(1);y_cf1=y_cf1(1);
            else
                x_cf1=x_cf1(2);y_cf1=y_cf1(2);
            end
        end
        if y_cf1<=-91 && y_cf1>=-120
            [x4,y4,flag]=reflect(x1,y1,x_cf1,y_cf1,-(x_cf1-center2(1))/(y_cf1-center2(2)),km2,jm2,9.2,9.2,1,lambda(j));
        end
    end
end
axis equal
set(gca,'color','none');
axis off;
hold off;
%% 主函数中用到的函数
%input:入射光路上一点x1y1 反射点x2y2 反射镜斜率km
% 目标线段(cd)的斜率ko和截距jo,端点横坐标CxDx,i用于选颜色,
% hua 1画(画图可在后面加i参数) 0不画 2最强光真颜色，3最强光黑色
%output:画出从入射光路上一点出发,经一次反射,打到cd上的完整光路
% x3 y3反射光与cd所在直线的交点,flag表示能否打中
function [x3,y3,flag]=reflect(x1,y1,x2,y2,km,ko,jo,Cx,Dx,hua,lambda)
    if Cx>Dx
        cddc=Cx;
        Cx=Dx;
        Dx=cddc;
    end
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
    if hua~=0
        if(~exist('lambda','var'))||(hua==3)
            color='k';
        else
            color=double(py.pycode.pyrgb.getRGB(lambda*1e9)); %lambda2RGB
        end
        if (hua==2) || (hua==3)
            plot([x1,x2],[y1,y2],'Color',color,'Linewidth',0.5,'LineStyle','--');%画入射光线
            plot([x2,x3],[y2,y3],'Color',color,'Linewidth',0.5,'LineStyle','--');%画反射光线
        else
            plot([x1,x2],[y1,y2],'Color',color,'Linewidth',0.5);%画入射光线
            plot([x2,x3],[y2,y3],'Color',color,'Linewidth',0.5);%画反射光线
        end
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