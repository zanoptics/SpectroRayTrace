clear;
x0=0;y0=0;
R=304.8;
F=0.5*R;
t=-2.388/180*pi:0.001*pi/180:2.388*pi/180;
center=[-R/2,-13];
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
Ax=10;Ay=-10; %平面镜左端
lm=52;
psi=-55.5;  %0级-42.5622   -1绿光 -55.5
Bx=Ax+lm*cosd(psi);By=Ay+lm*sind(psi); %平面镜右端
km1=(By-Ay)/(Bx-Ax); %镜斜率
jm1=Ay-km1*Ax; %镜截距
plot([Ax,Bx],[Ay,By],'r','Linewidth',2);


%抛物面镜
xp=127.1196;yp=30;
scatter(xp,yp,'k*');
%准线斜率截距(普通方程需要的量)
kp=0;bp=101.6; 
%换算为 参数方程 需要的量
p=bp/sqrt(kp^2+1); %(0,0)代入距离公式
phip=atand(kp)-90; %旋转角度 
tp=-0.563:0.001:-0.438;
xtp=2*p.*tp.^2-p/2;
ytp=2*p.*tp;
%对每个点,绕原点旋转变换
%向量化
temp=zeros(2,numel(xtp));
temp=[cosd(phip),-sind(phip);sind(phip),cosd(phip)]*[xtp;ytp];%旋转矩阵
xtp=temp(1,:);ytp=temp(2,:);
%向量化
plot(xtp+xp,ytp+yp,'r-','LineWidth',2);%抛物线

%探测器
Cx=155;Cy=29; %mirror2左
Dx=155;Dy=31; %mirror2右
km2=(Dy-Cy)/(Dx-Cx); %镜斜率
jm2=Cy-km2*Cx; %镜截距
% plot([Cx,Dx],[Cy,Dy],'k','Linewidth',2);

%遮挡
% plot([0,150],[0,20])
plot([127.12,127.12],[34,31],'r-','LineWidth',1);
plot([127.12,127.12],[29,26],'r-','LineWidth',1);
%光栅方程
%算psi 0.085301837270341
ki=0.085301837;
phi=90-atand(ki)+psi;
d=0.001/600; %1.67μm 即600条/mm
%波长    紫     蓝         绿
lambda=[400e-9 475e-9 550e-9 625e-9];
theta_z1=asind(sind(phi)-lambda/d);
theta_f1=asind(sind(phi)+lambda/d);
theta0=9.5;
syms x y;
for i=1:Nn
    kr=-(Px(i)-center(1))/(Py(i)-center(2));
    [x1,y1,flag]=reflect(x0,y0,Px(i),Py(i),kr,km1,jm1,Ax,Bx,1);
    if flag==0
        continue;
    end
    [x2,y2,flag]=reflect(Px(i),Py(i),x1,y1,km1,0,40,0,150,0);
    %求交点
    kg0=(y2-y1)/(x2-x1);  %衍射光斜率
    kg0=tand(atand(kg0)+2*theta0);
    if kg0^2==Inf
        s=solve((x-xp).^2+2.*p.*(y-yp-p/2)==0,x==x1,x,y);
    else
        s=solve((x-xp).^2+2.*p.*(y-yp-p/2)==0, ...
            (y-y1)==kg0.*(x-x1),x,y);
    end
    x_c0=double(s.x); 
    y_c0=double(s.y);
    if numel(x_c0)>1
        if abs(x_c0(1)-xp)<abs(x_c0(2)-xp)
            x_c0=x_c0(1);y_c0=y_c0(1);
        else
            x_c0=x_c0(2);y_c0=y_c0(2);
        end
    end
    k_c0=(xp-x_c0)/p; %抛物线斜率
    [x_d0,y_d0,flag]=reflect(x1,y1,x_c0,y_c0,k_c0,km2,jm2,Cx,Dx,1);

    for j=1:numel(lambda)
        %grating -1
        kg_f1=tand(90+theta_f1(j)+psi); %衍射光斜率
        if kg_f1^2==Inf
            s=solve((x-xp).^2+2.*p.*(y-yp-p/2)==0,x==x1,x,y);
        else
            s=solve((x-xp).^2+2.*p.*(y-yp-p/2)==0, ...
            y-y1==kg_f1.*(x-x1),x,y);
        end
        x_cf1=double(s.x);  
        y_cf1=double(s.y);  
        if numel(x_cf1)>1
            if abs(x_cf1(1)-xp)<abs(x_cf1(2)-xp)
                x_cf1=x_cf1(1);y_cf1=y_cf1(1);
            else
                x_cf1=x_cf1(2);y_cf1=y_cf1(2);
            end
        end
        k_cf1=(xp-x_cf1)/p; %抛物线斜率
        [x_df1,y_df1,flag]=reflect(x1,y1,x_cf1,y_cf1,k_cf1,km2,jm2,Cx,Dx,1,j);

%         %grating +1
%         kg_z1=tand(90+theta_z1(j)+psi); %衍射光斜率
%         if kg_z1^2==Inf
%             s=solve((x-xp).^2+2.*p.*(y-yp-p/2)==0,x==xr1,x,y);
%         else
%             s=solve((x-xp).^2+2.*p.*(y-yp-p/2)==0, ...
%             y-y1==kg_z1.*(x-x1),x,y);
%         end
%         x_cz1=double(s.x);  
%         y_cz1=double(s.y);  
%         if numel(x_cz1)>1
%             if abs(x_cz1(1)-xp)<abs(x_cz1(2)-xp)
%                 x_cz1=x_cz1(1);y_cz1=y_cz1(1);
%             else
%                 x_cz1=x_cz1(2);y_cz1=y_cz1(2);
%             end
%         end
%         k_cz1=(xp-x_cz1)/p; %抛物线斜率
%         [x_dz1,y_dz1,flag]=reflect(x1,y1,x_cz1,y_cz1,k_cz1,km2,jm2,Cx,Dx,1,j);
    end
end
axis equal
hold off;
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