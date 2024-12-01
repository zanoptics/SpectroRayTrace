%==========================================================================
%光线追迹 多束光线经过一个抛物面镜,反射到探测器上
%22.9.22
%汪雪
%% ==========================================================================
%定义光源S、抛物面镜M，探测器D
%==========================================================================
clear all;
hold on;
%抛物线:y=a(x-b)^2+c
% x=-5:0.001:25;
x=-100:0.001:-50;
a=-1/300; b=10; c=100;
y=a.*(x-b).^2+c;%抛物线方程
f=[b 1/(4*a)+c];%焦点坐标
x0=b; y0=1/(4*a)+c; %定义光源坐标
% x0=-20; y0=1/(4*a)+c; %定义光源坐标
%准线

% sita=30; %旋转角度（角度制）
sita=-10; %旋转角度（角度制）
%绕焦点f顺时针旋转α度
% T=[1 0 x0;0 1 y0;0 0 1];%平移矩阵  
% R=[cosd(alpha) sind(alpha) 0;sind(alpha) cosd(alpha) 0;0 0 1];%旋转矩阵

X=(x-f(1))*cosd(sita)+(y-f(2))*sind(sita)+f(1);
Y=-(x-f(1))*sind(sita)+(y-f(2))*cosd(sita)+f(2); %顺时针转35°
z=-50:0.001:50;
y1=tand(90-sita).*(z-f(1))+f(2); %抛物线对称轴 斜率=1.4281(35°)
Cx=-300; Cy=50; %C点坐标
Dx=-50; Dy=0; %D点坐标
%D=(Dy-Cy)/(Dx-Cx)*(x-Cx)+Cy;%用C、D点确定ccd探测面方程
kd=(Dy-Cy)/(Dx-Cx);%探测器方程斜率
jd=Cy-kd*Cx;%探测器方程截距
D=kd*x+jd;%探测器方程

%% ==========================================================================
%生成并画出多条入射光线，画出球面镜及探测器
%==========================================================================
% plot(x,y,'k','Linewidth',2);%画出抛物面镜
% hold on;
% axis equal
plot(X,Y,'b','Linewidth',2);
hold on;
% axis equal
daspect([1 1 1]);
plot(z,y1,'--r','Linewidth',1); %抛物线对称轴
plot(z,y1+90,'--r','Linewidth',1); %辅助线
scatter(f(1),f(2),'o','r');
plot([Cx,Dx],[Cy,Dy],'r','Linewidth',2);%画出探测器平面
nn=5;%随机取的光线个数
px=x(randperm(numel(x),nn)); %在AB上随机取五个点，作为光源出射光线与球面镜的交点
py=a.*(px-b).^2+c;
Px=(px-f(1))*cosd(sita)+(py-f(2))*sind(sita)+f(1);
Py=-(px-f(1))*sind(sita)+(py-f(2))*cosd(sita)+f(2);
for i=1:nn
   plot([x0,Px(i)],[y0,Py(i)],'Linewidth',1);%画出入射光线
end

%% ==========================================================================
%求入射光线斜率，求交点切线斜率，并求出反射光线斜率
%2*km=k1+k2  kr(i)=(-ki(i)*km^2-2*km+ki(i))/(km^2-2*km*ki(i)-1);
%求解反射光线与探测面交点，画出反射光线
%==========================================================================
ki=zeros(1,nn);%入射光线斜率
kr=zeros(1,nn);%反射光线斜率
jr=zeros(1,nn);%反射光线截距
Qx=zeros(1,nn);%反射光线与探测面交点
Qy=zeros(1,nn);
for i=1:nn
%     ki(i)=(Py(i)-y0)/(Px(i)-x0);
    %%对抛物线求导，可得出该点切线斜率y'=2*a*(x-b)
%     k=2*a*(Px(i)-b);%旋转前切线斜率
%     km=-tand(tand(k)-sita);
      %对Y求对x的偏导，来求解求线斜率dY(x,y)/dx=dY/dx+dY/dy*dy/dx;
      %Y=-(x-f(1))*sind(sita)+(y-f(2))*cosd(sita)+f(2);
      %km=-sind(sita)+cosd(sita)*(2*a*(px(i)-b));   
      km=tand(atand(2*a*(px(i)-b))-sita);      
%     kr(i)=2*km- ki(i);
%     kr(i)=(-ki(i)*km^2-2*km+ki(i))/(km^2-2*km*ki(i)-1);
    kr(i)=reflect(x0,y0,Px(i),Py(i),km);
    jr(i)=Py(i)-kr(i)*Px(i);%反射光线截距
    [Qx(i),Qy(i)]=linecross(kr(i),jr(i),kd,jd);%求反射光线与探测面交点 
    if Qx(i)>=min(Cx,Dx) & Qx(i)<=max(Cx,Dx)           
        plot([Px(i),Qx(i)],[Py(i),Qy(i)],'Linewidth',1);%画出反射光线
    end
end
hold off;

%% 求入射角=反射角的情况下的反射光线斜率
function k=reflect(x1,y1,x2,y2,km)
    %法线方程：y=-1/km*(x-x2)+y2
    kf=-1/km; jf=1/km*x2+y2;
    %过光源且垂直于法线的直线方程：y=km*(x-x1)+y1
    k1=km; j1=-km*x1+y1;
    %求关于法线与光源对称的点坐标
    [xm,ym]=linecross(kf,jf,k1,j1);%中点
    x3=2*xm-x1;
    y3=2*ym-y1;
    k=(y3-y2)/(x3-x2);
end

%% 已知两条直线的斜率和截距，求交点坐标
function [x,y]=linecross(k1,b1,k2,b2)
  x=[];
  y=[];
  if k1==k2&b1==b2
      disp('重合');
  elseif k1==k2&b1~=b2
      disp('无交点');
  else
     x=(b2-b1)/(k1-k2);
     y=k1*x+b1;
  end
end

