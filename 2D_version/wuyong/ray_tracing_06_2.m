%==========================================================================
%光线追迹 平行光经过一个反射光栅,反射到探测器上
%22.9.21
%汪雪
%% ==================================================================
%定义光源S、光栅G，探测器D
%==========================================================================
clear all;
alpha=60; %平行光入射的角度（角度制）  
ki=-tand(90-alpha);%入射光斜率
lamda=6.33e-4; %入射光红光波长为633nm
% lamda=5.32e-4; %入射光绿光波长为532nm
% d=0.8;%光栅的周期为d=800um
d=0.006;%光栅的周期为d=90um
Ax=20; Ay=20;
Bx=60; By=20; %光栅坐标
x=Ax:0.001:Bx;
kg=(By-Ay)/(Bx-Ax);%光栅斜率及截距
jg=Ay-kg*Ax;

Cx=20; Cy=60; %C点坐标
Dx=250; Dy=60; %D点坐标
kd=(Dy-Cy)/(Dx-Cx);%探测器方程斜率
jd=Cy-kd*Cx;%探测器方程截距

%% =================================================================
%生成并画出多条入射光线，画出光栅及探测器
%==========================================================================
plot([Ax,Bx],[Ay,By],'r','Linewidth',2)%画出光栅G
hold on;
daspect([1 1 1]);
plot([Cx,Dx],[Cy,Dy],'r','Linewidth',2);%画出探测器平面

nn=10;%随机取的光线个数
Px=x(randperm(numel(x),nn));  %在AB上随机取五个点，作为光源出射光线与抛面镜的交点
Py=kg.*Px+jg;
for i=1:nn
%     y=ki*(x-Px(i))
    yy=50;
    xx=yy/ki+Px(i);
   plot([xx,Px(i)],[yy,Py(i)],'k','Linewidth',1);%画出入射光线
end

%% =================================================================
%计算入射角
%mλ=d(sinθm?sinθi) θi为入射角，θm为m级反射角
%生成并画出多条反射光线
%==========================================================================
mm=3; %画0到3级反射光线
color=['m';'y';'g';'b'];
for j=0:mm
    for i=1:nn
        re=asin(j*lamda/d+sind(alpha));%出射角
        if ki>0
            kr=-tan(re-atan(-1/kg));
%            kr=-tan((pi/2-re)-atan(kg));%反射光线方程斜率
        else
%             kr=tan(re+atan(-1/kg));
           kr=tan(atan(kg)+(pi/2-re));%反射光线方程斜率
        end
        jr=Py(i)-kr*Px(i);%反射光线截距
        [X,Y]=linecross(kr,jr,kd,jd);%求反射光线与探测面交点 
        if X>=min(Cx,Dx) & X<=max(Cx,Dx)           
            plot([Px(i),X],[Py(i),Y],'Color',color(j+1),'Linewidth',1);%画出反射光线
        end 
    end
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