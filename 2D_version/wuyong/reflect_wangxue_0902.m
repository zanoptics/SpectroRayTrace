%==========================================================================
%光线追迹 多束光线经过一个平面反射镜,反射到探测器上
%22.8.31
%汪雪
%==========================================================================
%定义光源S、平面反射镜M，探测器D
%==========================================================================
% x=-100:0.0001:100;
tic

% x0=0;
% y0=100;%设置光源坐标S(x0,y0)
x0=-30;
y0=0;%设置光源坐标S(x0,y0)

% Ax=0; Ay=30; %A点坐标
% % Ax=0; Ay=10; %A点坐标
% Bx=50; By=20;%B点坐标
Ax=20; Ay=60; %A点坐标
Bx=70; By=40;%B点坐标
%M=(By-Ay)/(Bx-Ax)*(x-Ax)+Ay;%用A、B点确定平面反射镜方程
km=(By-Ay)/(Bx-Ax);%平面反射镜方程斜率
jm=Ay-km*Ax;%平面反射镜方程截距
z=Ax:0.001:Bx;
M=km*z+jm;%平面反射镜方程

% Cx=90; Cy=40; %C点坐标
% Cx=90; Cy=90; %C点坐标
% Dx=40; Dy=80;%D点坐标
Cx=100; Cy=-20; %C点坐标
Dx=40; Dy=-40; %D点坐标
%D=(Dy-Cy)/(Dx-Cx)*(x-Cx)+Cy;%用C、D点确定ccd探测面方程
kd=(Dy-Cy)/(Dx-Cx);%探测器方程斜率
jd=Cy-kd*Cx;%探测器方程截距

%==========================================================================
%生成并画出多条入射光线，画出平面镜及探测器
%==========================================================================
nn=10;%随机取的光线个数
Px=z(randperm(numel(z),nn));  %在AB上随机取五个点，作为光源出射光线与平面反射镜的交点
Py=km.*Px+jm; %在AB上随机取五个点，作为光源出射光线与平面反射镜的交点
hold on;
for i=1:nn
   plot([x0,Px(i)],[y0,Py(i)],'Linewidth',1);%画出5条入射光线
end
daspect([1 1 1]);
plot([Ax,Bx],[Ay,By],'k','Linewidth',2);%画出平面反射镜
plot([Cx,Dx],[Cy,Dy],'r','Linewidth',2);%画出探测器平面

%==========================================================================
%求解入射角，求出反射光线斜率及截距
%求解反射光线与探测面交点，画出反射光线
%==========================================================================
alpha=zeros(1,nn);%入射角余角
inan=zeros(1,nn);%入射角
for i=1:nn
    l=norm([x0-Ax,y0-Ay]);%求SA的模
    m=norm([x0-Px(i),y0-Py(i)]);%求SP_i的模
    n=norm([Ax-Px(i),Ay-Py(i)]);%求AP_i的模
    alpha(i)=acos((m^2+n^2-l^2)/(2*m*n));%用余弦定理求解入射角的余角
    inan(i)=(pi/2-alpha(i))/pi*180;%入射角，入射角=反射角
end
% inan %输出入射角

k_AD=(Dy-Ay)/(Dx-Ay);
k_AC=(Cy-Ay)/(Cx-Ax);%探测器可探测到的反射光线斜率范围为k_AC<k<k_AD

kr=zeros(1,nn);%反射光线斜率
jr=zeros(1,nn);%反射光线截距
X=zeros(1,nn);%反射光线与探测面交点
Y=zeros(1,nn);

% for i=1:nn
%     k1=(Py(i)-y0)/(Px(i)-x0);%求入射光线的斜率
%     if k1>0
%        kr(i)=-tan(alpha(i)-atan(km));%反射光线方程斜率
%     else
%        kr(i)=tan(atan(km)+alpha(i));%反射光线方程斜率
%     end
%        jr(i)=Py(i)-kr(i)*Px(i);%反射光线截距
%      if kr(i)>=k_AC & kr(i)<=k_AD
%             [X(i),Y(i)]=linecross(kr(i),jr(i),kd,jd);%求反射光线与探测面交点 
%             plot([Px(i),X(i)],[Py(i),Y(i)],'Linewidth',1);%画出反射光线
%         else
%             xx=Px(i)*3;
%             yy=kr(i)*xx+jr(i);
%             plot([Px(i),xx],[Py(i),yy],'Linewidth',1);%反射光线不能被探测面探测到
%      end
% end
 
for i=1:nn
    k1=(Py(i)-y0)/(Px(i)-x0);%求入射光线的斜率
    if k1>0
       kr(i)=-tan(alpha(i)-atan(km));%反射光线方程斜率
       jr(i)=Py(i)-kr(i)*Px(i);%反射光线截距
       if kr(i)>=-1/km & kr(i)<=k_AC
%          if kr(i)>=k_AC & kr(i)<=k_AD
            [X(i),Y(i)]=linecross(kr(i),jr(i),kd,jd);%求反射光线与探测面交点 
            plot([Px(i),X(i)],[Py(i),Y(i)],'Linewidth',1);%画出反射光线
        else
            xx=Px(i)*3;
            yy=kr(i)*xx+jr(i);
            plot([Px(i),xx],[Py(i),yy],'Linewidth',1);%反射光线不能被探测面探测到
        end
    else
        if km<0
            kr(i)=tan(atan(km)+alpha(i));%反射光线方程斜率
        else
            if alpha(i)>pi/2
                kr(i)=-tan(pi-alpha(i)-atan(km));%反射光线方程斜率
            else
                kr(i)=tan(atan(km)+alpha(i));%反射光线方程斜率
            end
        end
       jr(i)=Py(i)-kr(i)*Px(i);%反射光线截距
       if kr(i)>=k_AC & kr(i)<=k_AD
             [X(i),Y(i)]=linecross(kr(i),jr(i),kd,jd);%求反射光线与探测面交点 
             plot([Px(i),X(i)],[Py(i),Y(i)],'Linewidth',1);%画出反射光线
       else
           if alpha(i)>pi/2
               xx=-Px(i);
               yy=kr(i)*xx+jr(i);
               plot([Px(i),xx],[Py(i),yy],'Linewidth',1);%反射光线不能被探测面探测到    
           else
               xx=Px(i)*1.002;
             yy=kr(i)*xx+jr(i);
             plot([Px(i),xx],[Py(i),yy],'Linewidth',1);%反射光线不能被探测面探测到
           end
       end          
    end
end
hold off;
axis equal
toc
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


