%% 反射定律可视化
tic %开始计时
x0=0;y0=90; %光源
Ax=25;Ay=10; %平面镜左端
Bx=60;By=10; %平面镜右端
km=(By-Ay)/(Bx-Ax); %镜斜率
jm=Ay-km*Ax; %镜截距

Cx=90;Cy=80; %探测器右
Dx=50;Dy=70; %探测器左
kn=(Cy-Dy)/(Cx-Dx);
jn=Cy-kn*Cx; 

%画
hold on; %画在同一张
daspect([1,1,1]); %设置比例尺
%两点式
plot([Ax,Bx],[Ay,By],'k','Linewidth',2);%平面镜1
plot([Cx,Dx],[Cy,Dy],'k','Linewidth',2);%探测
Nn=2; %入射光线数
z=Ax:1:Bx; %设置离散点
%在平面镜上随机取Nn个点
Px=z(randperm(numel(z),Nn)); 
Py=km*Px+jm;

[x_,y_]=symmetry(x0,y0,km,jm); %光源的虚像
k_max=(Dy-y_)/(Dx-x_); %超过此斜率必打到外面去
k_min=(Cy-y_)/(Cx-x_); %小于此斜率打不中cd
coArray=['#0072BD';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#4DBEEE';'#A2142F'];
for i=1:Nn
    color=coArray(mod(i,7)+1,:);
    plot([x0,Px(i)],[y0,Py(i)],'Color',color,'Linewidth',1.5);%入射光线
    [x2,y2]=symmetry(x0,y0,-1/km,Py(i)+Px(i)/km);%法线的截距和斜率
    kr=(y2-Py(i))/(x2-Px(i)); %反射光斜率
    jr=y2-kr*x2; %反射光截距
    [x3,y3]=linecross(kr,jr,kn,jn); %求反射光与cd的交点
    if kr<=k_max && kr>=k_min
        plot([Px(i),x3],[Py(i),y3],'Color',color,'Linewidth',1.5); %打中
    else
        xx=1.03*x3; %取一个稍远一点的点表示没打中
        yy=kr*xx+jr;
        plot([Px(i),xx],[Py(i),yy],'Color',color,'Linewidth',1.5);
    end
end
hold off;
toc
%% 求已知点关于已知直线的对称点
%x1,y1 已知点
% k b 对称轴斜率、截距
function [x2,y2]=symmetry2(x1,y1,k,b)
    x2=[];
    y2=[];
    x2=((k-1/k)*x1-2*y1+2*b)/(-k-1/k);
    y2=-(x2-x1)/k+y1;
end
%推导: 1.中点在对称轴上 2. 连线与对称轴垂直
%% 求已知点关于已知直线的对称点
%x1,y1 已知点
% k b 对称轴斜率、截距
function [x2,y2]=symmetry(x1,y1,k,b)
    [x_c,y_c]=linecross(k,b,-1/k,y1+x1/k);
    x2=2*x_c-x1;
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