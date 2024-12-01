clear;
%比例尺
daspect([1 1 1]);
hold on;
%参数范围
t=-1.1:0.001:-0.3;
%焦点
x0=6;y0=8;
focus=[x0,y0];

%准线斜率截距(普通方程需要的量)
k=-2;b=5;
x=-0.5:1:2.5;
M=k.*x+b;
plot(x,M,'r--','LineWidth',1);
%换算为 参数方程 需要的量
p=b/sqrt(k^2+1); %(0,0)代入距离公式
phi=atand(k)-90; %旋转角度 从0°到-180°
%焦点在原点的抛物线的参数方程,普通方程为 y方=2p(x+p/2)
x1=2*p.*t.^2-p/2;
y1=2*p.*t;
plot(x1,y1,'k-','LineWidth',2);%抛物线
plot([-p,-p],[-5,5],'k--','LineWidth',1); %准线x=-p
scatter(0,0,'k*');%焦点
%对每个点,绕原点旋转变换
x2=zeros(1,numel(x1));y2=zeros(1,numel(y1));%初始化
for i=1:numel(x1)
    temp=num2cell([cosd(phi),-sind(phi);sind(phi),cosd(phi)]*[x1(i);y1(i)]);%旋转矩阵
    [x2(i),y2(i)]=deal(temp{:});%好像必须这样赋值(先变为cell数组，然后deal),有好方法告诉我
end
plot(x2,y2,'r-','LineWidth',2);%抛物线
%平移变换,变为以focus为焦点
plot(x2+focus(1),y2+focus(2),'-','Color',[0 0.45 0.74],'LineWidth',2)%抛物线
plot(x+focus(1),M+focus(2),'--','Color',[0 0.45 0.74],'LineWidth',1);%准线
scatter(x0,y0,'*','CData',[0 0.45 0.74]);%焦点

%与普通方程的效果对比
f3=@(x,y) ((x-x0)+k.*(y-y0)).^2-2.*k.*b.*(x-x0)+2.*b.*(y-y0)-b^2;
fimplicit(f3,'k.','LineWidth',2);
i=1;

% %尝试抛物面镜的反射
% for m1=-2:0.5:-0.2
%     i=i+1;
%     [x_c1,y_c1]=lpcross(k,b,m1,x0,y0,2);
% %     plot([x0,x_c1],[y0,y_c1],'Color','red','Linewidth',1.5);
%     kr1=pslope(x_c1,y_c1,x0,y0,k,b);
%     [x9,y9,flag]=reflect(x0,y0,x_c1,y_c1,kr1,-1,0,-5,2,i);
%     %                              让它打到y=-1*x+0上,x介于[-5,2]间
% end
hold off;
%% 求直线与抛物线的交点
%准线斜率k 截距b,入射光斜率m,焦点x0,yo,flag取哪一边的点1or2
%flag选1还是2，要试验，我们用的抛物面镜要么是抛物线对称轴左半部分的一小部分，
% 要么是右半部分的一小部分
function [x_c,y_c]=lpcross(k,b,m,x0,y0,flag)
    if flag==1
        x_c=(b*(k-m)+b*(k^2+m^2+k^2*m^2+1)^0.5)/(k*m+1)^2;
    elseif flag==2
        x_c=(b*(k-m)-b*(k^2+m^2+k^2*m^2+1)^0.5)/(k*m+1)^2;
    end
    y_c=m*x_c;
    x_c=x_c+x0;
    y_c=y_c+y0;
end
%% 抛物线上一点x,y的斜率
%抛物线焦点x0,y0移到原点后的准线的斜率k与截距b
function k=pslope(x,y,x0,y0,k,b)
    x=x-x0;y=y-y0; %先把点移到原点抛物线上
    k=(k*b-x-k*y)/(k^2*y+k*x+b);%原点抛物线上一点的切线斜率
end
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
    %km=0,法线斜率为Inf,故分类讨论
    if km==0
        x3=2*x2-x1;
        y3=y1;
    else
        [x3,y3]=symmetry(x1,y1,-1/km,y2+x2/km); %关于法线的对称点
    end
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