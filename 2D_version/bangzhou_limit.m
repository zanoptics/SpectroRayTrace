x0=0;y0=0;
%在画出的图上点出端点,然后将游标数据导出到工作区
[A,B]=cursor_info.Position; 
x1=A(1);y1=A(2);
x2=B(1);y2=B(2);
phi=banghzou(x0,y0,x1,y1,x2,y2)
zhijin=((x1-x2)^2+(y1-y2)^2)^0.5
%%
function phi=banghzou(x0,y0,x1,y1,x2,y2)
    %中点c
    xc=(x1+x2)./2;
    yc=(y1+y2)./2;
    if abs(x1-x2)<=1e-10
        d=abs(y0-yc);
        cf=abs(x0-x1);
    else
        k12=(y2-y1)./(x2-x1);
        %过0作一12的平行线,交12的中垂线于f
        [xf,yf]=linecross(k12,y0-k12.*x0,-1./k12,yc+xc/k12);
        d=((xf-x0)^2+(yf-y0)^2)^0.5; %0到f距离d
        cf=((xf-xc)^2+(yf-yc)^2)^0.5;
    end
    phi=atand(d/cf);
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