clc;format compact;clear;close all
PVs=[0.5,0.5,0];
PVn=[2.5,2.5,6];
PVo=[6.42592592592593,6.42592592592593,2.37037037037037];
PVc=[3.46296296296296,3.46296296296296,1.18518518518519];
A=PVs-PVn;
B=PVo-PVn;
C=PVc-PVn;
lhs=dot(A,C)
rhs=dot(B,C)
% lhs=sum(A.*A)
% rhs=sum(B.*B)
%%
%2023/03/01
clc;format compact;clear;close all
x0=0.5;y0=0.5;z0=0;
figure;plot3(x0,y0,z0,'o');hold on;grid on;
xspan_m1=linspace(1,4,21);
yspan_m1=xspan_m1;
A_m1=-1;B_m1=-1;C_m1=5;
xa_m1=0;ya_m1=0;za_m1=5;
D_m1=-(A_m1.*xa_m1+B_m1.*ya_m1+C_m1.*za_m1)
[X_m1,Y_m1]=meshgrid(xspan_m1,yspan_m1);
Z_m1=-(A_m1.*X_m1+B_m1.*Y_m1+D_m1)./C_m1;
surf(X_m1,Y_m1,Z_m1);
xs=2.5;ys=2.5;zs=6;
plot3([x0,xs],[y0,ys],[z0,zs])
[x3,y3,z3]=ptline_sym(x0,y0,z0,A_m1,B_m1,C_m1,xs,ys,zs)
plot3([xs,x3],[ys,y3],[zs,z3])
%某点关于线的中心对称点
%input: 点坐标x0,y0,z0，直线方向A,B,C和其上一点x2,y2,z2
function [x3,y3,z3]=ptline_sym(x0,y0,z0,A,B,C,x2,y2,z2)
    [xc,yc,zc]=lpcross(A,B,C,x2,y2,z2,A,B,C,-(A.*x0+B.*y0+C.*z0));
    
    x3=2.*xc-x0;
    y3=2.*yc-y0;
    z3=2.*zc-z0;
    plot3([x2,xc],[y2,yc],[z2,zc])

end
%任意直线line与任意平面plane的交点xc,yc,zc
%input: 直线方向A1,B1,C1和其上一点x1,y1,z1，平面方程4个参量A2,B2,C2,D2
function [xc,yc,zc]=lpcross(A1,B1,C1,x1,y1,z1,A2,B2,C2,D2)
    t=-(A2.*x1+B2.*y1+C2.*z1+D2)./(A1.*A2+B1.*B2+C1.*C2);
    xc=A1.*t+x1;
    yc=B1.*t+y1;
    zc=C1.*t+z1;
end