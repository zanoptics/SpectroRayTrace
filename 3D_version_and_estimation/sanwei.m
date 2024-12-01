%2023/03/01
%点p,线L,面m
clc;format compact;clear;close all
%光源source
xs=0.5;ys=0.5;zs=0;
figure;plot3(xs,ys,zs,'o');hold on;grid on;
%镜面1 mirror1
%范围
xspan_m1=linspace(1,4,21);
yspan_m1=xspan_m1;
A_m1=-1;B_m1=-1;C_m1=5;%法向量
x_m1=0;y_m1=0;z_m1=5;%任取一点
D_m1=-(A_m1.*x_m1+B_m1.*y_m1+C_m1.*z_m1)
[X_m1,Y_m1]=meshgrid(xspan_m1,yspan_m1);
Z_m1=-(A_m1.*X_m1+B_m1.*Y_m1+D_m1)./C_m1;
surf(X_m1,Y_m1,Z_m1);
%镜面上取一光线入射点,它确定了法线的位置所以用n表示
xn=2.5;yn=2.5;zn=6;
plot3([xs,xn],[ys,yn],[zs,zn])
[xr,yr,zr]=ptline_sym(xs,ys,zs,A_m1,B_m1,C_m1,xn,yn,zn)
plot3([xn,xr],[yn,yr],[zn,zr])
%某点关于线的对称点xr,yr,zr(reflect)
%input: 点坐标xi,yi,zi(i表示入射光线上一点)，直线方向Al,Bl,Cl和其上一点xl,yl,zl
function [xr,yr,zr]=ptline_sym(xi,yi,zi,Al,Bl,Cl,xl,yl,zl)
    [xc,yc,zc]=lmcross(Al,Bl,Cl,xl,yl,zl,Al,Bl,Cl,-(Al.*xi+Bl.*yi+Cl.*zi));
    xr=2.*xc-xi;
    yr=2.*yc-yi;
    zr=2.*zc-zi;
    plot3([xl,xc],[yl,yc],[zl,zc])%对称轴
end
%任意直线line与任意平面mirror的交点xc,yc,zc(cross)
%input: 直线方向Al,Bl,Cl和其上一点xl,yl,zl，平面方程4个参量Am,Bm,Cm,Dm
function [xc,yc,zc]=lmcross(Al,Bl,Cl,xl,yl,zl,Am,Bm,Cm,Dm)
    t=-(Am.*xl+Bm.*yl+Cm.*zl+Dm)./(Al.*Am+Bl.*Bm+Cl.*Cm);
    xc=Al.*t+xl;
    yc=Bl.*t+yl;
    zc=Cl.*t+zl;
end