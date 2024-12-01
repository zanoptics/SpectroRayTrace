%2023/03/08
%点p,线L,面m
%PV position vector 位矢
%DV direction vector 方向矢量
clc;format compact;clear;close all
%光源source
PVs=[0.5,0.5,0];
figure;plot3(PVs(1),PVs(2),PVs(3),'o');hold on;grid on
%镜面1 mirror1
%范围
xspan_m1=linspace(1,4,21);
yspan_m1=xspan_m1;
DVm1=[-1,-1,5];%法向量
PVm1=[0,0,5];%位置
Dm1=-dot(DVm1,PVm1); %公式上是这样
[X_m1,Y_m1]=meshgrid(xspan_m1,yspan_m1);
Z_m1=mirrorZ(DVm1,PVm1,X_m1,Y_m1);%纵坐标
surf(X_m1,Y_m1,Z_m1);
%镜面上取一光线入射点,它确定了法线的位置所以用n表示
PVn=[2.5,2.5,6];
plot3([PVs(1),PVn(1)],[PVs(2),PVn(2)],[PVs(3),PVn(3)]);%入射光线
PVr=reflect(PVs,DVm1,PVn);%光源坐标,法线方向(平面法向量),入射点(法线与镜子交点)
plot3([PVr(1),PVn(1)],[PVr(2),PVn(2)],[PVr(3),PVn(3)]);%反射光线

%函数
%1.用于计算平面在XmYm范围的纵坐标Zm
%input: DVm面的法向量，PVm面上一点坐标，XmYm绘制范围
function Zm=mirrorZ(DVm,PVm,Xm,Ym)
    A=DVm(1);B=DVm(2);C=DVm(3);D=-dot(DVm,PVm);
    Zm=-(A.*Xm+B.*Ym+D)./C;
end
%2.某点关于线的对称点PVr(reflect)
%input: 点坐标PVi(i表示入射光线上一点)，直线方向DVl和其上一点PVl
function PVr=reflect(PVi,DVl,PVl)
    DVm=DVl;%法线与辅助平面垂直,PVi是辅助平面上一点
    PVc=lmcross(DVl,PVl,DVm,PVi);%线方向、位置，面方向、位置
    PVr=2.*PVc-PVi; %c是i与r的中点
    plot3([PVc(1),PVl(1)],[PVc(2),PVl(2)],[PVc(3),PVl(3)],'--');%法线
end
%3.任意直线line与任意平面mirror的交点PVc(cross)
%input: 直线方向DVl和其上一点PVl，平面方向DVm,上面一点的位矢PVm
function PVc=lmcross(DVl,PVl,DVm,PVm)
    %几何意义:PVl到面的垂直分量是DVl方向向量的t倍
    t=dot(DVm,PVm-PVl)./dot(DVm,DVl);
    PVc=t.*DVl+PVl; %直线参数方程
end