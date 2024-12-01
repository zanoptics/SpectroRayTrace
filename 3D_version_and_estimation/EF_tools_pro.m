%% 对称式ct光路一键调参
clc;format compact
%参数
R1=304.8; %半径
f1=R1/2;
d=1/1200;
delta=10; %球面镜中心偏离90度的角度 记得去改球面镜姿态
PVs=[0,0,70]; %光源位置
b=R1/(2*cosd(delta));
PVn1=[0,R1*cosd(delta)-f1,70]  %球面镜1顶点位置
R1*(1-cosd(delta))
%球心位置计算
[xr,yr,zr]=ball2xyz(R1,90,90+delta);
PVb1=PVn1-[xr,yr,zr]

%确定光栅中心位置
PVg1=[PVb1(1),PVb1(2)+b,PVb1(3)] %xz坐标同球心,y坐标球心y+b
%屏幕位置
PVd1=[2*PVb1(1)-PVn1(1),PVs(2),PVs(3)]

%%
%1.给定两点坐标，算出两点确定向量的球坐标表示
function [r,theta,phi]=xyz2ball(PV1,PV2)
    V=PV2-PV1;
    r=dot(V,V,2).^0.5;
    theta=acosd(V(:,3)./r);
    phi=atan2d(V(:,2),V(:,1));
    phi(phi<0)=phi(phi<0)+360; %使phi的取值在0到2pi(有时不影响)
end
%2.
function [x,y,z]=ball2xyz(r,theta,phi)
    x=r.*sind(theta).*cosd(phi);
    y=r.*sind(theta).*sind(phi);
    z=r.*cosd(theta);
end