clc;format compact
%入射光线4方向向量的方向角，让光线垂直打在屏幕上，即屏幕姿态角
[r,theta,phi]=xyz2ball(PVs,PVn1)
%% 对称式ct光路一键调参
clc;format compact
%参数
R1=304.8; %半径
f1=R1/2;
d=1/1200;
delta=10; %球面镜中心偏离90度的角度 记得去改球面镜姿态
PVs=[0,0,70]; %光源位置
PVn1=[0,152.4,70] %球面镜1顶点位置

%球心位置计算
[xr,yr,zr]=ball2xyz(R1,90,90+delta);
PVb1=PVn1-[xr,yr,zr]
%确定光栅中心位置
b=R1/(2*cosd(delta));
PVg1=[PVb1(1),PVb1(2)+b,PVb1(3)] %xz坐标同球心,y坐标球心y+b
%屏幕位置
PVd1=[2*PVb1(1)-PVn1(1),PVs(2),PVs(3)]

%光栅旋转角度φ和要对准探测器的波长λ的关系
lambda=(6).*100e-6;%调整
zphi_g1=acosd(lambda/(-2*cosd(2*delta)*d))-90

% **可分的光波长上限**
%1.看公式sinθ=λ/d-sini 显然λ/d不得超过2
lambda_sup1=2*d*1e6 %nm

%2.能对准探测器的波长也有上限 由光路结构决定 δ=10°时最大1470左右
% φ-90<90-2δ 因为对准时该波长光线的反射角为2δ
lambda_sup2=2*d*sind(90-2*delta)*cosd(2*delta)*1e6 %nm

%% 光栅旋转角度φ和要对准探测器的波长λ的关系(画图)
clc;format compact;clear
d=1/1200;
delta=10; %球面镜中心偏离90度的角度
lambda=(8:0.001:15).*100e-6;
% lambda=9.*100e-6;
zphi_g1=acosd(lambda/(-2*cosd(2*delta)*d));
plot(lambda*1e6,zphi_g1,'-');title('波长与光栅旋转角度的关系');
cha=diff(zphi_g1);
mean(cha)
%%
clc;format compact;
d=1/1200;
delta=10; %球面镜中心偏离90度的角度
motor=0.0140625;
num=3059;
start=170;
ed=start-motor*(num-1);
ang_range=start:-motor:ed;
lambda_range=-2.*cosd(2.*delta).*d.*cosd(ang_range);
plot(ang_range,lambda_range*1e6,'-');title('波长与光栅旋转角度的关系');
%%
clc
d=1/1200;
theta_blazed=36+52/60;%闪耀角
i=36.1;
lambda_blazed=2.*d.*sind(theta_blazed).*cosd(theta_blazed-i)
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