clc;format compact;clear;close all
theta_i=5;phi_i=linspace(0,360,7);
[PHI_i,THETA_i]=meshgrid(phi_i,theta_i);
sx_i=sind(THETA_i).*cosd(PHI_i);
sy_i=sind(THETA_i).*sind(PHI_i);
sz_i=cosd(THETA_i);
DVi_set=(rotationmat(90,90)*[sx_i;sy_i;sz_i])';%方向
for jj=1:size(DVi_set,1)
    DVi1=DVi_set(jj,:)
end
% figure;plot3(DVi_set(:,1),DVi_set(:,2),DVi_set(:,3))






%1.旋转矩阵函数
function Rmat=rotationmat(zphi,ytheta)
    ythetamat=[
        cosd(ytheta),0,sind(ytheta);
        0,1,0;
        -sind(ytheta),0,cosd(ytheta);];%旋转矩阵Ry,先绕y轴转θ
    zphimat=[
        cosd(zphi),-sind(zphi),0;
        sind(zphi),cosd(zphi),0;
        0,0,1;];%旋转矩阵Rz,再绕z轴转φ
    Rmat=zphimat*ythetamat;%先把两个旋转变换乘起来
end
%8.计算任意位置任意朝向的球面镜坐标
%input:PVb球心,R半径,取值范围phi,theta,R与theta共同决定口径
%球面镜朝向角度zphi,ytheta
function [sx,sy,sz]=ball_transform(PVb,R,phi,theta,zphi,ytheta)
    [PHI,THETA]=meshgrid(phi,theta);
    %旋转对象(一个球心在原点的顶朝上的球面镜)
    sx=R.*sind(THETA).*cosd(PHI);
    sy=R.*sind(THETA).*sind(PHI);
    sz=R.*cosd(THETA);
%     surf(sx,sy,sz, ...
%     'FaceColor','interp','EdgeColor','none');
    %旋转到想要的方向
    Rmat=rotationmat(zphi,ytheta);
    %将画图时用的三个矩阵合并,每一列的三行对应一个点的xyz坐标
    object=[sx(1:end);sy(1:end);sz(1:end)];
    %旋转变换
    temp=Rmat*object;
    %更新坐标,同时平移
    sx(1:end)=temp(1,:)+PVb(1);
    sy(1:end)=temp(2,:)+PVb(2);
    sz(1:end)=temp(3,:)+PVb(3);
end