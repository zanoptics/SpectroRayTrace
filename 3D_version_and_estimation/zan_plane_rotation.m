clc;format compact;clear;close all
PVi=[0,0,0];%线上一点(光源)
figure;plot3(PVi(1),PVi(2),PVi(3),'o');
hold on;grid on;axis equal
xlabel('x');ylabel('y');zlabel('z');
lx=8;ly=4;
zphi_p=90;ytheta_p=90;%先绕y轴转theta,再绕z轴转phi
PVm=[5,5,5];plot3(PVm(1),PVm(2),PVm(3),'.');
[px,py,pz,DVm]=plane_transform(PVm,lx,ly,zphi_p,ytheta_p);
surf(px,py,pz, ...
    'FaceColor','interp','EdgeColor','none');
%球
PVb=[2,2,2];%球心,也是平移向量
plot3(PVb(1),PVb(2),PVb(3),'o');
R=1;
%先画一个球心在原点的顶朝上的球面镜,再旋转与平移
phi=linspace(0,360,31);%与x轴夹角范围
theta=linspace(0,45,31);%与z轴夹角范围,与口径相关
zphi_b=90;ytheta_b=45;%旋转角度
%计算球面镜坐标
[sx,sy,sz]=ball_transform(PVb,R,phi,theta,zphi_b,ytheta_b);%旋转+平移
%画球面镜
surf(sx,sy,sz, ...
    'FaceColor','interp','EdgeColor','none');
%函数
%-1.旋转矩阵函数
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
%0.画任意位置任意朝向的平面plane_transform
%input:平面中心坐标PVm,x方向长度lx,y方向长度ly,旋转角度zph,ytheta
%返回计算出来的坐标px,py,pz,和法向量用于计算直线与此平面的交点
function [px,py,pz,DVm]=plane_transform(PVm,lx,ly,zphi,ytheta)
    %思路是先画个躺在xoy平面的矩形,然后旋转再平移,与画球面镜做法相似
    [px,py]=meshgrid(linspace(-lx./2,lx./2,5),linspace(-ly./2,ly./2,5));
    pz=zeros(size(px));%纵坐标为0
    Rmat=rotationmat(zphi,ytheta);
    object=[px(1:end);py(1:end);pz(1:end)];
    temp=Rmat*object;
    px(1:end)=temp(1,:)+PVm(1);
    py(1:end)=temp(2,:)+PVm(2);
    pz(1:end)=temp(3,:)+PVm(3);
    DVm=(Rmat*[0;0;1])';%对z轴单位向量作变换
end
%5.计算任意位置任意朝向的球面镜坐标
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
%6.判断交点在不在镜子上
function flag=onBallMirror(PVn,PVb,R,theta_max,zphi,ytheta)
    %旋转矩阵是正交矩阵，它的逆就等于转置，我们把给交点平移，再左乘转置
    %就得到其在初始顶朝上球面镜上的对应点
    Rmat=rotationmat(zphi,ytheta);%旋转矩阵
    PVn_ori=(Rmat')*(PVn-PVb)';
%     plot3(PVn_ori(1),PVn_ori(2),PVn_ori(3),'.');
    if PVn_ori(3)>=R*cosd(theta_max)
        flag=1;
    else
        flag=0;
    end
end
%7.判断交点PVn在不在平面镜上
function flag=onPlaneMirror(PVn,PVm,lx,ly,zphi,ytheta)
    Rmat=rotationmat(zphi,ytheta);%旋转矩阵
    PVn_ori=(Rmat')*(PVn-PVm)';
    if abs(PVn_ori(1))<=lx/2 && abs(PVn_ori(2))<=ly/2
        flag=1;
    else
        flag=0;
    end
end