%2023/03/15
%球-尝试
clc;format compact;clear;close all
%入射光线
PVi=[0.5,0.5,0];%线上一点(光源)
DVi=[1,1.3,2];%随便取了一个方向
figure;plot3(PVi(1),PVi(2),PVi(3),'o');
hold on;grid on;axis equal
xlabel('x');ylabel('y');zlabel('z');
%球
PVb=[2,2,2];%球心,也是平移向量
plot3(PVb(1),PVb(2),PVb(3),'o');
R=1;
%先画一个球心在原点的顶朝上的球面镜,再旋转与平移
phi=linspace(0,360,31);%与x轴夹角范围
theta=linspace(0,45,31);%与z轴夹角范围,与口径相关
zphi=90;ytheta=45;%旋转角度
%计算球面镜坐标
[sx,sy,sz]=ball_transform(PVb,R,phi,theta,zphi,ytheta);%旋转+平移
%画球面镜
surf(sx,sy,sz, ...
    'FaceColor','interp','EdgeColor','none');
PVn=lbcross(DVi,PVi,PVb,R);%光线与球的交点
%判断是不是光线与镜子的交点,1在，0不在
flagb1=onBallMirror(PVn,PVb,R,theta(end),zphi,ytheta)
plot3(PVn(1),PVn(2),PVn(3),'r.');
DVn=balln(PVn,PVb);%法线
plot3([PVi(1),PVn(1)],[PVi(2),PVn(2)],[PVi(3),PVn(3)]);%入射光线
PVr=reflect(PVi,DVn,PVn);%套平面镜的反射
plot3([PVr(1),PVn(1)],[PVr(2),PVn(2)],[PVr(3),PVn(3)]);%反射光线

%函数
%0.旋转矩阵函数
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
%1.求直线与球的交点PVn
%input: 直线方向DVl,线上一点PVl,球心PVb,半径R
function PVn=lbcross(DVl,PVl,PVb,R)
    DVd=PVl-PVb;
    a=dot(DVl,DVl);
    b=2.*dot(DVd,DVl);
    c=dot(DVd,DVd)-R.^2;
%     %数值求解at^2+bt+c=0
%     syms x
%     S=double(vpasolve(a.*x.^2+b.*x+c == 0, x));
    S=zeros(2,1);
    delta=b.^2-4.*a.*c;
    if delta<0
        disp('与球无交点');
    end
    S(1)=(-b+delta.^0.5)./(2.*a);
    S(2)=(-b-delta.^0.5)./(2.*a);
    %取绝对值大的s
    if abs(S(1))>abs(S(2))
        t=S(1);
    else
        t=S(2);
    end
    PVn=t.*DVl+PVl;
end
%2.球上一点的法向量
%input: 球上一点,球心PVb
function DVn=balln(PVn,PVb)
    DVn=PVb-PVn;
end
%平面镜反射的函数
%3.某点关于线的对称点PVr(reflect)
%input: 点坐标PVi(i表示入射光线上一点)，直线方向DVl和其上一点PVl
function PVr=reflect(PVi,DVl,PVl)
    DVm=DVl;%法线与辅助平面垂直,PVi是辅助平面上一点
    PVc=lmcross(DVl,PVl,DVm,PVi);%线方向、位置，面方向、位置
    PVr=2.*PVc-PVi; %c是i与r的中点
    plot3([PVc(1),PVl(1)],[PVc(2),PVl(2)],[PVc(3),PVl(3)],'--');%法线
end
%4.任意直线line与任意平面mirror的交点PVc(cross)
%input: 直线方向DVl和其上一点PVl，平面方向DVm,上面一点的位矢PVm
function PVc=lmcross(DVl,PVl,DVm,PVm)
    %几何意义:PVl到面的垂直分量是DVl方向向量的t倍
    t=dot(DVm,PVm-PVl)./dot(DVm,DVl);
    PVc=t.*DVl+PVl; %直线参数方程
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
