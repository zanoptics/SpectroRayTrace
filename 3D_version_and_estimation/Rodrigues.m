close all;clc;clear
%光线选项1：六角形光斑
%画球
PVb1=[1,1,0];
R1=0.5;
theta_i0=linspace(0,360,51);
phi_i0=linspace(0,360,51);
[PHI_i0,THETA_i0]=meshgrid(phi_i0,theta_i0);
X=R1.*sind(THETA_i0).*cosd(PHI_i0)+PVb1(1);
Y=R1.*sind(THETA_i0).*sind(PHI_i0)+PVb1(2);
Z=R1.*cosd(THETA_i0)+PVb1(3);
set(gcf,'position',[140 80 800 500]); %后两个是分辨率，再大好像没用
set(gca,'position',[0.07,0.1,0.9,0.85]); %后两个是图片占figure的比例
plot3(PVb1(1),PVb1(2),PVb1(3),'o');
hold on;grid on;axis equal
xlabel('x');ylabel('y');zlabel('z');
% surf(X,Y,Z)
%旋转
AxisPhi=45;
AxisTheta=90-atand((1/2)^0.5);
AxisDV=rotationmat(AxisPhi,AxisTheta)*[0,0,1]';
AxisPV=[-1,1,1];
plot3([AxisPV(1)-2*AxisDV(1),AxisPV(1)+2*AxisDV(1)],[AxisPV(2)-2*AxisDV(2),AxisPV(2)+2*AxisDV(2)],[AxisPV(3)-2*AxisDV(3),AxisPV(3)+2*AxisDV(3)],'-');
[X_rot,Y_rot,Z_rot]=axisrot(100*AxisDV,AxisPV,66.7,PVb1(1),PVb1(2),PVb1(3));
plot3(X_rot,Y_rot,Z_rot,'o')
%%
clc
%球参数
PVb1=[1,1,0];
R1=0.5;
theta_i0=linspace(0,360,41);
phi_i0=linspace(0,360,41);
[PHI_i0,THETA_i0]=meshgrid(phi_i0,theta_i0);
X=R1.*sind(THETA_i0).*cosd(PHI_i0)+PVb1(1);
Y=R1.*sind(THETA_i0).*sind(PHI_i0)+PVb1(2);
Z=R1.*cosd(THETA_i0)+PVb1(3);
%rot参数
AxisPhi=45;
AxisTheta=atand((1/2)^0.5);
AxisDV=rotationmat(AxisPhi,AxisTheta)*[0,0,1]';
AxisPV=[-1,1,1];
%动画
filename ='ball.gif ';
num=37;
m=moviein(num);
AngRot=linspace(0,360,num);
for idx=1:num
    fig=figure('Visible','off');hold on;view(134.69,54.51);daspect([1 1 1]);axis([-4 2 -1.8 3.8 -1 3]); 
    plot3([AxisPV(1)-2*AxisDV(1),AxisPV(1)+2*AxisDV(1)],[AxisPV(2)-2*AxisDV(2),AxisPV(2)+2*AxisDV(2)],[AxisPV(3)-2*AxisDV(3),AxisPV(3)+2*AxisDV(3)],'-');
    [X_rot,Y_rot,Z_rot]=axisrot(AxisDV,AxisPV,AngRot(idx),X,Y,Z);
    surf(X_rot,Y_rot,Z_rot)
    [X_rot,Y_rot,Z_rot]=axisrot(AxisDV,AxisPV,-AngRot(idx),X,Y,Z);
    surf(X_rot,Y_rot,Z_rot)
    m(:,idx)=getframe(fig);
    img =frame2im(m(:, idx ));
    [ A , cmap ]=rgb2ind( img ,256);
    if idx ==1
    imwrite( A , cmap , filename ,'gif','loopcount', Inf ,'DelayTime',0.1);
    else
    imwrite( A , cmap , filename ,'gif', 'WriteMode','append','DelayTime',0.1);
    end 
end
fig.Visible='on';
movie(m,1)
hold off;

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
%14.绕轴旋转
%参考：罗德里格旋转公式（Rodrigues' rotation formula）
%DV旋转轴方向向量，PV旋转轴上一点的坐标，theta旋转角度，后三个是被旋转物体的xyz
function [obj1x,obj1y,obj1z]=axisrot(DV,PV,theta,obj1x,obj1y,obj1z)
    DV=DV./dot(DV,DV)^0.5;%转为单位向量
    Rk=[0,-DV(3),DV(2);
        DV(3),0,-DV(1);
        -DV(2),DV(1),0]; %叉乘矩阵
    M=eye(3)+(1-cosd(theta)).*(Rk*Rk)+sind(theta).*Rk; %绕轴转矩阵
    obj1=[obj1x(:)'-PV(1);obj1y(:)'-PV(2);obj1z(:)'-PV(3)];
    temp=M*obj1;
    obj1x(:)=temp(1,:)+PV(1);
    obj1y(:)=temp(2,:)+PV(2);
    obj1z(:)=temp(3,:)+PV(3);
end