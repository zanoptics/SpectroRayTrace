%2023/03/08
%球-尝试
clc;format compact;clear;close all
PVi=[0 0 0];
figure;plot3(PVi(1),PVi(2),PVi(3),'o');hold on;grid on;axis equal
xlabel('x');ylabel('y');zlabel('z');
xlim([-2,3]);ylim([-2,3]);zlim([-1,3]);
%球
PVb=[1,1,1];
plot3(PVb(1),PVb(2),PVb(3),'o');
R=1;
phi=linspace(0,360,31);%与x轴夹角
theta=linspace(0,30,31);%与z轴夹角
[PHI,THETA]=meshgrid(phi,theta);
sx=R.*sind(THETA).*cosd(PHI);
sy=R.*sind(THETA).*sind(PHI);
sz= R.*cosd(THETA);
zphi=90;ytheta=45;
zphimat=[
    cosd(zphi),-sind(zphi),0;
    sind(zphi),cosd(zphi),0;
    0,0,1;];
ythetamat=[
    cosd(ytheta),0,sind(ytheta);
    0,1,0;
    -sind(ytheta),0,cosd(ytheta);];
rotationmat=zphimat*ythetamat
object=[sx(1:end);sy(1:end);sz(1:end)];
temp=rotationmat*object;
sx(1:end)=temp(1,:);
sy(1:end)=temp(2,:);
sz(1:end)=temp(3,:);
surf(sx,sy,sz, ...
    'FaceColor','b','EdgeColor','none');
surf(sx+PVb(1),sy+PVb(2),sz+PVb(3), ...
    'FaceColor','y','EdgeColor','none');