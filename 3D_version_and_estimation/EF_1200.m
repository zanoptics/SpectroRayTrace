%2023/05/04
%两个球面镜CT光路但是等同于EF的光路
%规定：光线方向向量与光线传播方向一致
%长度单位均为mm
clc;format compact;clear;close all
%仪器布置
%1.光源
PVs=[0,0,70];
figure;
set(gcf,'position',[140 80 800 500]); %后两个是分辨率，再大好像没用
set(gca,'position',[0.07,0.1,0.9,0.85]); %后两个是图片占figure的比例
plot3(PVs(1),PVs(2),PVs(3),'o');
hold on;grid on;axis equal
xlabel('x');ylabel('y');zlabel('z');
%2.球面镜1(准直)
PVb1=[52.9280 -147.7694   70.0000];
% plot3(52.9280,-147.7694,70.0000,'ro');
R1=304.8;%球心与半径
aperture1=25.4;theta_b1M=asind(aperture1./(2.*R1));%口径与张角
%先画一个球心在原点的顶朝上球面镜
phi_b1=linspace(0,360,31);
theta_b1=linspace(0,theta_b1M,31);
zphi_b1=100;ytheta_b1=90; %姿态参数
%旋转+平移
[bx1,by1,bz1]=ball_transform(PVb1,R1,phi_b1,theta_b1,zphi_b1,ytheta_b1);
surf(bx1,by1,bz1, ...
    'FaceColor','interp','EdgeColor','none');
%3.光栅
PVg1=[52.9280    6.9816   70.0000];%中心位置
PVg1_axis=[52.9280    6.9816   40.0000];%转轴底部
lx_g1=12.7;ly_g1=12.7;%平铺时的长宽尺寸
zphi_g1=90;ytheta_g1=90; %姿态参数
ytheta_g1_delta=0;
d=1/1200;%光栅常数,300条/mm
theta_blazed=36+52/60;%闪耀角


%误差分析v%
%旋转+平移
[px1,py1,pz1,DVg1]=plane_transform(PVg1,lx_g1,ly_g1,zphi_g1,ytheta_g1);%正常的情况

%使得光栅微微倾斜
[px1,py1,pz1,DVg1,PVg1]=axisrot([1,0,0],PVg1_axis,ytheta_g1_delta,px1,py1,pz1,DVg1,PVg1);

%求出倾斜了的电机旋转轴
motor_DV=cross([1,0,0],DVg1);
%光栅旋转角度
AngRot=52.9522;
[px1,py1,pz1,DVg1,PVg1]=axisrot(motor_DV,PVg1_axis,AngRot,px1,py1,pz1,DVg1,PVg1);
%误差分析^%

surf(px1,py1,pz1, ...
    'FaceColor','interp','EdgeColor','none');
%4.球面镜2(会聚)
PVb2=PVb1;R2=304.8;%球心与半径
aperture2=25.4;theta_b2M=asind(aperture2./(2.*R2));%口径与张角
%先画一个球心在原点的顶朝上球面镜
phi_b2=linspace(0,360,31);
theta_b2=linspace(0,theta_b2M,31);
zphi_b2=80;ytheta_b2=90; %姿态参数
%旋转+平移
[bx2,by2,bz2]=ball_transform(PVb2,R2,phi_b2,theta_b2,zphi_b2,ytheta_b2);
surf(bx2,by2,bz2, ...
    'FaceColor','interp','EdgeColor','none');
%5.ccd或狭缝或观察屏detector
PVd1=[105.8559         0   70.0000];%中心位置
plot3(PVd1(1),PVd1(2),PVd1(3),'wo','LineWidth',1);
lx_d1=25;ly_d1=25;%平铺时的长宽尺寸
zphi_d1=90;ytheta_d1=90; %姿态参数(根据入射光线4+xyz2ball函数确定)
%旋转+平移
[px2,py2,pz2,DVd1]=plane_transform(PVd1,lx_d1,ly_d1,zphi_d1,ytheta_d1);
surf(px2,py2,pz2, ...
    'FaceColor','interp','EdgeColor','none');
%挡光板位置
plot3([52.928,52.928],[70,150],[70,70])
% plot3([90,90],[0,30],[70,70])
%----------------------------------------------
%光线追迹
% 1.入射光线1
%光线选项1：六角形光斑
theta_i0=2.2;%发散角
phi_i0=linspace(0,360,9);phi_i0=phi_i0(1:end-1);%一圈上等角度间隔取点
[PHI_i0,THETA_i0]=meshgrid(phi_i0,theta_i0);
sx_i=sind(THETA_i0).*cosd(PHI_i0);
sy_i=sind(THETA_i0).*sind(PHI_i0);
sz_i=cosd(THETA_i0);
DVi_set=(rotationmat(90,90)*[sx_i;sy_i;sz_i])';%旋转
% %光线选项2：矩形光斑(包括狭缝)
% d_i=15;l_i=15;%缝宽和缝长
% x_i=linspace(-d_i/2,d_i/2,3);
% y_i=linspace(-l_i/2,l_i/2,3);
% [X_i,Y_i]=meshgrid(x_i,y_i); Z_i=zeros(size(Y_i))+150;
% DVi_set=(rotationmat(90,90)*[X_i(:)';Y_i(:)';Z_i(:)'])';%旋转
%复色光颜色参数
% lambda=[400e-6 500e-6 600e-6 700e-6 1000e-6 1100e-6 1200e-6 1300e-6 1400e-6]';%列向量
% lambda=([14.7,14.02,14.04])'.*100e-6;
% lambda=(8:0.125:9)'.*100e-6;
% lambda=[1200e-6 1201e-6]';
lambda=(12:0.125:13.5)'.*100e-6;
lambda_blazed_set=zeros(size(lambda));
for jj=1:size(DVi_set,1)
    DVi1=DVi_set(jj,:);
    %2.入射光线1与球1交点
    PVn1=lbcross(DVi1,PVs,PVb1,R1);%交点
    %判断是否在球面镜上
    flag_n1=onBallMirror(PVn1,PVb1,R1,theta_b1M,zphi_b1,ytheta_b1);
    plot3([PVs(1),PVn1(1)],[PVs(2),PVn1(2)],[PVs(3),PVn1(3)], ...
        'Color','k','LineWidth',1,'Marker','.');%入射光线1
    DVn1=balln(PVn1,PVb1);%法向量
    %3.反射光线1到入射光线2, DVr1-->DVi2,PVn1-->PVi2
    PVr1=reflect(PVs,DVn1,PVn1);%反射光线上一点
    DVi2=PVr1-PVn1;% DVr1(略)-->DVi2
    PVi2=PVn1;
    %4.入射光线2与光栅的交点
    PVn2=lmcross(DVi2,PVi2,DVg1,PVg1);
%     flag_n2=onPlaneMirror(PVn2,PVg1,lx_g1,ly_g1,zphi_g1,ytheta_g1);
    plot3([PVi2(1),PVn2(1)],[PVi2(2),PVn2(2)],[PVi2(3),PVn2(3)], ...
        'Color','k','LineWidth',1,'Marker','.');%入射光线2
    %5.反射光线2到入射光线3,DVr2-->DVi3,PVn2-->PVi3
    for m=[1] %干涉级次
        [PVr2,theta_i2]=grating_reflect(m,d,lambda,PVi2,DVg1,PVn2);%矩阵每行是不同波长对应的反射光线上一点坐标
        lambda_blazed_set(jj,1)=grating_blaze(m,d,theta_blazed,theta_i2);%计算闪耀波长
        DVi3=PVr2-PVn2;% DVr2(略)-->DVi3
        PVi3=PVn2;
        %6.入射光线3与球3交点
        PVn3=zeros(size(DVi3));DVn3=PVn3;PVr3=PVn3;
        DVi4=PVn3;PVn4=PVn3;
        for ii=1:length(lambda)
            if isnan(DVi3(ii,:))
                continue
            end
            [PVn3(ii,:),flag_b2]=lbcross(DVi3(ii,:),PVi3,PVb2,R2);%暂时没考虑没有交点的情况
            flag_n3=onBallMirror(PVn3(ii,:),PVb2,R2,theta_b2M,zphi_b2,ytheta_b2);
    %         if flag_n3==0
    %             continue
    %         end
            lightcolor=getRGB(lambda(ii)*1e6);
            plot3([PVi3(1),PVn3(ii,1)],[PVi3(2),PVn3(ii,2)],[PVi3(3),PVn3(ii,3)], ...
                'Color',lightcolor,'LineWidth',1,'Marker','.');%入射光线3
            DVn3(ii,:)=balln(PVn3(ii,:),PVb2);
            %7.反射光线3到入射光线4,DVr3-->DVi4,PVn3-->PVi4
            PVr3(ii,:)=reflect(PVi3,DVn3(ii,:),PVn3(ii,:));
            DVi4(ii,:)=PVr3(ii,:)-PVn3(ii,:); %DVr3(略)-->DVi4
            PVi4=PVn3;
            %8.入射光线4与探测器的交点
            PVn4(ii,:)=lmcross(DVi4(ii,:),PVi4(ii,:),DVd1,PVd1);
            flag_n4=onPlaneMirror(PVn4(ii,:),PVd1,lx_d1,ly_d1,zphi_d1,ytheta_d1);
            plot3([PVi4(ii,1),PVn4(ii,1)],[PVi4(ii,2),PVn4(ii,2)],[PVi4(ii,3),PVn4(ii,3)], ...
                'Color',lightcolor,'LineWidth',1,'Marker','.');%入射光线4
        end
    end
end
%----------------------------------------------
%函数
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
%2画任意位置任意朝向的平面plane_transform
%input:平面中心坐标PVm,x方向长度lx,y方向长度ly,旋转角度zph,ytheta
%返回计算出来的坐标px,py,pz,和法向量用于计算直线与此平面的交点
function [px,py,pz,DVm]=plane_transform(PVm,lx,ly,zphi,ytheta)
    %思路是先画个躺在xoy平面的矩形,然后旋转再平移,与画球面镜做法相似
    [px,py]=meshgrid(linspace(-lx./2,lx./2,5),linspace(-ly./2,ly./2,5));
    pz=zeros(size(px));%纵坐标为0
    Rmat=rotationmat(zphi,ytheta);
    object=[px(:)';py(:)';pz(:)'];
    temp=Rmat*object;
    px(:)=temp(1,:)+PVm(1);
    py(:)=temp(2,:)+PVm(2);
    pz(:)=temp(3,:)+PVm(3);
    DVm=(Rmat*[0;0;1])';%对z轴单位向量作变换
end
%3判断交点PVn在不在平面镜上
function flag=onPlaneMirror(PVn,PVm,lx,ly,zphi,ytheta)
    %旋转矩阵是正交矩阵，它的逆就等于转置，我们把给交点平移，再左乘转置
    Rmat=rotationmat(zphi,ytheta);%旋转矩阵
    PVn_ori=(Rmat')*(PVn-PVm)';
    if abs(PVn_ori(1))<=lx/2 && abs(PVn_ori(2))<=ly/2
        flag=1;
    else
        flag=0;
    end
end
%4.任意直线line与任意平面mirror的交点PVc(cross)
%input: 直线方向DVl和其上一点PVl，平面方向DVm,上面一点的位矢PVm
function PVc=lmcross(DVl,PVl,DVm,PVm)
    %几何意义:PVl到面的垂直分量是DVl方向向量的t倍
    t=dot(DVm,PVm-PVl)./dot(DVm,DVl);
    PVc=t.*DVl+PVl; %直线参数方程
end
%平面镜反射的函数
%5.某点关于线的对称点PVr(reflect)
%input: 点坐标PVi(i表示入射光线上一点)，直线方向DVl和其上一点PVl
function PVr=reflect(PVi,DVl,PVl)
    %法线与辅助平面垂直(线方向和面方向相同),PVi是辅助平面上一点
    PVc=lmcross(DVl,PVl,DVl,PVi);%线方向、位置，面方向、位置
    PVr=2.*PVc-PVi; %c是i与r的中点
%     plot3([PVc(1),PVl(1)],[PVc(2),PVl(2)],[PVc(3),PVl(3)],'--');%法线
end
%6.求直线与球的交点PVn
%input: 直线方向DVl,线上一点PVl,球心PVb,半径R
function [PVn,flag]=lbcross(DVl,PVl,PVb,R)
    flag=1;
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
        flag=0;
        return
    end
    S(1)=(-b+delta.^0.5)./(2.*a);
    S(2)=(-b-delta.^0.5)./(2.*a);
    %取大于0且绝对值最小的(前提是光线方向向量是其传播方向)
    t=min(S(S>0));
    PVn=t.*DVl+PVl;
end
%7.球上一点的法向量
%input: 球上一点,球心PVb
function DVn=balln(PVn,PVb)
    DVn=PVb-PVn;
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
    object=[sx(:)';sy(:)';sz(:)'];
    %旋转变换
    temp=Rmat*object;
    %更新坐标,同时平移
    sx(:)=temp(1,:)+PVb(1);
    sy(:)=temp(2,:)+PVb(2);
    sz(:)=temp(3,:)+PVb(3);
end
%9.判断交点在不在球面镜上
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
%10.光栅反射函数grating_reflect
%input:级数m,光栅常数d,波长λ(可以是一个列向量),入射光线上一点PVi,法线方向DVl,交界点PVl
function [PVgr,theta_i]=grating_reflect(m,d,lambda,PVi,DVl,PVl)
    %法线与辅助平面垂直DVl就是mirror法向量,且PVi也是辅助平面上一点
    PVc=lmcross(DVl,PVl,DVl,PVi);%线方向、位置，面方向、位置
    %求
    theta_i=acosd(dot(PVi-PVl,PVc-PVl)./(dot(PVi-PVl,PVi-PVl).*dot(PVc-PVl,PVc-PVl)).^0.5);
    theta_r=asind(sind(theta_i)-m.*lambda./d);%列向量
    for jj=1:length(lambda)
        if imag(theta_r(jj))~=0
            disp([num2str(m),'级衍射中',num2str(lambda(jj)*1e6),'nm的衍射角超过了90度']);
        end
    end
    scalar=tand(theta_r)./tand(theta_i);%列向量
    PVgr=scalar.*(PVc-PVi)+PVc;% 列向量.*行向量 +行向量 
end
%11.闪耀光栅闪耀波长grating_blaze(公式跟用法很相关,但结果挺接近的)
function lambda_b=grating_blaze(m,d,theta_b,theta_i)
    lambda_b=2.*d.*cosd(theta_b-theta_i).*sind(theta_b)./m;
end
%12.波长转RGB(网上抄的,有点问题,将就着用吧)
%可以研究章佳杰色度学科普文章，自己写一个，应用光学中也有相关知识
function RGB=getRGB(dWave)
    maxPix=1;gamma=1;
    waveArea = [380,440,490,510,580,645,780];
    minusWave = [0,440,440,510,510,645,780];
    deltWave = [1,60,50,20,70,65,35];
    for p=1:length(waveArea)
        if dWave<waveArea(p)
            break
        end
    end
    pVar=abs(minusWave(p)-dWave)/deltWave(p);
    rgbs=[[0,0,0];[pVar,0,1];[0,pVar,1];[0,1,pVar];[pVar,1,0];[1,pVar,0];[1,0,0];[0,0,0]];
    %在光谱边缘处颜色变暗
    if (dWave>=380) && (dWave<420)
        alpha = 0.3+0.7*(dWave-380)/(420-380);
    elseif (dWave>=420) && (dWave<701)
        alpha = 1.0;
    elseif (dWave>=701) && (dWave<780)
        alpha = 0.3+0.7*(780-dWave)/(780-700);
    else
        alpha = 0;       %非可见区
        %给红外光加点颜色，自780开始由黑变灰再变白。
        if dWave>780
            a=5e-5;
            grey=255-1/(a*(dWave-780)+1/255);
            RGB=(zeros(1,3)+grey)/255;
            return
        end
    end
    RGB=maxPix.*(rgbs(p,:).*alpha).^gamma;
end
%13.球坐标到xyz坐标
function [x,y,z]=ball2xyz(r,theta,phi)
    x=r.*sind(theta).*cosd(phi);
    y=r.*sind(theta).*sind(phi);
    z=r.*cosd(theta);
end

%14.绕轴旋转
%参考：罗德里格旋转公式（Rodrigues' rotation formula）
%DV旋转轴方向向量，PV旋转轴上一点的坐标，theta旋转角度，后三个是被旋转物体的xyz
function [obj1x,obj1y,obj1z,objDV,objPV]=axisrot(DV,PV,theta,obj1x,obj1y,obj1z,objDV,objPV)
    DV=DV./(dot(DV,DV).^0.5);%转为单位向量
    Rk=[0,-DV(3),DV(2);
        DV(3),0,-DV(1);
        -DV(2),DV(1),0]; %叉乘矩阵
    M=eye(3)+(1-cosd(theta)).*(Rk*Rk)+sind(theta).*Rk; %绕轴转矩阵
    obj1=[obj1x(:)'-PV(1);obj1y(:)'-PV(2);obj1z(:)'-PV(3)];
    temp=M*obj1;
    obj1x(:)=temp(1,:)+PV(1);
    obj1y(:)=temp(2,:)+PV(2);
    obj1z(:)=temp(3,:)+PV(3);
    objDV=(M*objDV')'; %向量没有位置，不用变换前先-PV，最后再+PV的操作
    objPV=(M*(objPV-PV)')'+PV;
end