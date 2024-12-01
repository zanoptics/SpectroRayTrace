clear;
k0=0.6;
phi=90-abs(atand(-1/k0));
d=0.001/600; %1.67μm
%红光波长
lambda_red=760e-9;
theta_z1=asind(sind(phi)-lambda_red/d);
theta_f1=asind(sind(phi)+lambda_red/d);
hold on
plot([-5,5],[0,0],'k',LineWidth=2)
y=-50:05:50;
scatter(0,y,'k.')
fplot(@(x) tand(phi+90).*x,Color='#7E2F8E');
fplot(@(x) tand(-phi+90).*x,Color='#7E2F8E');%紫
fplot(@(x) tand(-theta_f1+90).*x,Color='#77AC30');%绿
fplot(@(x) tand(-theta_z1+90).*x,Color='#EDB120');%红
axis equal
%% 
clear;
k0=0.6;
phi=90-abs(atand(-1/k0));
d=0.001/600; %1.67μm
%红光波长
lambda_red=400e-9 : 0.5e-9: 550e-9;
theta_z1=asind(sind(phi)-lambda_red/d);
theta_f1=asind(sind(phi)+lambda_red/d);
hold on
plot([400e-9,760e-9],[0,0])
plot(lambda_red,theta_z1,'b'); %+1 蓝
plot(lambda_red,theta_f1,'k'); %-1 黑
%% 寻找k0下限，+1不能有小于0的角度 0.5124
clear;
for k0=1:-0.0001:0.25
    phi=90-abs(atand(-1/k0));
    d=0.001/600; %1.67μm
    %红光波长
    lambda_red=400e-9 : 0.5e-9: 760e-9;
    theta_z1=asind(sind(phi)-lambda_red/d);
    index=logical(find(theta_z1<=0));
    if(index)
        k0
        break;
    end
end
%% 寻找k0上限，-1不能有大于90的衍射角,大于90会出虚数 0.6483
clear;
for k0=0.25:0.0001:1
    phi=90-abs(atand(-1/k0));
    d=0.001/600; %1.67μm
    %红光波长
    lambda_red=400e-9 : 0.5e-9: 760e-9;
    theta_f1=asind(sind(phi)+lambda_red/d);
    index=logical(find(imag(theta_f1)~=0));
    if(index)
        k0
        break;
    end
end

