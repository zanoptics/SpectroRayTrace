clc;format compact;clear;close all

%11.闪耀光栅闪耀波长grating_blaze(公式跟用法很相关,但结果挺接近的)
function lambda_b=grating_blaze(m,d,theta_b,theta_i)
    lambda_b=-2.*d.*cosd(theta_i+theta_b).*sind(theta_b)./m;
end