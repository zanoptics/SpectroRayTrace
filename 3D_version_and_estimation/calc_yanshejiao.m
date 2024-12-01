clc
m=1;
d=1/1200;
lambda=1000e-6
% theta_r=asind(m.*lambda./d-sin(50))%列向量
theta_r=asind(m.*lambda./d-sind(50))
%使用同侧的衍射光
% i       r
% 50    25.72
% 36.87 36.87
% 15    70.25
%%
%看资料计算区
%单位m
c=3e8;
lambda=[780e-9,2500e-9,1000000e-9];
f=c./lambda
%%
clc;format compact
% wave_num=[3425.6,2974.8,1055.8]; %波数 cm^-1
wave_num=[1000:1000:11000]; %波数 cm^-1
wave_len=1./wave_num *1e7  %波长nm
show=[wave_num',wave_len'];