clear; clc
ki=0.144356955380577;
a=atand(ki);
lambda=497e-9;
d=0.001/300; %3.33μm 即300条/mm
psi=acosd(lambda./(2.*d.*cosd(a)))
%为什么求出来不是个负的
% 85.6802  497