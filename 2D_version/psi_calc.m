% %输入为lambda,计算ψ为多少时,可把对应波长的光聚焦在焦点
clear;
k0=0.6;
lambda=[400e-9 475e-9 550e-9];
d=0.001/600; %1.67μm 即600条/mm
a=atand(1/k0);
b=acotd(1/k0);
A=cosd(b)+sind(a);
B=sind(b)-cosd(a);
C=lambda./(d*(A^2+B^2)^0.5);
w=atand(B/A);
psi=asind(C)-w;
% k0=0.6;
% lambda=400e-9;
% d=0.001/600;
% C=lambda/d;
% a=atand(1/k0);
% b=acotd(1/k0);
% A=sind(a)+cosd(b);
% B=sind(b)-cosd(a);
% psi=asind(C/sqrt(A^2+B^2))-atand(B/A);

