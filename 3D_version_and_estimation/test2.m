clc
%直角坐标转球坐标
PVn=PVn-[2,2,2]
theta=acosd(PVn(3))
phi=acosd(PVn(1)./sind(theta))
[a,b,c]=cart2sph(PVn(1),PVn(2),PVn(3))
a.*180./pi
90-b.*180/pi