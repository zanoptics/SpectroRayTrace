clc;format compact;
%判断四点共面
%球面镜
% common_plane=[PVi,1;PVn,1;PVb,1;PVr,1;];%不太懂原理
% det(common_plane)
det([PVn-PVi;PVb-PVi;PVr-PVi])
%入射角等于反射角
angle_eq=@(A,B) dot(A,B)./(dot(A,A).*dot(B,B)).^0.5;
angle_eq(PVi-PVn,PVb-PVn).*180./pi==angle_eq(PVr-PVn,PVb-PVn).*180./pi
