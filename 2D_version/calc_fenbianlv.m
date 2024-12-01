clc;format long
f=1.524e-1; %焦距 m
d=0.001/300; %每mm 300刻线 m
N=25*300; %7500
theta=15;
dl_dlambda=f/(d*cosd(theta)) %单位波长(m)在焦平面上相距多远(m)
dlambda=0.1e-9 %=0.1nm
xiansesan=dl_dlambda*dlambda*1e3 %线色散 波长相差0.1nm, 光谱上相距0.0047mm 
huatu=15.341/320*0.1
%色散本领
%波长相差0.1nm, 光谱上相距0.0047mm
%我们相机可以分辨大概0.0020mm 由已拍摄到的数据和仿真出来的结果对比得
%分辨本领,色散开了不一定能分辨开
%能分辨的最小波长差
delta_lambda=780/N  %能分辨的最小波长差就大概0.1nm

%%
clc
%买错了
huatu=3/320*0.1 %波长相差0.1nm, 光谱上相距0.009 mm 
%%
clc;
f=50.8e-3; %焦距 m
d=0.001/300; %每mm 300刻线 m
N=25*300; %7500
theta=25;
dl_dlambda=f/(d*cosd(theta)) %单位波长(m)在焦平面上相距多远(m)
dlambda=0.1e-9 %=0.1nm
xiansesan=dl_dlambda*dlambda*1e3 %线色散 波长相差0.1nm, 光谱上相距0.0047mm 
huatu=15.341/320*0.1
%色散本领
%波长相差0.1nm, 光谱上相距0.0047mm
%我们相机可以分辨大概0.0020mm 由已拍摄到的数据和仿真出来的结果对比得
%分辨本领,色散开了不一定能分辨开
%能分辨的最小波长差
delta_lambda=780/N  %能分辨的最小波长差就大概0.1nm