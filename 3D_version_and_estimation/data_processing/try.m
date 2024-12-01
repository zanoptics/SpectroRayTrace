clc; clear
load('酒精.txt')
load('无容器测量数据.txt')

alcohol=X__(98:1:3184);
air=X_______(98:1:3184);
figure;plot(air,'DisplayName',"air");hold on;
plot(alcohol,'DisplayName',"alcohol")
absorbance=log(air./alcohol);
plot(absorbance,'k','DisplayName',"absorb");
legend
%%
clear;clc;
load 酒精2.txt
sz=0.4;bt=0.57;lf=0.05:0.5:1; %大小、距下边距离，距左边距离
figure(1);set(gcf,'position',[192/3 108/3 192*5.5 108*5.5]);
axes('Position',[lf(1) bt sz sz]);
m=X__2;
[peaks,locs]=findpeaks(m,'minpeakdistance',length(m)/5);
plot(m);
hold on;
plot(locs,peaks,'r*'); 
T=round(sum(diff(locs))/4)
%%
clear;clc
load 酒精3.txt
plot(X__3);hold on 
load 水.txt
X_=X_(108:end);
plot(X_)
water=X_(1:20212);
alcohol=X__3(1:20212);
figure;plot(water,'DisplayName',"water");hold on;
plot(alcohol,'DisplayName',"alcohol");
absorbance=log(water./alcohol);
plot(absorbance,'k','DisplayName',"absorb");
legend
%% 水（3组）.txt
clc;format compact;clear;
water=importdata("水（3组）.txt");
alcohol=importdata("酒精（3组）.txt");
% figure;plot(water,'DisplayName',"water");hold on;
% plot(alcohol,'DisplayName',"alcohol");
% legend;
waterbar=avg(water,3);
alcoholbar=avg(alcohol,3);
[waterbar,alcoholbar]=align(waterbar,alcoholbar);
absorbance=log(waterbar./alcoholbar);
figure;plot(waterbar,'DisplayName',"water");hold on;
plot(alcoholbar,'DisplayName',"alcohol");
plot(absorbance,'k','DisplayName',"absorb");
legend;
%%
clear;clc;close all;format compact;
water=importdata("水（五组数据）.txt");
alcohol=importdata("酒精（五组数据）.txt");
figure;plot(water,'DisplayName',"water");hold on;
plot(alcohol,'DisplayName',"alcohol");
[water,alcohol]=align(water,alcohol);
absorbance=log(water./alcohol);
plot(absorbance,'k','DisplayName',"absorb");
legend
%%
clear;clc;close all;format compact;
water=importdata("水（3组）.txt");
alcohol=importdata("酒精（3组）.txt");
figure;plot(water,'DisplayName',"water");hold on;
plot(alcohol,'DisplayName',"alcohol");
[water,alcohol]=align(water,alcohol);
absorbance=log(water./alcohol);
plot(absorbance,'k','DisplayName',"absorb");
%%
% function waterbar=avg(water,n)
%     t=round(length(water)/n);
%     [mw1,imw1]=min(water(1:t));
%     [mw2,imw2]=min(water(t+1:2*t));
%     imw2=imw2+t;
%     T=imw2-imw1;
%     waterbar=zeros(T,1);
%     for ii=1:T
%         waterbar(ii)=(water(imw1+ii)+water(imw2+ii))/2;
%     end
% end

function [waterbar,alcoholbar]=align(waterbar,alcoholbar)
    lw=length(waterbar);
    la=length(alcoholbar);
    if lw>la
        waterbar=waterbar(1:la);
    else
        alcoholbar=alcoholbar(1:lw);
    end
end