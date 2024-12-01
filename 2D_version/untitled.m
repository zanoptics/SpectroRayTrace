close all;clc;clear;
img=imread('Image_20221212222829948.bmp'); %读取图像
[height,width,d]=size(img); % 获取位图的高/宽/维数 2048*3076*3
imgray=rgb2gray(img);  %rgb图像转为灰度图像 
figure;imshow(imgray) % 显示灰度图
[X,Y]=meshgrid(1:width,1:height); % 产生供三维绘图的X,Y数据
figure;
surf(Y,X,imgray); shading interp;colorbar;
xlabel('x');ylabel('y')
daspect([1 1 1])
%%
clc
N=1600;
d=6.0535;
d/320*0.1;
%%
close all;clc;clear;
img=imread('Image_20221212222829948.bmp'); %读取图像
sp=1350; %样本
% figure(1);imshow(img)
imgray=double(rgb2gray(img));
I=mean(imgray(sp,:),1);
figure(2);plot(1:3072,I)
axis([0 3100 0 255])
for i=1:3072
    if I(i)==1
        i
        break;
    end
end
for i=3072:-1:1
    if I(i)==1
        i
        break;
    end
end
%%
close all;clc;clear;
img=imread('Image_20221212222829948.bmp'); %读取图像
sp=1140; %样本
%1590:1610
%光强
imgray=double(rgb2gray(img));
I=mean(imgray(sp,:),1);
%颜色
mid=mean(double(img(sp,:,:)),1);
cmap=reshape(mid,3072,3)./255;
%画
hold on
for i=700:2500
%     stem(473-(i-1921)*2.4/17,I(i),Color=cmap(i,:),LineWidth=1.5,Marker='.')%1550
    stem(i,I(i),Color=cmap(i,:),LineWidth=1.5,Marker='.')
end
set(gca,'color','none');
%%
close all;clc;clear;
img=imread('Image_20221212222829948.bmp'); %读取图像
imgray=rgb2gray(img);
filename='abc.gif';
num=100;
m=moviein(num);
for idx=1:num
    sp=8*idx+950;
    fig=figure('Visible','off');hold on;axis([0 3100 0 300]);
    I=mean(imgray(sp,:),1);
    plot(1:3072,I)
    m(:,idx)=getframe(fig);
    im=frame2im(m(:,idx));
    [A,cmap]=rgb2ind(im,256);
    if idx==1
        imwrite(A,cmap,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,cmap,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end
fig.Visible='on';
movie(m,4)
hold off;