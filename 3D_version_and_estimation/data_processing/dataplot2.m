clc;clear;close all;format compact
figure(1)
% for ii=[29,32]
%     eval(['plot(readox("光电无东西数据2(',num2str(ii),').txt"))'])
%     legend
%     hold on
% end
% 
% for ii=6:8
%     alcohol=readox(['光电酒精(',num2str(ii),').txt']);
%     alcohol=alcohol(1:3059);
%     plot(alcohol,"-","LineWidth",1,"DisplayName",['99alcohol-',num2str(ii)]);hold on;legend;
%     try
%         alcohol99=alcohol99+alcohol./3;
%     catch
%         alcohol99=alcohol./3;
%     end
% end

for ii=12:16
    alcohol=readox(['光电酒精(',num2str(ii),').txt']);
    alcohol=alcohol(1:3059);
    plot(alcohol,"-","LineWidth",1,"DisplayName",['50alcohol-',num2str(ii)]);hold on;legend;
    try
        alcohol50=alcohol50+alcohol./3;
    catch
        alcohol50=alcohol./3;
    end
end

% for ii=9:11
%     water=readox(['光电水(',num2str(ii),').txt']);
%     water=water(1:3059);
%     plot(water,"-","LineWidth",1,"DisplayName",['water-',num2str(ii)]);hold on;legend;
%     try
%         water100=water100+water./3;
%     catch
%         water100=water./3;
%     end
% end
%%
clc;clear;close all;format compact
figure(1)
% alc99_mean=BatchPlot('酒精','99',6:8);
alc50_mean=BatchPlot('酒精','50',[18:20]);
% wat100_mean=BatchPlot('水','100',9:11);
%%
figure(2)
absorb_alc99=log(wat100_mean./alc99_mean);
absorb_alc50=log(wat100_mean./alc50_mean);
absorb_alc0=log(wat100_mean./wat100_mean);
absorb_alc99=absorb_alc99(end:-1:1);
absorb_alc50=absorb_alc50(end:-1:1);
plot(absorb_alc99,'DisplayName','99酒精');hold on;legend
plot(absorb_alc50,'DisplayName','50酒精');
plot(absorb_alc0,'DisplayName','0酒精');
title('吸光度')
%%
clc;clear;close all;format compact
figure(2)
alcohol=readox("光电酒精(6).txt");
water=readox("光电水(11).txt");
absorbance=log(water./alcohol);
plot(alcohol,"DisplayName","alcohol","LineWidth",1,"MarkerSize",1);hold on
plot(water,'DisplayName','water');legend
figure(3)
plot(absorbance)


% hold on
% %of 创建一个低通滤波器，保留低频成分，去除高频成分
% datafft=fft(data);
% filter = ones(size(datafft)); % 创建与信号频谱相同大小的滤波器
% cutoff_frequency = 5; % 设置截止频率
% filter(cutoff_frequency:end-cutoff_frequency+2) = 0; % 去除高频成分
% filtered_signal_fft = datafft .* filter;
% filtered_signal = ifft(filtered_signal_fft);
% % 绘制原始信号
% plot(data,'--');
% title('Original Signal');


function data=readox(filename)
    fileID = fopen(filename);
    a=cell2mat(textscan(fileID,"%x"));
    data=a(1:2:end)+a(2:2:end).*(16.^2);
    data=double(data);
%     size(data)
end

function data_avg=BatchPlot(SampleName,conc,index)
    len=length(index);
    for ii=index
        data=readox(['光电',SampleName,'(',num2str(ii),').txt']);
        data=data(1:3059);
        plot(data,"-","LineWidth",1,"DisplayName",[conc,SampleName,num2str(ii)]);hold on;legend;
        try
            data_avg=data_avg+data./len;
        catch
            data_avg=data./len;
        end
    end
end









