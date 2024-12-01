clc;clear
a=double(py.pyrgb.getRGB(550))
x=linspace(-pi,pi,51);
plot(x,sin(x),'Color',a,'LineWidth',2)