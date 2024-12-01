%==========================================================================
%����׷�� �������߾���һ��ƽ�淴�侵,���䵽̽������
%22.8.31
%��ѩ
%==========================================================================
%�����ԴS��ƽ�淴�侵M��̽����D
%==========================================================================
% x=-100:0.0001:100;
tic

% x0=0;
% y0=100;%���ù�Դ����S(x0,y0)
x0=-30;
y0=0;%���ù�Դ����S(x0,y0)

% Ax=0; Ay=30; %A������
% % Ax=0; Ay=10; %A������
% Bx=50; By=20;%B������
Ax=20; Ay=60; %A������
Bx=70; By=40;%B������
%M=(By-Ay)/(Bx-Ax)*(x-Ax)+Ay;%��A��B��ȷ��ƽ�淴�侵����
km=(By-Ay)/(Bx-Ax);%ƽ�淴�侵����б��
jm=Ay-km*Ax;%ƽ�淴�侵���̽ؾ�
z=Ax:0.001:Bx;
M=km*z+jm;%ƽ�淴�侵����

% Cx=90; Cy=40; %C������
% Cx=90; Cy=90; %C������
% Dx=40; Dy=80;%D������
Cx=100; Cy=-20; %C������
Dx=40; Dy=-40; %D������
%D=(Dy-Cy)/(Dx-Cx)*(x-Cx)+Cy;%��C��D��ȷ��ccd̽���淽��
kd=(Dy-Cy)/(Dx-Cx);%̽��������б��
jd=Cy-kd*Cx;%̽�������̽ؾ�

%==========================================================================
%���ɲ���������������ߣ�����ƽ�澵��̽����
%==========================================================================
nn=10;%���ȡ�Ĺ��߸���
Px=z(randperm(numel(z),nn));  %��AB�����ȡ����㣬��Ϊ��Դ���������ƽ�淴�侵�Ľ���
Py=km.*Px+jm; %��AB�����ȡ����㣬��Ϊ��Դ���������ƽ�淴�侵�Ľ���
hold on;
for i=1:nn
   plot([x0,Px(i)],[y0,Py(i)],'Linewidth',1);%����5���������
end
daspect([1 1 1]);
plot([Ax,Bx],[Ay,By],'k','Linewidth',2);%����ƽ�淴�侵
plot([Cx,Dx],[Cy,Dy],'r','Linewidth',2);%����̽����ƽ��

%==========================================================================
%�������ǣ�����������б�ʼ��ؾ�
%��ⷴ�������̽���潻�㣬�����������
%==========================================================================
alpha=zeros(1,nn);%��������
inan=zeros(1,nn);%�����
for i=1:nn
    l=norm([x0-Ax,y0-Ay]);%��SA��ģ
    m=norm([x0-Px(i),y0-Py(i)]);%��SP_i��ģ
    n=norm([Ax-Px(i),Ay-Py(i)]);%��AP_i��ģ
    alpha(i)=acos((m^2+n^2-l^2)/(2*m*n));%�����Ҷ����������ǵ����
    inan(i)=(pi/2-alpha(i))/pi*180;%����ǣ������=�����
end
% inan %��������

k_AD=(Dy-Ay)/(Dx-Ay);
k_AC=(Cy-Ay)/(Cx-Ax);%̽������̽�⵽�ķ������б�ʷ�ΧΪk_AC<k<k_AD

kr=zeros(1,nn);%�������б��
jr=zeros(1,nn);%������߽ؾ�
X=zeros(1,nn);%���������̽���潻��
Y=zeros(1,nn);

% for i=1:nn
%     k1=(Py(i)-y0)/(Px(i)-x0);%��������ߵ�б��
%     if k1>0
%        kr(i)=-tan(alpha(i)-atan(km));%������߷���б��
%     else
%        kr(i)=tan(atan(km)+alpha(i));%������߷���б��
%     end
%        jr(i)=Py(i)-kr(i)*Px(i);%������߽ؾ�
%      if kr(i)>=k_AC & kr(i)<=k_AD
%             [X(i),Y(i)]=linecross(kr(i),jr(i),kd,jd);%���������̽���潻�� 
%             plot([Px(i),X(i)],[Py(i),Y(i)],'Linewidth',1);%�����������
%         else
%             xx=Px(i)*3;
%             yy=kr(i)*xx+jr(i);
%             plot([Px(i),xx],[Py(i),yy],'Linewidth',1);%������߲��ܱ�̽����̽�⵽
%      end
% end
 
for i=1:nn
    k1=(Py(i)-y0)/(Px(i)-x0);%��������ߵ�б��
    if k1>0
       kr(i)=-tan(alpha(i)-atan(km));%������߷���б��
       jr(i)=Py(i)-kr(i)*Px(i);%������߽ؾ�
       if kr(i)>=-1/km & kr(i)<=k_AC
%          if kr(i)>=k_AC & kr(i)<=k_AD
            [X(i),Y(i)]=linecross(kr(i),jr(i),kd,jd);%���������̽���潻�� 
            plot([Px(i),X(i)],[Py(i),Y(i)],'Linewidth',1);%�����������
        else
            xx=Px(i)*3;
            yy=kr(i)*xx+jr(i);
            plot([Px(i),xx],[Py(i),yy],'Linewidth',1);%������߲��ܱ�̽����̽�⵽
        end
    else
        if km<0
            kr(i)=tan(atan(km)+alpha(i));%������߷���б��
        else
            if alpha(i)>pi/2
                kr(i)=-tan(pi-alpha(i)-atan(km));%������߷���б��
            else
                kr(i)=tan(atan(km)+alpha(i));%������߷���б��
            end
        end
       jr(i)=Py(i)-kr(i)*Px(i);%������߽ؾ�
       if kr(i)>=k_AC & kr(i)<=k_AD
             [X(i),Y(i)]=linecross(kr(i),jr(i),kd,jd);%���������̽���潻�� 
             plot([Px(i),X(i)],[Py(i),Y(i)],'Linewidth',1);%�����������
       else
           if alpha(i)>pi/2
               xx=-Px(i);
               yy=kr(i)*xx+jr(i);
               plot([Px(i),xx],[Py(i),yy],'Linewidth',1);%������߲��ܱ�̽����̽�⵽    
           else
               xx=Px(i)*1.002;
             yy=kr(i)*xx+jr(i);
             plot([Px(i),xx],[Py(i),yy],'Linewidth',1);%������߲��ܱ�̽����̽�⵽
           end
       end          
    end
end
hold off;
axis equal
toc
%% ��֪����ֱ�ߵ�б�ʺͽؾ࣬�󽻵�����
function [x,y]=linecross(k1,b1,k2,b2)
  x=[];
  y=[];
  if k1==k2&b1==b2
      disp('�غ�');
  elseif k1==k2&b1~=b2
      disp('�޽���');
  else
     x=(b2-b1)/(k1-k2);
     y=k1*x+b1;
  end
end


