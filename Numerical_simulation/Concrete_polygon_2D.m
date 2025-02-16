clear
clc
close all

L=400;
W=100;
n1=[0.075  0.1213 0.0777 0.0939 0.123 0.256 0.124];
d=[15 10 7 5 4 3 2 1]./2;

c=1;
zarea=L*W;
area=zarea.*n1;
mj=zeros(size(n1,2),1);
geshu=zeros(size(n1,2),1);
for k=1:size(n1,2)
    while mj(k)<=area(k)
        jzb.n(c,1)=randi([6 10],1,1);
        jzb.th(c,1)=0;
        a=0;
        for j=1:jzb.n(c,1)
            b=2*pi/jzb.n(c,1)+(2*rand(1)-1)*0.05*2*pi/jzb.n(c,1);
            jzb.th(c,j+1)=b+a;
            a=jzb.th(c,j+1);
            jzb.r(c,j)=(d(k)+d(k+1))/2+(d(k)-(d(k)+d(k+1))/2)*(2*rand(1)-1)/2;
            jzb.R(c,j)=jzb.r(c,j)*1.05;
        end
        mji=0;
        for l=1:jzb.n(c,1)
            if l==jzb.n(c,1)
                jzb.mj(c,l)=1/2*jzb.R(c,1)*jzb.R(c,l)*sin(jzb.th(c,l+1)-jzb.th(c,l));
            else
          jzb.mj(c,l)=1/2*jzb.R(c,l)*jzb.R(c,l+1)*sin(jzb.th(c,l+1)-jzb.th(c,l));
            end
          mji= jzb.mj(c,l)+mji;
        end
        mj(k)=mj(k)+mji;
        geshu(k)=geshu(k)+1;
        jzb.th(c,jzb.n(c,1)+1)=2*pi;
        jzb.r(c,jzb.n(c,1)+1)=jzb.r(c,1);
        jzb.rmax(c,1)=max(jzb.r(c,:));
        jzb.R(c,jzb.n(c,1)+1)=jzb.R(c,1);
        jzb.Rmax(c,1)=max(jzb.R(c,:));
        c=c+1;
    end
end

range=[0 L;0 W];
plot([range(1,:) range(1,2) range(1,1) range(1,1)],[range(2,1) range(2,1) range(2,2) range(2,2) range(2,1)]);
hold on;
fill([range(1,:) range(1,2) range(1,1) range(1,1)],[range(2,1) range(2,1) range(2,2) range(2,2) range(2,1)],[0.75 0.75 0.75]); 
cum=0;
s1=c-1;
for i=1:1E4
    if cum==s1
        break;
    end
    sjzb.x(cum+1,1)=rand(1)*(range(1,2)-range(1,1))+range(1,1);
    sjzb.y(cum+1,1)=rand(1)*(range(2,2)-range(2,1))+range(2,1);
    r=jzb.Rmax(cum+1,1);
    if (sjzb.x(cum+1,1)-r>range(1,1) && sjzb.x(cum+1,1)+r<range(1,2)) &&...
            (sjzb.y(cum+1,1)-r>range(2,1) && sjzb.y(cum+1,1)+r<range(2,2)) 
        if cum==0
            cum=cum+1;
            for j=1:jzb.n(cum,1)+1
                xx1(cum,j)=jzb.r(cum,j)*cos(jzb.th(cum,j))+sjzb.x(cum,1);
                yy1(cum,j)=jzb.r(cum,j)*sin(jzb.th(cum,j))+sjzb.y(cum,1);
            end
        else
            sum=0;
            for j=1:cum
                dir=sqrt((sjzb.x(j,1)-sjzb.x(cum+1,1))^2+(sjzb.y(j,1)-sjzb.y(cum+1,1))^2);
                if dir<jzb.Rmax(j,1)+jzb.Rmax(cum+1,1)
                    break;
                else
                    sum=sum+1;
                end
            end
            if sum==cum
                cum=cum+1;
                for j=1:jzb.n(cum,1)+1
                xx1(cum,j)=jzb.r(cum,j)*cos(jzb.th(cum,j))+sjzb.x(cum,1);
                yy1(cum,j)=jzb.r(cum,j)*sin(jzb.th(cum,j))+sjzb.y(cum,1);
                end
            end
        end
    end
end

clear sum
for i=1:cum
    X1=xx1(i,:);
    Y1=yy1(i,:);
    plot(X1(X1~=0),Y1(Y1~=0),'-r');
    fill(X1(X1~=0),Y1(Y1~=0),[0.75 0.75 0.75]);
    x(i)=sum(xx1(i,:))/length(X1);  %圆心横坐标
    y(i)=sum(yy1(i,:))/length(Y1);  %圆心纵坐标
    radius(i)=sqrt((xx1(i,1)-x(i)).^2+(yy1(i,1)-y(i)).^2); %计算半径
    Rmax=max(radius);
    
    
     if radius(i)>0.8*Rmax
           plot(X1(X1~=0),Y1(Y1~=0),'-r');
           fill(X1(X1~=0),Y1(Y1~=0),[ 0.16 0.14 0.13]); %象牙黑
      else
             if radius(i)<=0.8*Rmax && radius(i)>0.6*Rmax
                plot(X1(X1~=0),Y1(Y1~=0),'-g'); 
                fill(X1(X1~=0),Y1(Y1~=0),[ 0.5 0.54 0.53]); %冷灰
             
              else
                if radius(i)<=0.6*Rmax && radius(i)>0.4*Rmax
                    plot(X1(X1~=0),Y1(Y1~=0),'-y'); 
                    fill(X1(X1~=0),Y1(Y1~=0),[0.44 0.5 0.41]); %石板灰
               else
                    if radius(i)<=0.4*Rmax && radius(i)>0.2*Rmax
                        plot(X1(X1~=0),Y1(Y1~=0),'-k'); 
                        fill(X1(X1~=0),Y1(Y1~=0),'w'); %白色
                  else
                    plot(X1(X1~=0),Y1(Y1~=0),'-c'); 
                    fill(X1(X1~=0),Y1(Y1~=0),[0.5 0.5 0.41]); %暖灰色
                    end
               end
            end
      end
    
end
axis image;
xlabel('\itx(mm)')
ylabel('\ity(mm)')
set(gca,'fontsize',8,'fontname','times')
print('混凝土多边形骨料.tif','-dtiff2','-r300')

fid=fopen('多边形骨料.txt','wt');
fprintf(fid,'pline\n');
fprintf(fid,'%.1f',0);
fprintf(fid,',');
fprintf(fid,'%.1f',0);
fprintf(fid,'\n');
fprintf(fid,'%.1f',L);
fprintf(fid,',');
fprintf(fid,'%.1f',0);
fprintf(fid,'\n');
fprintf(fid,'%.1f',L);
fprintf(fid,',');
fprintf(fid,'%.1f',W);
fprintf(fid,'\n');
fprintf(fid,'%.1f',0);
fprintf(fid,',');
fprintf(fid,'%.1f',W);
fprintf(fid,'\n');
fprintf(fid,'%.1f',0);
fprintf(fid,',');
fprintf(fid,'%.1f',0);
fprintf(fid,'\n');
fprintf(fid,'\n');


for i=1:cum
    X1=xx1(i,:);
    X1=X1(X1~=0);
    [m,n]=size(X1);
    fprintf(fid,'pline\n');
    for j=1:n-1
      fprintf(fid,'%.4f',xx1(i,j));
      fprintf(fid,',');
      fprintf(fid,'%.4f',yy1(i,j));
      fprintf(fid,'\n');
    end
    fprintf(fid,'%.4f',xx1(i,n));
    fprintf(fid,',');
    fprintf(fid,'%.4f',yy1(i,n));
    fprintf(fid,'\n');
    fprintf(fid,'\n');
end

