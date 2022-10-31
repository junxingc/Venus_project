function [totalforce,totalslab]=videodensity(length,spacestep,thickness,sinkrate,timestepnum,angle,name)
%length,spacestep,thickness (m)
%sinkrate mm/yr
%timestepnum
%angle rad
load('D_mesh.mat');
avi_object = VideoWriter(append(name,'.avi'));
avi_object.FrameRate = timestepnum/5;
open(avi_object);

[X,Y,H,T,real_length,real_thickness]=gridgenerating(length,spacestep,thickness);
radrate = angle/(length*sin(angle)/(sinkrate*1E-3/365/24/3600));
time1 = angle/radrate/timestepnum;
n=[];
n1=[];
m=[];
D_delta=zeros(76,201);
totalforce=zeros(timestepnum,2);
totalslab=zeros(timestepnum,2);
h=waitbar(0,'please wait');
for i=1:timestepnum
    delete(n);
    delete(n1);
    delete(m);
    [Xnew,Ynew,Hnew]=rotateslab(X,Y,H,radrate,time1*i);
    T=conductstep(T,Hnew,spacestep,time1);
    for j=1:76
        for k=1:201
            P_Earth(j,k)=Earth_P_condition(Hnew(j,k));
            P_Venus(j,k)=Venus_P_condition(Hnew(j,k));
        end
    end
    D_delta=density_data(T,Hnew,P_Earth,length,spacestep,thickness);
    [D_delta,totalslab1]=Crust_density(T,Hnew,P_Earth,length,spacestep,thickness,D_delta);
    totalslab(i,1)=i;
    totalslab(i,2)=totalslab1;
    imagesc([0 length/1000],[0 thickness/1000],D_delta,'Interpolation','bilinear');hold on;
    [m,n,n1] = Crust_area(T,length,spacestep);hold on;
    colorbar;
    caxis([-500 500]);
    colormap(map);
    set(gcf,'position',[100,100,length/500,thickness/500]);
    M=getframe(gcf);
    writeVideo(avi_object,M);
    str=[num2str(i*100/timestepnum),'%'];
    waitbar(i/timestepnum,h,str);
    totalforce(i,1)=i;
    totalforce(i,2)=sum(sum(D_delta));
end
close(avi_object);
delete(h);
end
function [P]=Earth_P_condition(H)
    P = 9.8*3100*H/1E8;
end

function [P]=Venus_P_condition(H)
    P = 8.87*3100*H/1E8;
end