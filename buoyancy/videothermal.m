function []=videothermal(length,spacestep,thickness,sinkrate,timestepnum,angle,name)
%length,spacestep,thickness (m)
%sinkrate mm/yr
%timestepnum
%angle rad
avi_object = VideoWriter(append(name,'.avi'));
avi_object.FrameRate = timestepnum/5;
open(avi_object);

[X,Y,H,T,real_length,real_thickness]=gridgenerating(length,spacestep,thickness);
radrate = angle/(length*sin(angle)/(sinkrate*1E-3/365/24/3600));
time1 = angle/radrate/timestepnum;
n1=[];
c=[];
h=waitbar(0,'please wait');
for i=1:timestepnum
    delete(n1);
    delete(c);
    [Xnew,Ynew,Hnew]=rotateslab(X,Y,H,radrate,time1*i);
    T=conductstep(T,Hnew,spacestep,time1);
    T_C=T-273.15;
    c=imagesc([0 length/1000],[0 thickness/1000],T_C,'Interpolation','bilinear');hold on;
    [n1] = Crust_area(T,length,spacestep);hold on;
    colorbar;
    caxis([273.15 1600]);
    colormap turbo;
    set(gca,'YDir','reverse');
    set(gcf,'position',[100,100,length/500,thickness/500]);
    xlim([0 length/1000]);
    M=getframe(gcf);
    writeVideo(avi_object,M);
    str=[num2str(i*100/timestepnum),'%'];
    waitbar(i/timestepnum,h,str);
end
close(avi_object);
delete(h);


