function [T_save,H_save]=TEST(length,spacestep,thickness,sinkrate,timestepnum,angle)
%length,spacestep,thickness (m)
%sinkrate mm/yr
%timestepnum
%angle rad

[X,Y,H,T,real_length,real_thickness]=gridgenerating(length,spacestep,thickness);
radrate = angle/(length*sin(angle)/(sinkrate*1E-3/365/24/3600));
time1 = angle/radrate/timestepnum;
for i=1:timestepnum/2
    [Xnew,Ynew,Hnew]=rotateslab(X,Y,H,radrate,time1*i);
    T=conductstep(T,Hnew,spacestep,time1);
end
T_C=T-273.15;
imagesc([0 length/1000],[0 thickness/1000],T_C,'Interpolation','bilinear');
colorbar;
caxis([273.15 1600]);
colormap turbo;
set(gcf,'position',[100,100,length/500,thickness/500]);
T_save = T;
H_save = Hnew;
