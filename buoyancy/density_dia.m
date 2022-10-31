length=400000;
spacestep=2000;
thickness=150000;
for j=1:76
        for k=1:201
            P_Earth(j,k)=Earth_P_condition(H_10(j,k));
            P_Venus(j,k)=Venus_P_condition(H_10(j,k));
        end
end
D_delta=density_data(T_10,H_10,P_Earth,length,spacestep,thickness);
[D_delta,totalslab1]=Crust_density(T_10,H_10,P_Earth,length,spacestep,thickness,D_delta);
imagesc([0 length/1000],[0 thickness/1000],D_delta,'Interpolation','bilinear');hold on;
[m,n,n1] = Crust_area(T_10,length,spacestep);hold on;
colorbar;
caxis([-500 500]);
colormap(map);
set(gcf,'position',[100,100,length/500,thickness/500]);
function [P]=Earth_P_condition(H)
    P = 9.8*3100*H/1E8;
end

function [P]=Venus_P_condition(H)
    P = 8.87*3100*H/1E8;
end