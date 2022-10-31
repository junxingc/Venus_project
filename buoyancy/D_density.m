function []=D_density(T,H,P,length,spacestep,thickness)
T_mantle = zeros(size(T,1),size(T,2));
D_mantle = zeros(size(T,1),size(T,2));
D_slab = zeros(size(T,1),size(T,2));
D_delta = zeros(size(T,1),size(T,2));
for i=1:size(T,2)
    for j=1:size(T,1)    
        T_mantle(j,i) = geotherm(H(j,i));
        D_mantle(j,i) = thermalD(T_mantle(j,i),P(j,i));
        D_slab(j,i)  = thermalD(T(j,i),P(j,i));
        D_delta(j,i) = D_slab(j,i)-D_mantle(j,i);
    end
end
imagesc([0 length/1000],[0 thickness/1000],D_delta,'Interpolation','bilinear');
colorbar;
caxis([3000 3500]);
colormap turbo;
set(gcf,'position',[100,100,length/500,thickness/500]); hold on;
end

function density=thermalD(T,P)
    density_T = 3330*exp(-((2.832E-5)*(T-273.15)+0.5*(3.79E-8)*(T^2-273.15^2)));
    K_T=(1250*0.9+1330*0.1)+(-0.1875*0.9-0.1995*0.1)*(T-298);                              
    density = density_T/((1-4*P/(K_T+4*P))^0.25);
end

