function [D_delta,total]=Crust_density(T,H,P,length,spacestep,thickness,D_delta)
T_mantle = zeros(size(T,1),size(T,2));
D_mantle = zeros(size(T,1),size(T,2));
D_slab = zeros(size(T,1),size(T,2));
total=0;
    for j=1:size(T,1) 
        for i=1:size(T,2)
        	if (j-1)*spacestep/1000<=45
            	if T(j,i)-750-273.15<0
                    T_mantle(j,i) = geotherm(H(j,i));
                    D_mantle(j,i) = thermalD(T_mantle(j,i),P(j,i));
                    D_slab(j,i)  = Phase_density(T(j,i),P(j,i));
                    D_delta(j,i) = D_slab(j,i)-D_mantle(j,i);
                    total = total+D_delta(j,i);
                end
            end
        end           
    end
end

function density=thermalD(T,P)
    density_T = 3330*exp(-((2.832E-5)*(T-273.15)+0.5*(3.79E-8)*(T^2-273.15^2)));
    K_T=(1250*0.9+1330*0.1)+(-0.1875*0.9-0.1995*0.1)*(T-298);                              
    density = density_T/(1-4*P/(K_T+4*P))^0.25;
end
