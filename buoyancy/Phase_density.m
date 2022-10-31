function [D]=Phase_density(T,P)    
    load('D_mesh.mat');
    if T<673.15
        D=3100;
    else
        V(1)=D_result_Earth(floor(P/5),floor((T-350-273.15)/50));
        V(2)=D_result_Earth(floor(P/5),ceil((T-350-273.15)/50));
        V(3)=D_result_Earth(ceil(P/5),ceil((T-350-273.15)/50));
        V(4)=D_result_Earth(ceil(P/5),floor((T-350-273.15)/50));
        y=P-floor(P/5)*5;
        x=T-floor((T-350-273.15)/50)*50-350-273.15;
        D = V(1)*(1-x/50)*(1-y/5)+V(2)*x/50*(1-y/5)+V(4)*(1-x/50)*y/5+V(3)*x/50*y/5;
    end
end