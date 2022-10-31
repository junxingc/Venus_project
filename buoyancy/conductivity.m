function alpha=conductivity(T)
% k thermal conductivity, c heat capacity, density
alpha = zeros(size(T,1),size(T,2));
k = zeros(size(T,1),size(T,2));
c = zeros(size(T,1),size(T,2));

for i=1:size(T,2)
    for j=1:size(T,1)
        k(j,i) = 5.3/(1+0.0015*(T(j,i)-273.15))+(-1.0365E-4)*T(j,i)+(2.2451E-7)*(T(j,i)^2)+(-3.4071E-11)*(T(j,i)^3);
        c(j,i) = (233.18-1801.6*(T(j,i)^(-1/2))-26.794E7*(T(j,i)^(-3)))*(890/140.7) + (252-2013.7*(T(j,i)^(-1/2))-6.219E7*(T(j,i)^(-3)))*(110/203.8);
        alpha(j,i) = k(j,i)/c(j,i)/3300;
    end
end
end