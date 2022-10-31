function [Tnew]=conductstep(T,H,spacestep,time)
alpha = conductivity(T);
%boundary condition
Tnew = zeros(size(T,1),size(T,2));
for i=1:(size(T,2))
    for j=1:(size(T,1))
        Tnew(j,i) = T(j,i);
    end
end

for i=1:size(T,2)
    T(1,i) = geotherm(H(1,i));
    T(end,i) = geotherm(H(end,i));
end

for j=2:(size(Tnew,1)-1)
    T(j,end) = geotherm(H(j,end));
    %T(j,1) = geotherm(H(j,1));
end

for i=1:size(Tnew,2)
    Tnew(1,i) = geotherm(H(1,i));
    Tnew(end,i) = geotherm(H(end,i));
end

for j=2:(size(Tnew,1)-1)
    Tnew(j,end) = geotherm(H(j,end));
    %Tnew(j,1) = geotherm(H(j,1));
end


for j=2:(size(Tnew,1)-1)
    Tnew(j,1) = T(j,1)+alpha(j,1)*time*(T(j+1,1)+T(j-1,1)+2*T(j,2)-4*T(j,1))/(spacestep^2);
    %Tnew(j,1) = T(j,1)+alpha*time*(T(j+1,1)+T(j-1,1)+2*T(j,2)-4*T(j,1))/(spacestep^2);
end


for i=2:(size(Tnew,2)-1)
    for j=2:(size(Tnew,1)-1)
        Tnew(j,i) = T(j,i)+alpha(j,i)*time*(T(j+1,i)+T(j-1,i)+T(j,i-1)+T(j,i+1)-4*T(j,i))/(spacestep^2);
        %Tnew(j,i) = T(j,i)+alpha*time*(T(j+1,i)+T(j-1,i)+T(j,i-1)+T(j,i+1)-4*T(j,i))/(spacestep^2);
    end
end

end