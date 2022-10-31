function [Xnew,Ynew,Hnew]=rotateslab(X,Y,H,speed,time)
%This is the function to rotate the grid X,Y; speed is in radian/s;time is in s
degree = speed*time;
real_thickness=Y(1,1);
Xnew=zeros(size(Y,1),size(X,2));
Ynew=zeros(size(Y,1),size(X,2));
Hnew=zeros(size(Y,1),size(X,2));

for i=2:size(X,2)
    for j=1:size(Y,1)
        a(j,i) = atan(Y(j,i)/X(j,i));
        b(j,i) = a(j,i)-degree;
        r(j,i) = sqrt(X(j,i)^2+Y(j,i)^2); 
        Xnew(j,i) = r(j,i)*cos(b(j,i));
        Ynew(j,i) = r(j,i)*sin(b(j,i));
        Hnew(j,i) = -Ynew(j,i)+real_thickness;
    end
end

for j=1:size(Y,1)
        a(j,1) = pi/2;
        b(j,1) = a(j,1)-degree;
        r(j,1) = sqrt(X(j,1)^2+Y(j,1)^2); 
        Xnew(j,1) = r(j,1)*cos(b(j,1));
        Ynew(j,1) = r(j,1)*sin(b(j,1));
        Hnew(j,1) = -Ynew(j,1)+real_thickness;
end


Xnew(end,1) = 0;
Ynew(end,1) = 0;
Hnew(end,1) = real_thickness;

end