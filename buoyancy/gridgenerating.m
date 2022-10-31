function [X,Y,H,T,real_length,real_thickness]=gridgenerating(length,spacestep,thickness)
% This is the function to form grid. The length is the closest mutiple
% of the scacestep. parameters are in meter; The output the grid and the real length and real thickness 

%form grid
    length_n = floor(length/spacestep);
    thickness_n = floor(thickness/spacestep);
    real_length = length_n*spacestep;
    real_thickness = thickness_n*spacestep;
    X = zeros(thickness_n+1,length_n+1);
    Y = zeros(thickness_n+1,length_n+1);
    H = zeros(thickness_n+1,length_n+1);
    T = zeros(thickness_n+1,length_n+1);
%grid x,y
for i=1:length_n+1
    for j=1:thickness_n+1
        X(j,i) = (i-1)*spacestep; 
        Y(j,i) = real_thickness-(j-1)*spacestep;
        H(j,i) = -Y(j,i)+real_thickness;
        T(j,i) = geotherm(H(j,i));
    end
end

end