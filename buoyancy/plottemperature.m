function []=plottemperature(X,Y,T)
    for i=1:size(T,2)
        for j=1:size(T,1)
            Temperature(j,i)=T(j,i);
            xvalues=X(j,i);
            yvalues=Y(j,i);
        end
    end   
    h=heatmap(T);
end