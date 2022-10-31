function []=plotslabedge(X,Y,color)
    x1= [X(1,1) X(1,end)  X(end,end) X(end,1) X(1,1)];hold on;
    y1= [Y(1,1) Y(1,end)  Y(end,end) Y(end,1) Y(1,1)];hold on;
    line(x1,y1,'Color',color);hold on;    
end