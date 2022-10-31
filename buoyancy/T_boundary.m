function [j_saved,i_saved]=T_boundary(T,spacestep,T_limit)
    j_saved = zeros(1);
    i_saved = zeros(1);    
    for j=1:size(T,1) 
            if T(j,end)-T_limit-273.15>0
                for i=2:size(T,2)
                    if (j-1)*spacestep/1000<=45
                           if (T(j,i-1)-T_limit-273.15)*(T(j,i)-T_limit-273.15)<=0
                               j_saved(end+1,1) = j;
                               i_saved(end+1,1) = i-1;
                           end
                    end
                end
            else
                if (j-1)*spacestep/1000<=45
                    j_saved(end+1,1) = j;
                    i_saved(end+1,1) = size(T,2);
                end
            end
    end
end