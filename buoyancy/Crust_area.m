function [m,n,n1]=Crust_area(T,length,spacestep)
    j_saved = zeros(1);
    i_saved = zeros(1);    
    for j=1:size(T,1) 
            if T(j,end)-750-273.15>0
                for i=2:size(T,2)
                    if (j-1)*spacestep/1000<=45
                           if (T(j,i-1)-750-273.15)*(T(j,i)-750-273.15)<=0
                               j_saved(end+1,1) = (j-1)*spacestep/1000;
                               i_saved(end+1,1) = (i-2-(T(j,i-1)-750-273.15)/(T(j,i)-T(j,i-1)))*spacestep/1000;
                           end
                    end
                        if (j-1)*spacestep/1000<=45+spacestep/1000 && (j-1)*spacestep/1000>45
                           if (T(j,i-1)-750-273.15)*(T(j,i)-750-273.15)<=0
                               j_temp = (j-1)*spacestep/1000;
                               i_temp = (i-2-(T(j,i-1)-750-273.15)/(T(j,i)-T(j,i-1)))*spacestep/1000;
                               j_saved(end+1,1) = 45;
                               i_saved(end+1,1) = i_temp-(i_temp-i_saved(end,1))*(j_temp-45)/(j_temp-j_saved(end-1,1));
                           end
                        end
                end
            else
                if (j-1)*spacestep/1000<=45
                    j_saved(end+1,1) = (j-1)*spacestep/1000;
                    i_saved(end+1,1) = length/1000;
                end
                j_saved(end+1,1) = 45;
                i_saved(end+1,1) = length/1000;
            end
    end
        j_saved(1,1) = j_saved(2,1);
        i_saved(1,1) = 0;
        j_saved(end+1,1) = 45;
        i_saved(end+1,1) = 0;
        %{
        j_saved = j_saved';
        i_saved = i_saved';
        values = spcrv([[i_saved(1) i_saved i_saved(end)];[j_saved(1) j_saved j_saved(end)]],3);
        plot(values(1,:),values(2,:)); hold on;
        %}        
        m=patch([0 (i_saved(1:end-1))' length/1000 length/1000],[0 (j_saved(1:end-1))' 45 0],'white','LineStyle','none'); hold on;
        alpha(m,0);
        hatch(m,[30 5 1],'k-');
        n = plot(i_saved,j_saved,'black','LineWidth',2'); hold on;
        n1 = plot([0 (i_saved(1:end-1))' length/1000 length/1000 0],[0 (j_saved(1:end-1))' 45 0 0],'black','LineWidth',2'); hold on;
        %n = fill(i_saved,j_saved,'white','LineStyle','none'); hold on;
end