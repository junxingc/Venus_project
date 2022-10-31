for j=1:size(P_crust_Venus,1)
    for i=1:size(P_crust_Venus,2)
        if any(P_crust_Venus(j,i))
            scatter((T_20(j,i)-273.15),P_crust_Venus(j,i));hold on;
        end
    end
end