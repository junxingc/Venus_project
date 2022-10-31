for j=1:76
    for i=1:201
        P_Earth_20(j,i)=Earth_P_condition(H_20(j,i));
        P_Venus_20(j,i)=Venus_P_condition(H_20(j,i));
    end
end

function [P]=Earth_P_condition(H)
    P = 9.8*3100*H/1E8;
end

function [P]=Venus_P_condition(H)
    P = 8.87*3100*H/1E8;
end