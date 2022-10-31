function [T]=geotherm(H)
% Geotherm use liner simpled from McKenzie & Bickle, 1988
%T in K;H in m
    if H<150000
        T = 1550*H/150000+273.15;
    elseif H>=150000
        T = 1550+273.15;
    end
end