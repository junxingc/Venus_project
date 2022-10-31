function V = CalcV(P,T,A)
% Calculates the volume of an endmember
% ===============================================
% A is the structure of endmember
% P is the pressure [Kbar]
% T is the temperature [C]

%if strcmp(A.Sym,'MAKE')
%   global dataset; %#ok<TLEV>
%    V=0;
%    for i=1:length(A.H)
%        V=V+A.H(i)*CalcV(P,T,dataset(A.S(i)));
%    end
%    return
%end

V=A.V;alpha0=A.alpha0;alpha1=A.alpha1;
K=A.K; dKdt=A.dKdT;dKdp=A.dKdp;

T=T+273;                                             % temperature in Kelvin

% integral of VdP
if strcmp(A.Sym,'H2O')||strcmp(A.Sym,'CO2')||strcmp(A.Sym,'CO2')||strcmp(A.Sym,'CO')||strcmp(A.Sym,'CH4')||strcmp(A.Sym,'H2')
    [G,V]=CalcG_f(P,T-273,A);
    return;
elseif strcmp(A.Sym,'O2')
    V=0;
else
    V_T=V*exp(alpha0*(T-298)-2*alpha1*alpha0*(T^0.5-298^0.5));     % volume at T
    K_T=K+dKdt*(T-298);                              % bulk modulus at T
    V = V_T*(1-4*P/(K_T+4*P))^0.25;
    %V = V_T*(1-4*P/(K_T+4*P));
end
end