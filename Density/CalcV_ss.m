function V = CalcV_ss(P,T,A,X)
% Calculates the volume of a phase (solid solution)
% ===============================================
% A is the structure of endmember
% P is the pressure [Kbar]
% T is the temperature [C]
%global HP98;
load('Venus_NCKFMASHTO.mat')
n=length(X)+1;

a=CalcA(A,P,T,X);
V=0;
for i=1:n
    V=V+a(2,i).*CalcV(P,T,HP98(A.endmember(1,i)));
%   V = V+a(2,i).*HP98(A.endmember(1,i)).V;
end

% volume due to nonideal mixing
%dP=1e-2*P+1e-2;
%a1=CalcA(A,P-dP,T,X);
%a2=CalcA(A,P+dP,T,X);

%V=V+sum(a(2,:).*((a2(1,:)-a1(1,:)))/2/dP);
 
end