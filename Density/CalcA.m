function RTlna=CalcA(SS,P,T,X)

% calculate the chemical potentials of endmembers in a solid solution, in
% terms of RTlna [KJ]
% =========================
% SS is the structer of solid solution
% p is the proportions of endmembers in the solid solution
% P is the pressure [Kbar]
% T is the temperature [C]
% X is the composational variables: [x1,x2,...,xn]
% format of RT ln a: [RTlna_1,RTlna_2,RTlna_3;#1,#2,#3], where #1 #2 #3 are
% the index number in the dataset

n=SS.endnum;
RTlna=zeros(3,n);
RTlna(3,:)=SS.endmember;
T=T+273.15;%R=8.314472e-3;
R=8.31e-3;

if n==1
    if strcmp(SS.sys,'PM') % polymorph
        global dataset;
        T=T-273.15;
        RTlna(2,1)=1;
        m=length(SS.alpha)+1;
        G=zeros(m,1);
        G(1)=CalcG(P,T,dataset(SS.endmember));
        for i=2:length(G)
            G(i)=CalcG(P,T,dataset(SS.alpha(i-1)));
        end
        [RTlna(1,1),a]=min(G-G(1));
    elseif strcmp(SS.minrl,'H2O')
        RTlna(1,1)=R*T*log(0.7);
        RTlna(2,1)=1;
    else
        RTlna(1,1)=0;
        RTlna(2,1)=1;
    end
    return;
end
% =======================================================
% check if the compositional variables are reasonable to avoid solvi
%for i=1:n-1
%    if X(i)<SS.vari_lower(i)||X(i)>SS.vari_upper(i)
%        fprintf('unreasonable compositional variables, try another starting guess\n');
%        fprintf(SS.minrl);fprintf('\n');
%        fprintf(SS.variables{i});fprintf('= %g\n',X(i));
%        fprintf('That variable should lie in [%g,%g]\n',SS.vari_lower(i),SS.vari_upper(i))
%        RTlna=[];
%        return;
%    end
%end

% =======================================================
% proportion of endmembers
p=ones(1,n);
for i=1:n
    f=SS.p{i};
    p(i)=f(X);
end
clear f;
RTlna(2,:)=p;

% if abs(sum(p)-1)>1e-14
%     fprintf(SS.minrl);
%     fprintf('Error: the total proportions of endmembers do not make 100%\n\n');
%     return;
% end

if strcmp(SS.minrl,'liq')
%% Aquous Fluid
    if strcmp(SS.minrl,'liq')
        p_f=zeros(1,5);
        p_f(1:3)=p;
        p_f(4:5)=0;
        
        A_f=CalcA_f(P,T-273,p_f);
        RTlna(1,:)=A_f(1:3);
        RTlna(1,:)=R*T*log(RTlna(1,:));
    end
    
elseif strcmp(SS.model,'rh_ge')
%% rhombohedral oxide model by Ghiorso and Evans 2008
    
    X_temp=SS.a_ideal(X);
    if ~prod((X_temp>=0)*1)||sum(X_temp)>1
        RTlna(1,:)=1i;
        return
    end

    RTlna(1,:)=CalcA_rhoxide(T-273.15,X_temp);
    RTlna(1,:)=R*T*log(real(RTlna(1,:)));
    
%newly added carbonate
elseif strcmp(SS.model,'cab')
    
%    a_ideal = zeros(1,n);
%    for i=1:n
%        f=SS.a_ideal{i};
%        a_ideal(i)=f(X);
%    end
%    RTlna_ideal = R*T*log(a_ideal);
    
%    if ~prod((X_temp>=0)*1)||sum(X_temp)>1
%        RTlna(1,:)=1i;
%        return
%    end
    RTlna(1,:)=R*T*log(CalcA_cab(T,X));
%newly added carbonate end 

else
%% solid minerals
% =======================================================
% ideal activity
    a_ideal=zeros(1,n);
    for i=1:n
        f=SS.a_ideal{i};
        a_ideal(i)=f(X);
    end
    RTlna_ideal=R*T*log(a_ideal);

% =======================================================
% Margular coefficient

    RTlna_W=zeros(1,n);
    W=SS.W+T*SS.dWdT+P*SS.dWdP;    

    if strcmp(SS.model,'ideal')
        RTlna_W=zeros(1,n);

    elseif strcmp(SS.model,'sf')
        for l=1:n
            gamma=0;
            for i=1:n-1
                if i==l
                    qi=1-p(i);
                else
                    qi=-p(i);
                end
                for j=i+1:n
                    if j==l
                        qj=1-p(j);
                    else
                        qj=-p(j);
                    end
                    gamma=gamma-qi*qj*W(i,j-1);
                end
            end
            RTlna_W(1,l)=gamma;
        end
    
    elseif strcmp(SS.model,'asf')
        alpha=SS.alpha+T*SS.dalphadT+P*SS.dalphadP;
        s=sum(p.*alpha);
        phi=p.*alpha/s;
    
        for l=1:n
            gamma=0;
            for i=1:n-1
                if i==l
                    qi=1-phi(i);
                else
                    qi=-phi(i);
                end
                for j=i+1:n
                    if j==l
                        qj=1-phi(j);
                    else
                    qj=-phi(j);
                    end
                    gamma=gamma-qi*qj*W(i,j-1)*2*alpha(l)/(alpha(i)+alpha(j));
                end
            end
            RTlna_W(1,l)=gamma;
        end
    elseif strcmp(SS.model,'OMG_nonideal')
        for l=1:n
            f=SS.a_nonideal{l};
            RTlna_W(1,l)=f(p,X);
        end
    end

% DQF
    RTlna_DQF=SS.DQF+T*SS.dDQFdT+P*SS.dDQFdP;

% =======================================================
    RTlna(1,:)=RTlna_ideal+RTlna_W+RTlna_DQF;

end


end

    

%%
function a=CalcA_rhoxide(T,x)
% Calculate activities of rhombohedral oxide solid solutions (Ghiorso and
% Evans, 2008)

% endmembers: hematite, ilmenite, geikielite, pyrophanite, corundum
% compositional parameters

X=x(1);Y=x(2);Z=x(3);C=x(4);
H=1-X-Y-Z-C;
T=T+273.15;R=8.31e-3;
if Z==0&&Y==0
    m=1;
elseif Z==0||Y==0
    m=2;
else
    m=3;
end

% model parameters
W_hm_il  = 22.536;
DW_hm_il = -5.627;
% DW_hm_il = -5.064;          % Fe-Ti subsystem
W_hm_il2 = -0.833;
% W_hm_il2 = -1.960;          % Fe-Ti subsystem
DW_hm_il3= 0;
W_il_cr  = 80;%W_hm_il;
DW_il_cr = DW_hm_il;
W_il_il  = 69.908;
D2W_il_il= 12.756;
alpha    = 0.073021;

W_hm_gk  = 22.536;
DW_hm_gk = -5.627;
W_hm_gk2 = W_hm_il2;
DW_hm_gk3= 0;
W_gk_cr  = W_il_cr;
DW_gk_cr = DW_il_cr;
W_il_gk  = 2.6;
W_il_gk_T= 88.1;
W_gk_gk  = W_il_il;
D2W_gk_gk= D2W_il_il;

W_hm_py  = 22.536;
DW_hm_py = -5.627;
W_hm_py2 = W_hm_il2;
DW_hm_py3= 0;
W_py_cr  = W_il_cr;
DW_py_cr = DW_il_cr;
W_il_py  = 2.2;
W_il_py_T= 30.244;
W_py_py  = W_il_il;
D2W_py_py= D2W_il_il;

W_gk_py  = W_il_gk;
W_gk_py_T= W_il_gk_T;

W_hm_cr  = 116-32e-3*T;     % Majzlan, Navrotsky, and Evans, 2002;
DW_hm_cr = 0;

% ordering parameters

s=X-1e-8;t=Y-1e-8;u=Z-1e-8;    % ordered

dGds = @(s,t,u)...
    +R*T/2*(log(X+s)-log(X-s)+log(X+Y+Z+s+t+u)-log(X+Y+Z-s-t-u))...
    +1/2*t*(W_il_gk-W_il_gk_T)...
    +1/2*u*(W_il_py-W_il_py_T)...
    -2*s*(W_il_il/4+H*(DW_hm_il/2+W_hm_il2/4)-1/2*C*DW_il_cr)...
    +D2W_il_il/2*s*(X^2-2*s^2);
dGdt = @(s,t,u)....
    +R*T/2*(log(Y+t)-log(Y-t)+log(X+Y+Z+s+t+u)-log(X+Y+Z-s-t-u))...
    +1/2*s*(W_il_gk-W_il_gk_T)...
    +1/2*u*(W_gk_py-W_gk_py_T)...
    -2*t*(W_gk_gk/4+H*(DW_hm_gk/2+W_hm_gk2/4)-1/2*C*DW_gk_cr)...
    +D2W_gk_gk/2*t*(Y^2-2*t^2);
dGdu = @(s,t,u)....
    +R*T/2*(log(Z+u)-log(Z-u)+log(X+Y+Z+s+t+u)-log(X+Y+Z-s-t-u))...
    +1/2*s*(W_il_py-W_il_py_T)...
    +1/2*t*(W_gk_py-W_gk_py_T)...
    -2*u*(W_py_py/4+H*(DW_hm_py/2+W_hm_py2/4)-1/2*C*DW_py_cr)...
    +D2W_py_py/2*u*(Z^2-2*u^2);
deriG={dGds,dGdt,dGdu};

F=zeros(m,1);
A=zeros(m,m);
check=1;
x=[s,t,u]';
while check>1e-16
    x_old=x;
    for i=1:m
        F(i)=deriG{i}(x(1),x(2),x(3));
    end
    for i=1:m
        dx=x;
        dx(i)=-x(i)*1e-8-1e-8+x(i);
                
        for j=1:m
            A(j,i)=(deriG{j}(dx(1),dx(2),dx(3))-F(j))/(dx(i)-x(i));
        end
    end
    
    dx=-A\F;
    
    x(1:m)=x(1:m)+dx*0.5;
    check=sum((0.5*dx).^2);
end

for i=1:3
    if x(i)<0;
        x(i)=0;
    end
end

s=x(1);t=x(2);u=x(3);
    
% excess gibbs free energies

Ge_hm = ...
    +(W_hm_il+(3*H-X)*DW_hm_il)*X*(1-H)...
    +(W_hm_gk+(3*H-Y)*DW_hm_gk)*Y*(1-H)...
    +(W_hm_py+(3*H-Z)*DW_hm_py)*Z*(1-H)...
    +(W_hm_cr+(3*H-C)*DW_hm_cr)*C*(1-H)...
    ...
    -1/2*X*Y*(W_il_gk+W_il_gk_T)...
    -1/2*s*t*(W_il_gk-W_il_gk_T)...
    -1/2*X*Z*(W_il_py+W_il_py_T)...
    -1/2*s*u*(W_il_py-W_il_py_T)...
    -1/2*Y*Z*(W_gk_py+W_gk_py_T)...
    -1/2*t*u*(W_gk_py-W_gk_py_T)...
    ...
    -(W_il_cr+2*(X-C)*DW_il_cr)*X*C...
    -(W_gk_cr+2*(Y-C)*DW_gk_cr)*Y*C...
    -(W_py_cr+2*(Z-C)*DW_py_cr)*Z*C...
    ...
    -(W_il_il/4+H*(DW_hm_il/2+W_hm_il2/4)-1/2*C*DW_il_cr)*(X^2-s^2)...
    -(W_gk_gk/4+H*(DW_hm_gk/2+W_hm_gk2/4)-1/2*C*DW_gk_cr)*(Y^2-t^2)...
    -(W_py_py/4+H*(DW_hm_py/2+W_hm_py2/4)-1/2*C*DW_py_cr)*(Z^2-u^2)...
    ...
    +(DW_hm_il/2+W_hm_il2/4)*(1-H)*(X^2-s^2)-1/2*DW_il_cr*C*(X^2-s^2)...
    +(DW_hm_gk/2+W_hm_gk2/4)*(1-H)*(Y^2-t^2)-1/2*DW_gk_cr*C*(Y^2-t^2)...
    +(DW_hm_py/2+W_hm_py2/4)*(1-H)*(Z^2-u^2)-1/2*DW_py_cr*C*(Z^2-u^2)...
    ...
    -3*D2W_il_il/4*(X^2-s^2)*s^2 ...
    -3*D2W_gk_gk/4*(Y^2-t^2)*t^2 ...
    -3*D2W_py_py/4*(Z^2-u^2)*u^2;


Ge_il = Ge_hm...
    + (W_hm_il+(H-X)*DW_hm_il)*(H-X) - 2*DW_hm_il*H*X...
    - (W_hm_gk+(H-Y)*DW_hm_gk)*Y     -   DW_hm_gk*H*Y...
    - (W_hm_py+(H-Z)*DW_hm_py)*Z     -   DW_hm_py*H*Z...
    - (W_hm_cr+(H-C)*DW_hm_cr)*C     -   DW_hm_cr*H*C...
    ...
    + 1/2*Y*(W_il_gk+W_il_gk_T)...
    + 1/2*Z*(W_il_py+W_il_py_T)...
    ...
    + (W_il_cr+(X-C)*DW_il_cr)*C+DW_il_cr*X*C...
    - (DW_hm_il/2+W_hm_il2/4)*(X^2-s^2)...
    + (W_il_il/2+H*(DW_hm_il+W_hm_il2/2)-C*DW_il_cr)*X...
    + D2W_il_il/2*X*s^2;

Ge_gk = Ge_hm...
    + (W_hm_gk+(H-Y)*DW_hm_gk)*(H-Y) - 2*DW_hm_gk*H*Y...
    - (W_hm_il+(H-X)*DW_hm_il)*X     -   DW_hm_il*H*X...
    - (W_hm_py+(H-Z)*DW_hm_py)*Z     -   DW_hm_py*H*Z...
    - (W_hm_cr+(H-C)*DW_hm_cr)*C     -   DW_hm_cr*H*C...
    ...
    + 1/2*X*(W_il_gk+W_il_gk_T)...
    + 1/2*Z*(W_gk_py+W_gk_py_T)...
    ...
    + (W_gk_cr+(Y-C)*DW_gk_cr)*C+DW_gk_cr*Y*C...
    - (DW_hm_gk/2+W_hm_gk2/4)*(Y^2-t^2)...
    + (W_gk_gk/2+H*(DW_hm_gk+W_hm_gk2/2)-C*DW_gk_cr)*Y...
    + D2W_gk_gk/2*Y*t^2;

Ge_py = Ge_hm...
    + (W_hm_py+(H-Z)*DW_hm_py)*(H-Z) - 2*DW_hm_py*H*Z...
    - (W_hm_il+(H-X)*DW_hm_il)*X     -   DW_hm_il*H*X...
    - (W_hm_gk+(H-Y)*DW_hm_gk)*Y     -   DW_hm_gk*H*Y...
    - (W_hm_cr+(H-C)*DW_hm_cr)*C     -   DW_hm_cr*H*C...
    ...
    + 1/2*X*(W_il_py+W_il_py_T)...
    + 1/2*Y*(W_gk_py+W_gk_py_T)...
    ...
    + (W_py_cr+(Z-C)*DW_py_cr)*C+DW_py_cr*Z*C...
    - (DW_hm_py/2+W_hm_py2/4)*(Z^2-u^2)...
    + (W_py_py/2+H*(DW_hm_py+W_hm_py2/2)-C*DW_py_cr)*Z...
    + D2W_py_py/2*Z*u^2;

Ge_cr = Ge_hm...
    - (W_hm_il+(2*H-X)*DW_hm_il)*X...
    - (W_hm_gk+(2*H-Y)*DW_hm_gk)*Y...
    - (W_hm_py+(2*H-Z)*DW_hm_py)*Z...
    - (W_hm_cr+(  H-C)*DW_hm_cr)*(H-C)-2*DW_hm_cr*H*C...
    ...
    + (W_il_cr+(X-2*C)*DW_il_cr)*X...
    + (W_gk_cr+(Y-2*C)*DW_gk_cr)*Y...
    + (W_py_cr+(Z-2*C)*DW_py_cr)*Z...
    - (DW_hm_il/2+DW_il_cr/2+W_hm_il2/4)*(X^2-s^2)...
    - (DW_hm_gk/2+DW_gk_cr/2+W_hm_gk2/4)*(Y^2-t^2)...
    - (DW_hm_py/2+DW_py_cr/2+W_hm_py2/4)*(Z^2-u^2);

% ideal activities

ai_hm = H^(2*(1-alpha));
ai_il = ((X+s)/2)^(1/2)*((X-s)/2)^(1/2)*(X/2)^(-alpha)...
        *((X+Y+Z+s+t+u)/2)^(1/2)*((X+Y+Z-s-t-u)/2)^(1/2)...
        *((X+Y+Z)/2)^(-alpha);
ai_gk = ((Y+t)/2)^(1/2)*((Y-t)/2)^(1/2)*(Y/2)^(-alpha)...
        *((X+Y+Z+s+t+u)/2)^(1/2)*((X+Y+Z-s-t-u)/2)^(1/2)...
        *((X+Y+Z)/2)^(-alpha);
ai_py = ((Z+u)/2)^(1/2)*((Z-u)/2)^(1/2)*(Z/2)^(-alpha)...
        *((X+Y+Z+s+t+u)/2)^(1/2)*((X+Y+Z-s-t-u)/2)^(1/2)...
        *((X+Y+Z)/2)^(-alpha);
ai_cr = C^(2*(1-alpha));

a_hm = ai_hm * exp(Ge_hm/R/T);
a_il = ai_il * exp(Ge_il/R/T);
a_gk = ai_gk * exp(Ge_gk/R/T);
a_py = ai_py * exp(Ge_py/R/T);
a_cr = ai_cr * exp(Ge_cr/R/T);

a_temp=[a_hm,a_il,a_gk,a_py,a_cr]';
p_temp=[   H,   X,   Y,   Z,   C]';
s_temp=[   0,   s,   t,   t,   0]';   

k=1;
for i=1:5
    if p_temp(i)>0
        p(k)=p_temp(i);
        a(k)=a_temp(i)+1i*s_temp(i);
        k=k+1;
    end
end

end


%newly added carbonate start
function act=CalcA_cab(T,x)
a=x(1);b=x(2);c=x(3);d=x(4);e=x(5);
X_mg_dol = a/e ;
X_fe_dol = b/e ;
X_ca_dol = 1-a/e-b/e;
X_mg_cc = c/(1-e) ;
X_fe_cc = d/(1-e) ;
X_ca_cc = 1-c/(1-e)-d/(1-e);
%refine input
R=8.31e-3;

% coefficient
W_camg_cc = 23.24;
W_camg_dol = -96.85+T*0.03623;
W_mgca_cc = 24.3-T*0.007743;
W_mgca_dol = -55.48-T*0.02285;
W_cafe_cc = 27.3132-5796.51/T-(1.43147e-6)*(T^2); 
W_cafe_dol = -1.04-T*0.2124+(2.027e-4)*(T^2);
W_feca_cc = -70.7514+0.0918848*T+27924.4/T-(3.59552e-5)*(T^2);
W_feca_dol = -86.74+T*0.07118+(9.184e-5)*(T^2);
W_mgfe_cc = 0;
W_femg_cc = 0;
W_mgfe_dol = 0;
W_femg_dol = 0;
C_cc = ((W_mgca_cc-W_camg_cc)+(W_cafe_cc-W_feca_cc)+(W_femg_cc-W_mgfe_cc))/2;
C_dol = -185.45-T*0.3692+(1.418e-4)*(T^2);

% standard miu
miu_dol_mg = 0;
miu_dol_fe = 0;
miu_dol_ca = 0;
miu_cc_mg = 0;
miu_cc_fe = 0;
miu_cc_ca = 0;


% acitivity
a_dol_mg = exp(((R*T*log(X_mg_dol)+((X_fe_dol)^2*(W_mgfe_dol+2*X_mg_dol*(W_femg_dol-W_mgfe_dol))...
+ (X_ca_dol)^2*(W_mgca_dol+2*X_mg_dol*(W_camg_dol-W_mgca_dol))...
+ (X_fe_dol*X_ca_dol)*((W_femg_dol+W_mgfe_dol+W_camg_dol+W_mgca_dol-W_feca_dol-W_cafe_dol)/2 + ...
 X_mg_dol*(W_femg_dol-W_mgfe_dol+W_camg_dol-W_mgca_dol)+(X_fe_dol-X_ca_dol)*(W_feca_dol-W_cafe_dol)...
- (1-2*X_mg_dol)*C_dol)))-miu_dol_mg)/R/T);

a_fe_dol = exp(((R*T*log(X_fe_dol)+((X_ca_dol)^2*(W_feca_dol+2*X_fe_dol*(W_cafe_dol-W_feca_dol))...
+ (X_mg_dol)^2*(W_femg_dol+2*X_fe_dol*(W_mgfe_dol-W_femg_dol))...
+ (X_ca_dol*X_mg_dol)*((W_cafe_dol+W_feca_dol+W_mgfe_dol+W_femg_dol-W_camg_dol-W_mgca_dol)/2 + ...
 X_fe_dol*(W_cafe_dol-W_feca_dol+W_mgfe_dol-W_femg_dol)+(X_ca_dol-X_mg_dol)*(W_camg_dol-W_mgca_dol)...
- (1-2*X_fe_dol)*C_dol)))-miu_dol_fe)/R/T);



a_ca_dol = exp(((R*T*log(X_ca_dol)+((X_mg_dol)^2*(W_camg_dol+2*X_ca_dol*(W_mgca_dol-W_camg_dol))...
+ (X_fe_dol)^2*(W_cafe_dol+2*X_ca_dol*(W_feca_dol-W_cafe_dol))...
+ (X_mg_dol*X_fe_dol)*((W_mgca_dol+W_camg_dol+W_feca_dol+W_cafe_dol-W_mgfe_dol-W_femg_dol)/2 + ...
 X_ca_dol*(W_mgca_dol-W_camg_dol+W_feca_dol-W_cafe_dol)+(X_mg_dol-X_fe_dol)*(W_mgfe_dol-W_femg_dol)...
- (1-2*X_ca_dol)*C_dol)))-miu_dol_ca)/R/T);



a_mg_cc = exp(((R*T*log(X_mg_cc)+((X_fe_cc)^2*(W_mgfe_cc+2*X_mg_cc*(W_femg_cc-W_mgfe_cc))...
+ (X_ca_cc)^2*(W_mgca_cc+2*X_mg_cc*(W_camg_cc-W_mgca_cc))...
+ (X_fe_cc*X_ca_cc)*((W_femg_cc+W_mgfe_cc+W_camg_cc+W_mgca_cc-W_feca_cc-W_cafe_cc)/2 + ...
 X_mg_cc*(W_femg_cc-W_mgfe_cc+W_camg_cc-W_mgca_cc)+(X_fe_cc-X_ca_cc)*(W_feca_cc-W_cafe_cc)...
- (1-2*X_mg_cc)*C_cc)))-miu_cc_mg)/R/T);


a_fe_cc = exp(((R*T*log(X_fe_cc) + ((X_ca_cc)^2*(W_feca_cc+2*X_fe_cc*(W_cafe_cc-W_feca_cc))...
+ (X_mg_cc)^2*(W_femg_cc+2*X_fe_cc*(W_mgfe_cc-W_femg_cc))...
+ (X_ca_cc*X_mg_cc)*((W_cafe_cc+W_feca_cc+W_mgfe_cc+W_femg_cc-W_camg_cc-W_mgca_cc)/2 + ...
 X_fe_cc*(W_cafe_cc-W_feca_cc+W_mgfe_cc-W_femg_cc)+(X_ca_cc-X_mg_cc)*(W_camg_cc-W_mgca_cc)...
- (1-2*X_fe_cc)*C_cc)))-miu_cc_fe)/R/T);


a_ca_cc = exp(((R*T*log(X_ca_cc) + ((X_mg_cc)^2*(W_camg_cc+2*X_ca_cc*(W_mgca_cc-W_camg_cc))...
+ (X_fe_cc)^2*(W_cafe_cc+2*X_ca_cc*(W_feca_cc-W_cafe_cc))...
+ (X_mg_cc*X_fe_cc)*((W_mgca_cc+W_camg_cc+W_feca_cc+W_cafe_cc-W_mgfe_cc-W_femg_cc)/2 +...
 X_ca_cc*(W_mgca_cc-W_camg_cc+W_feca_cc-W_cafe_cc)+(X_mg_cc-X_fe_cc)*(W_mgfe_cc-W_femg_cc)...
- (1-2*X_ca_cc)*C_cc)))-miu_cc_ca)/R/T);

act = [a_dol_mg,a_fe_dol,a_ca_dol,a_mg_cc,a_fe_cc,a_ca_cc];
end
%newly added carbonate end