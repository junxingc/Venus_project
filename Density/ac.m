load('HP98.mat')
% make a(241) = 3/7 cumm (66) + 4/7 grun(67)
HP98(241).minrl='a';
HP98(241).Sym='MAKE';
HP98(241).H     = 3/7*HP98(66).H + 4/7*HP98(67).H;
HP98(241).S     = 3/7*HP98(66).S + 4/7*HP98(67).S;
HP98(241).V     = 3/7*HP98(66).V + 4/7*HP98(67).V;
HP98(241).a     = 3/7*HP98(66).a + 4/7*HP98(67).a;
HP98(241).b     = 3/7*HP98(66).b + 4/7*HP98(67).b;
HP98(241).c     = 3/7*HP98(66).c + 4/7*HP98(67).c;
HP98(241).d     = 3/7*HP98(66).d + 4/7*HP98(67).d;
HP98(241).alpha0= 1/HP98(241).V * (3/7*HP98(66).alpha0*HP98(66).V + 4/7*HP98(67).alpha0*HP98(67).V);
HP98(241).alpha1= 10;
HP98(241).K     =   HP98(241).V / (3/7*HP98(66).V     /HP98(66).K + 4/7*HP98(67).V     /HP98(67).K);
HP98(241).dKdp  = 4  ;
HP98(241).dKdT  = -1.5e-4*HP98(241).K  ;
HP98(241).Tc    = NaN;
HP98(241).Smax  = NaN;
HP98(241).Vmax  = NaN;
HP98(241).Comp  = 3/7*HP98(66).Comp + 4/7*HP98(67).Comp;

% make b(242) = 2/7 cumm (66) + 5/7 grun(67)
HP98(242).minrl='b';
HP98(242).Sym='MAKE';
HP98(242).H     = 2/7*HP98(66).H + 5/7*HP98(67).H;
HP98(242).S     = 2/7*HP98(66).S + 5/7*HP98(67).S;
HP98(242).V     = 2/7*HP98(66).V + 5/7*HP98(67).V;
HP98(242).a     = 2/7*HP98(66).a + 5/7*HP98(67).a;
HP98(242).b     = 2/7*HP98(66).b + 5/7*HP98(67).b;
HP98(242).c     = 2/7*HP98(66).c + 5/7*HP98(67).c;
HP98(242).d     = 2/7*HP98(66).d + 5/7*HP98(67).d;
HP98(242).alpha0= 1/HP98(242).V * (2/7*HP98(66).alpha0*HP98(66).V + 5/7*HP98(67).alpha0*HP98(67).V);
HP98(242).alpha1= 10;
HP98(242).K     =   HP98(242).V / (2/7*HP98(66).V     /HP98(66).K + 5/7*HP98(67).V     /HP98(67).K);
HP98(242).dKdp  = 4  ;
HP98(242).dKdT  = -1.5e-4*HP98(242).K  ;
HP98(242).Tc    = NaN;
HP98(242).Smax  = NaN;
HP98(242).Vmax  = NaN;
HP98(242).Comp  = 2/7*HP98(66).Comp + 5/7*HP98(67).Comp;

% make mrb(243) = 1 gl(61) - 2 jd(50) + 2 acm(51)
HP98(243).minrl='mrb';
HP98(243).Sym='MAKE';
HP98(243).H     = HP98(61).H - 2*HP98(50).H + 2*HP98(51).H;
HP98(243).S     = HP98(61).S - 2*HP98(50).S + 2*HP98(51).S;
HP98(243).V     = HP98(61).V - 2*HP98(50).V + 2*HP98(51).V;
HP98(243).a     = HP98(61).a - 2*HP98(50).a + 2*HP98(51).a;
HP98(243).b     = HP98(61).b - 2*HP98(50).b + 2*HP98(51).b;
HP98(243).c     = HP98(61).c - 2*HP98(50).c + 2*HP98(51).c;
HP98(243).d     = HP98(61).d - 2*HP98(50).d + 2*HP98(51).d;
HP98(243).alpha0= 1/HP98(243).V * (HP98(61).alpha0*HP98(61).V - 2*HP98(50).alpha0*HP98(50).V + 2*HP98(51).alpha0*HP98(51).V);
HP98(243).alpha1= 10;
HP98(243).K     =   HP98(243).V / (HP98(61).V     /HP98(61).K - 2*HP98(50).V     /HP98(50).K + 2*HP98(51).V     /HP98(51).K);
HP98(243).dKdp  = 4;
HP98(243).dKdT  = -1.5e-4*HP98(243).K  ;
HP98(243).Tc    = NaN;
HP98(243).Smax  = NaN;
HP98(243).Vmax  = NaN;
HP98(243).Comp  = HP98(61).Comp - 2*HP98(50).Comp + 2*HP98(51).Comp;

% make om = 1/2 jd (50) + 1/2 di(48)
n=231;n1=50;n2=48;
HP98(n).minrl='om';
HP98(n).Sym='MAKE';
HP98(n).H     = 1/2*HP98(n1).H + 1/2*HP98(n2).H;
HP98(n).S     = 1/2*HP98(n1).S + 1/2*HP98(n2).S;
HP98(n).V     = 1/2*HP98(n1).V + 1/2*HP98(n2).V;
HP98(n).a     = 1/2*HP98(n1).a + 1/2*HP98(n2).a;
HP98(n).b     = 1/2*HP98(n1).b + 1/2*HP98(n2).b;
HP98(n).c     = 1/2*HP98(n1).c + 1/2*HP98(n2).c;
HP98(n).d     = 1/2*HP98(n1).d + 1/2*HP98(n2).d;
HP98(n).alpha0= 1/HP98(n).V * (1/2*HP98(n1).alpha0*HP98(n1).V + 1/2*HP98(n2).alpha0*HP98(n2).V);
HP98(n).alpha1= 10;
HP98(n).K     =   HP98(n).V / (1/2*HP98(n1).V     /HP98(n1).K + 1/2*HP98(n2).V     /HP98(n2).K);
HP98(n).dKdp  = 4  ;
HP98(n).dKdT  = -1.5e-4*HP98(n).K  ;
HP98(n).Tc    = NaN;
HP98(n).Smax  = NaN;
HP98(n).Vmax  = NaN;
HP98(n).Comp  = 1/2*HP98(n1).Comp + 1/2*HP98(n2).Comp;

% make cfm = 1/2 di(48) + 1/2 hed(49)
n=232;n1=48;n2=49;
HP98(n).minrl='cfm';
HP98(n).Sym='MAKE';
HP98(n).H     = 1/2*HP98(n1).H + 1/2*HP98(n2).H;
HP98(n).S     = 1/2*HP98(n1).S + 1/2*HP98(n2).S;
HP98(n).V     = 1/2*HP98(n1).V + 1/2*HP98(n2).V;
HP98(n).a     = 1/2*HP98(n1).a + 1/2*HP98(n2).a;
HP98(n).b     = 1/2*HP98(n1).b + 1/2*HP98(n2).b;
HP98(n).c     = 1/2*HP98(n1).c + 1/2*HP98(n2).c;
HP98(n).d     = 1/2*HP98(n1).d + 1/2*HP98(n2).d;
HP98(n).alpha0= 1/HP98(n).V * (1/2*HP98(n1).alpha0*HP98(n1).V + 1/2*HP98(n2).alpha0*HP98(n2).V);
HP98(n).alpha1= 10;
HP98(n).K     =   HP98(n).V / (1/2*HP98(n1).V     /HP98(n1).K + 1/2*HP98(n2).V     /HP98(n2).K);
HP98(n).dKdp  = 4  ;
HP98(n).dKdT  = -1.5e-4*HP98(n).K  ;
HP98(n).Tc    = NaN;
HP98(n).Smax  = NaN;
HP98(n).Vmax  = NaN;
HP98(n).Comp  = 1/2*HP98(n1).Comp + 1/2*HP98(n2).Comp;

% make jac = 1/2 jd(50) + 1/2 acm(51)
n=233;n1=50;n2=51;
HP98(n).minrl='jac';
HP98(n).Sym='MAKE';
HP98(n).H     = 1/2*HP98(n1).H + 1/2*HP98(n2).H;
HP98(n).S     = 1/2*HP98(n1).S + 1/2*HP98(n2).S;
HP98(n).V     = 1/2*HP98(n1).V + 1/2*HP98(n2).V;
HP98(n).a     = 1/2*HP98(n1).a + 1/2*HP98(n2).a;
HP98(n).b     = 1/2*HP98(n1).b + 1/2*HP98(n2).b;
HP98(n).c     = 1/2*HP98(n1).c + 1/2*HP98(n2).c;
HP98(n).d     = 1/2*HP98(n1).d + 1/2*HP98(n2).d;
HP98(n).alpha0= 1/HP98(n).V * (1/2*HP98(n1).alpha0*HP98(n1).V + 1/2*HP98(n2).alpha0*HP98(n2).V);
HP98(n).alpha1= 10;
HP98(n).K     =   HP98(n).V / (1/2*HP98(n1).V     /HP98(n1).K + 1/2*HP98(n2).V     /HP98(n2).K);
HP98(n).dKdp  = 4  ;
HP98(n).dKdT  = -1.5e-4*HP98(n).K  ;
HP98(n).Tc    = NaN;
HP98(n).Smax  = NaN;
HP98(n).Vmax  = NaN;
HP98(n).Comp  = 1/2*HP98(n1).Comp + 1/2*HP98(n2).Comp;

% make obi = 2/3 phl + 1/3 ann
HP98(201).minrl='obi';
HP98(201).Sym='MAKE';
HP98(201).H     = 2/3*HP98(80).H  +1/3*HP98(81).H;
HP98(201).S     = 2/3*HP98(80).S  +1/3*HP98(81).S;
HP98(201).V     = 2/3*HP98(80).V  +1/3*HP98(81).V;
HP98(201).a     = 2/3*HP98(80).a  +1/3*HP98(81).a;
HP98(201).b     = 2/3*HP98(80).b  +1/3*HP98(81).b;
HP98(201).c     = 2/3*HP98(80).c  +1/3*HP98(81).c;
HP98(201).d     = 2/3*HP98(80).d  +1/3*HP98(81).d;
HP98(201).alpha0=     HP98(80).alpha0;
HP98(201).alpha1=     HP98(80).alpha1;
HP98(201).K     =     HP98(80).K     ;
HP98(201).dKdp  =     HP98(80).dKdp  ;
HP98(201).dKdT  =     HP98(80).dKdT  ;
HP98(201).Tc    =NaN;
HP98(201).Smax  =NaN;
HP98(201).Vmax  =NaN;
HP98(201).Comp  = 2/3*HP98(80).Comp+1/3*HP98(81).Comp;

% make tbi = phl(80) - br(135) + ru(120)
HP98(202).minrl='tbi';
HP98(202).Sym='MAKE';
HP98(202).H     = HP98(80).H - HP98(135).H + HP98(120).H;
HP98(202).S     = HP98(80).S - HP98(135).S + HP98(120).S;
HP98(202).V     = HP98(80).V - HP98(135).V + HP98(120).V;
HP98(202).a     = HP98(80).a - HP98(135).a + HP98(120).a;
HP98(202).b     = HP98(80).b - HP98(135).b + HP98(120).b;
HP98(202).c     = HP98(80).c - HP98(135).c + HP98(120).c;
HP98(202).d     = HP98(80).d - HP98(135).d + HP98(120).d;
HP98(202).alpha0= 1/HP98(202).V * (HP98(80).alpha0*HP98(80).V - HP98(135).alpha0*HP98(135).V + HP98(120).alpha0*HP98(120).V);
HP98(202).alpha1= 10;
HP98(202).K     = HP98(202).V   / (HP98(80).V/HP98(80).K      - HP98(135).V/HP98(135).K      + HP98(120).V/HP98(120).K);
HP98(202).dKdp  = 4  ;
HP98(202).dKdT  = -1.5e-4*HP98(202).K  ;
HP98(202).Tc    = NaN;
HP98(202).Smax  = NaN;
HP98(202).Vmax  = NaN;
HP98(202).Comp  = HP98(80).Comp - HP98(135).Comp + HP98(120).Comp;

% make fbi = east(83) - 1/2 cor(123) + 1/2 hem(124)
HP98(203).minrl='fbi';
HP98(203).Sym='MAKE';
HP98(203).H     = -5895.64;
HP98(203).S     = HP98(83).S - 1/2* HP98(123).S + 1/2* HP98(124).S;
HP98(203).V     = HP98(83).V - 1/2* HP98(123).V + 1/2* HP98(124).V;
HP98(203).a     = HP98(83).a - 1/2* HP98(123).a + 1/2* HP98(124).a;
HP98(203).b     = HP98(83).b - 1/2* HP98(123).b + 1/2* HP98(124).b;
HP98(203).c     = HP98(83).c - 1/2* HP98(123).c + 1/2* HP98(124).c;
HP98(203).d     = -HP98(83).d - 1/2* HP98(123).d + 1/2* HP98(124).d;
HP98(203).alpha0= 1/HP98(203).V * (HP98(83).alpha0*HP98(83).V - 1/2* HP98(123).alpha0*HP98(123).V + 1/2* HP98(124).alpha0*HP98(124).V);
HP98(203).alpha1= 10;
HP98(203).K     = HP98(203).V / (HP98(83).V/HP98(83).K      - 1/2* HP98(123).V/HP98(123).K      + 1/2* HP98(124).V/HP98(124).K);
HP98(203).dKdp  = 4  ;
HP98(203).dKdT  = -0.0769  ;
HP98(203).Tc    = NaN;
HP98(203).Smax  = NaN;
HP98(203).Vmax  = NaN;
HP98(203).Comp  = HP98(83).Comp - 1/2* HP98(123).Comp + 1/2* HP98(124).Comp;

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


amount_1 = length(HP98);
for i=1:amount_1
    if length(HP98(i).Comp) ~= 0;
        sum_mass = HP98(i).Comp(1,1)*60.084 + HP98(i).Comp(1,2)*79.866 + HP98(i).Comp(1,3)*101.96 + HP98(i).Comp(1,4)*71.844...
         + HP98(i).Comp(1,5)*40.3044 + HP98(i).Comp(1,6)*86.9368 + HP98(i).Comp(1,7)*56.0774 + HP98(i).Comp(1,8)* 61.9789...
         + HP98(i).Comp(1,9)*94.2 + HP98(i).Comp(1,10)*15.999 + HP98(i).Comp(1,11)*18.01528 + HP98(i).Comp(1,12)*44.01;
        HP98(i).mass=sum_mass;
    end
end

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


% 57:tr 59:ts 60:parg 61:gl 66:cumm 67:grun 241:a 242:b 243:mrb
hb.minrl='hb';
hb.sys='NCFMASHO';
hb.endmember=[57,59,60,61,66,67,241,242,243];
hb.endmass=[HP98(57).mass,HP98(59).mass,HP98(60).mass,HP98(61).mass,HP98(66).mass,HP98(67).mass,HP98(241).mass,HP98(242).mass,HP98(243).mass];
hb.endnum=9;
hb.ensym={'tr','ts','parg','gl','cumm','grun','a','b','mrb'};
hb.variables={'x(hb)','y(hb)','z(hb)','a(hb)','c(hb)','f(hb)','Q1(hb)','Q2(hb)'};
hb.p={@(X)X(5)+X(3)-X(2)-X(6)-X(4)/2;
               @(X)    -X(3)+X(2)+X(6)-X(4)/2;
               @(X)                    X(4);
               @(X)X(3)-X(6);
               @(X)(1-X(1))*(1-X(5)-X(3)          ) -   X(8)*(1-X(2)-X(6)) - 3/2*X(7);
               @(X)   X(1) *(1+X(5)+X(3)-X(2)-X(6)) - 2*X(8)*(1-X(2)-X(6)) - 5/2*X(7);
               @(X)  -X(1) *(  X(5)+X(3)          ) +   X(8)*(1-X(2)-X(6)) + 5/2*X(7);
               @(X)  -X(1) *(  X(5)+X(3)-X(2)-X(6)) + 2*X(8)*(1-X(2)-X(6)) + 3/2*X(7);
               @(X)                           X(6)}'; 
hb.model='asf';
hb.W=[   20,   25,   65,   45,   75,   57,   63,   65;
                 NaN,  -40,   25,   70,   80,   70, 72.5,   25;
                 NaN,  NaN,   50,   90,106.7, 94.8, 94.8,   50;
                 NaN,  NaN,  NaN,  100,113.5,  100,111.2,    0;
                 NaN,  NaN,  NaN,  NaN,   33,   18,   23,  100;
                 NaN,  NaN,  NaN,  NaN,  NaN,   12,    8,113.5;
                 NaN,  NaN,  NaN,  NaN,  NaN,  NaN,   20,  100;
                 NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,111.2];
hb.dWdT=[  0,  0,  0,  0,  0,  0,  0,  0;
                  NaN,  0,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,NaN,  0];
hb.dWdP=[  0,  0,  0,  0,  0,  0,  0,  0;
                  NaN,  0,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,NaN,  0];
hb.a_ideal={@(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))^2    *(            X(5)                                      )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)2*(1-X(4)) * (1-X(1)+X(7))^3 * (                 X(2)      )^2    *(            X(5)                                      )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2)^(1/2)*(  X(4)/4+X(6)/2+X(2)/2-X(3)/2)^(1/2);
                     @(X)8*(  X(4)) * (1-X(1)+X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))*X(2) *(            X(5)                                      )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2)^(1/2)*(  X(4)/4+X(6)/2+X(2)/2-X(3)/2)^(1/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * (                 X(2)      )^2    *(                 X(3)                                 )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))^2    *((1-X(1))*(1-X(5)-X(3)) - X(8)*(1-X(2)-X(6)) - 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (  X(1)-X(7))^3 * ((  X(1)-X(8))*(1-X(2)-X(6)))^2    *(   X(1) *(1-X(5)-X(3)) + X(8)*(1-X(2)-X(6)) + 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * ((  X(1)-X(8))*(1-X(2)-X(6)))^2    *(   X(1) *(1-X(5)-X(3)) + X(8)*(1-X(2)-X(6)) + 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (  X(1)-X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))^2    *(   X(1) *(1-X(5)-X(3)) + X(8)*(1-X(2)-X(6)) + 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * (                      X(6) )^2    *(                 X(3)                                 )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2)}';
hb.alpha=[  1,1.5,1.7,0.8,  1,  1,  1,  1,0.8];
hb.dalphadT=[0,0,0,0,0,0,0,0,0];
hb.dalphadP=[0,0,0,0,0,0,0,0,0];
hb.DQF=[0,10,15,3,-6.4,-5,-15.1,-17.1,8];
hb.dDQFdT=[0,0,0,0,0,0,0,0,0];
hb.dDQFdP=[0,0,0,0,0,0,0,0,0];
           
gl.minrl='gl';
gl.sys='NCFMASHO';
gl.endmember=[57,59,60,61,66,67,241,242,243];
gl.endmass=[HP98(57).mass,HP98(59).mass,HP98(60).mass,HP98(61).mass,HP98(66).mass,HP98(67).mass,HP98(241).mass,HP98(242).mass,HP98(243).mass];
gl.endnum=9;
gl.ensym={'tr','ts','parg','gl','cumm','grun','a','b','mrb'};
gl.variables={'x(gl)','y(gl)','z(gl)','a(gl)','c(gl)','f(gl)','Q1(gl)','Q2(gl)'};
gl.p={@(X)X(5)+X(3)-X(2)-X(6)-X(4)/2;
               @(X)    -X(3)+X(2)+X(6)-X(4)/2;
               @(X)                    X(4);
               @(X)X(3)-X(6);
               @(X)(1-X(1))*(1-X(5)-X(3)          ) -   X(8)*(1-X(2)-X(6)) - 3/2*X(7);
               @(X)   X(1) *(1+X(5)+X(3)-X(2)-X(6)) - 2*X(8)*(1-X(2)-X(6)) - 5/2*X(7);
               @(X)  -X(1) *(  X(5)+X(3)          ) +   X(8)*(1-X(2)-X(6)) + 5/2*X(7);
               @(X)  -X(1) *(  X(5)+X(3)-X(2)-X(6)) + 2*X(8)*(1-X(2)-X(6)) + 3/2*X(7);
               @(X)               X(6)}'; 
gl.a_ideal={@(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))^2    *(            X(5)                                      )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)2*(1-X(4)) * (1-X(1)+X(7))^3 * (                 X(2)      )^2    *(            X(5)                                      )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2)^(1/2)*(  X(4)/4+X(6)/2+X(2)/2-X(3)/2)^(1/2);
                     @(X)8*(  X(4)) * (1-X(1)+X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))*X(2) *(            X(5)                                      )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2)^(1/2)*(  X(4)/4+X(6)/2+X(2)/2-X(3)/2)^(1/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * (                 X(2)      )^2    *(                 X(3)                                 )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))^2    *((1-X(1))*(1-X(5)-X(3)) - X(8)*(1-X(2)-X(6)) - 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (  X(1)-X(7))^3 * ((  X(1)-X(8))*(1-X(2)-X(6)))^2    *(   X(1) *(1-X(5)-X(3)) + X(8)*(1-X(2)-X(6)) + 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * ((  X(1)-X(8))*(1-X(2)-X(6)))^2    *(   X(1) *(1-X(5)-X(3)) + X(8)*(1-X(2)-X(6)) + 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (  X(1)-X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))^2    *(   X(1) *(1-X(5)-X(3)) + X(8)*(1-X(2)-X(6)) + 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * (                      X(6) )^2    *(                 X(3)                                 )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2)}';
gl.model='asf';
gl.W=[   20,   25,   65,   45,   75,   57,   63,   65;
                 NaN,  -40,   25,   70,   80,   70, 72.5,   25;
                 NaN,  NaN,   50,   90,106.7, 94.8, 94.8,   50;
                 NaN,  NaN,  NaN,  100,113.5,  100,111.2,    0;
                 NaN,  NaN,  NaN,  NaN,   33,   18,   23,  100;
                 NaN,  NaN,  NaN,  NaN,  NaN,   12,    8,113.5;
                 NaN,  NaN,  NaN,  NaN,  NaN,  NaN,   20,  100;
                 NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,111.2];
gl.dWdT=[  0,  0,  0,  0,  0,  0,  0,  0;
                  NaN,  0,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,NaN,  0];
gl.dWdP=[  0,  0,  0,  0,  0,  0,  0,  0;
                  NaN,  0,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,NaN,  0];
gl.alpha=[  1,1.5,1.7,0.8,  1,  1,  1,  1,0.8];
gl.dalphadT=[0,0,0,0,0,0,0,0,0];
gl.dalphadP=[0,0,0,0,0,0,0,0,0];
gl.DQF=[0,10,15,3,-6.4,-5,-15.1,-17.1,8];
gl.dDQFdT=[0,0,0,0,0,0,0,0,0];
gl.dDQFdP=[0,0,0,0,0,0,0,0,0];
           
           
act.minrl='act';
act.sys='NCFMASHO';
act.endmember=[57,59,60,61,66,67,241,242,243];
act.endmass=[HP98(57).mass,HP98(59).mass,HP98(60).mass,HP98(61).mass,HP98(66).mass,HP98(67).mass,HP98(241).mass,HP98(242).mass,HP98(243).mass];
act.endnum=9;
act.ensym={'tr','ts','parg','gl','cumm','grun','a','b','mrb'};
act.variables={'x(act)','y(act)','z(act)','a(act)','c(act)','f(act)','Q1(act)','Q2(act)'};
act.p={@(X)X(5)+X(3)-X(2)-X(6)-X(4)/2;
               @(X)    -X(3)+X(2)+X(6)-X(4)/2;
               @(X)                    X(4);
               @(X)X(3)-X(6);
               @(X)(1-X(1))*(1-X(5)-X(3)          ) -   X(8)*(1-X(2)-X(6)) - 3/2*X(7);
               @(X)   X(1) *(1+X(5)+X(3)-X(2)-X(6)) - 2*X(8)*(1-X(2)-X(6)) - 5/2*X(7);
               @(X)  -X(1) *(  X(5)+X(3)          ) +   X(8)*(1-X(2)-X(6)) + 5/2*X(7);
               @(X)  -X(1) *(  X(5)+X(3)-X(2)-X(6)) + 2*X(8)*(1-X(2)-X(6)) + 3/2*X(7);
               @(X)               X(6)}';
act.a_ideal={@(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))^2    *(            X(5)                                      )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)2*(1-X(4)) * (1-X(1)+X(7))^3 * (                 X(2)      )^2    *(            X(5)                                      )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2)^(1/2)*(  X(4)/4+X(6)/2+X(2)/2-X(3)/2)^(1/2);
                     @(X)8*(  X(4)) * (1-X(1)+X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))*X(2) *(            X(5)                                      )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2)^(1/2)*(  X(4)/4+X(6)/2+X(2)/2-X(3)/2)^(1/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * (                 X(2)      )^2    *(                 X(3)                                 )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))^2    *((1-X(1))*(1-X(5)-X(3)) - X(8)*(1-X(2)-X(6)) - 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (  X(1)-X(7))^3 * ((  X(1)-X(8))*(1-X(2)-X(6)))^2    *(   X(1) *(1-X(5)-X(3)) + X(8)*(1-X(2)-X(6)) + 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * ((  X(1)-X(8))*(1-X(2)-X(6)))^2    *(   X(1) *(1-X(5)-X(3)) + X(8)*(1-X(2)-X(6)) + 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (  X(1)-X(7))^3 * ((1-X(1)+X(8))*(1-X(2)-X(6)))^2    *(   X(1) *(1-X(5)-X(3)) + X(8)*(1-X(2)-X(6)) + 3/2*X(7))^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2);
                     @(X)1*(1-X(4)) * (1-X(1)+X(7))^3 * (                      X(6) )^2    *(                 X(3)                                 )^2 *(1-X(4)/4-X(6)/2-X(2)/2+X(3)/2)}';
act.model='asf';
act.W=[   20,   25,   65,   45,   75,   57,   63,   65;
                 NaN,  -40,   25,   70,   80,   70, 72.5,   25;
                 NaN,  NaN,   50,   90,106.7, 94.8, 94.8,   50;
                 NaN,  NaN,  NaN,  100,113.5,  100,111.2,    0;
                 NaN,  NaN,  NaN,  NaN,   33,   18,   23,  100;
                 NaN,  NaN,  NaN,  NaN,  NaN,   12,    8,113.5;
                 NaN,  NaN,  NaN,  NaN,  NaN,  NaN,   20,  100;
                 NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,111.2];
act.dWdT=[  0,  0,  0,  0,  0,  0,  0,  0;
                  NaN,  0,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,NaN,  0];
act.dWdP=[  0,  0,  0,  0,  0,  0,  0,  0;
                  NaN,  0,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,  0,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,  0,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,  0,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,  0,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,  0,  0;
                  NaN,NaN,NaN,NaN,NaN,NaN,NaN,  0];
act.alpha=[  1,1.5,1.7,0.8,  1,  1,  1,  1,0.8];
act.dalphadT=[0,0,0,0,0,0,0,0,0];
act.dalphadP=[0,0,0,0,0,0,0,0,0];
act.DQF=[0,10,15,3,-6.4,-5,-15.1,-17.1,8];
act.dDQFdT=[0,0,0,0,0,0,0,0,0];
act.dDQFdP=[0,0,0,0,0,0,0,0,0];

% -----------------------------------------------------------------------

% 50:jd;48:di;hed:49;acm:51;om:231;cfm:232;jac:233
di.minrl='di';
di.sys='NCFMASO';
di.endmember=[50,48,49,51,231,232,233];
di.endmass=[HP98(50).mass,HP98(48).mass,HP98(49).mass,HP98(51).mass,HP98(231).mass,HP98(232).mass,HP98(233).mass];
di.endnum=7;
di.ensym={'jd','di','hed','acm','om','cfm','jac'};
di.variables={'x(di)';'j(di)';'f(di)';'Q(di)';'Qa(di)';'Qf(di)'};
di.p={@(X)X(2)*(1-X(3))-X(4)-X(5);                      % jd
              @(X)(1-X(2)-X(4))*(1-X(1)+X(6))-2*X(4)*X(1);      % di
              @(X)(1-X(2)-X(4))*(  X(1)+X(6));                  % hed
              @(X)X(3)*X(2)-X(5);                               % acm
              @(X)2*X(4);                                       % om
              @(X)2*X(4)*X(1)+2*X(6)*(X(2)-1+X(4));             % cfm
              @(X)2*X(5)}';                                     % jac
di.a_ideal={@(X)(                X(2)*(1-X(3))-X(4)+X(5))^(1/2) * (    X(2)*(1-X(3))+X(4)-X(5))^(1/2) * (  X(2)-X(4))^(1/2) * (  X(2)+X(4))^(1/2);
                    @(X)((1-X(2)+X(4))*(1-X(1)+X(6))-2*X(4)*X(6))^(1/2) * ((1-X(2)-X(4))*(1-X(1)-X(6)))^(1/2) * (1-X(2)+X(4))^(1/2) * (1-X(2)-X(4))^(1/2);
                    @(X)((1-X(2)+X(4))*(  X(1)-X(6))+2*X(4)*X(6))^(1/2) * ((1-X(2)-X(4))*(  X(1)+X(6)))^(1/2) * (1-X(2)+X(4))^(1/2) * (1-X(2)-X(4))^(1/2);
                    @(X)(                X(2)*   X(3)      -X(5))^(1/2) * (    X(2)*   X(3)      +X(5))^(1/2) * (  X(2)-X(4))^(1/2) * (  X(2)+X(4))^(1/2);
                    @(X)((1-X(2)+X(4))*(1-X(1)+X(6))-2*X(4)*X(6))^(1/2) * (    X(2)*(1-X(3))+X(4)-X(5))^(1/2) * (1-X(2)+X(4))^(1/2) * (  X(2)+X(4))^(1/2);
                    @(X)((1-X(2)+X(4))*(  X(1)-X(6))+2*X(4)*X(6))^(1/2) * ((1-X(2)-X(4))*(1-X(1)-X(6)))^(1/2) * (1-X(2)+X(4))^(1/2) * (1-X(2)-X(4))^(1/2);
                    @(X)(                X(2)*(1-X(3))-X(4)+X(5))^(1/2) * (    X(2)*   X(3)      +X(5))^(1/2) * (  X(2)-X(4))^(1/2) * (  X(2)+X(4))^(1/2)}';
di.model='sf';
di.W=[   26,   24,    5, 15.5, 25.2,    3;
                NaN,    4,   15,15.75,    2,21.05;
                NaN,  NaN,   14, 17.2,    2, 20.1;
                NaN,  NaN,  NaN, 12.8, 15.5,    3;
                NaN,  NaN,  NaN,  NaN,18.45, 19.3;
                NaN,  NaN,  NaN,  NaN,  NaN,21.05];
di.dWdT=[  0,  0,  0,  0,  0,  0;
                 NaN,  0,  0,  0,  0,  0;
                 NaN,NaN,  0,  0,  0,  0;
                 NaN,NaN,NaN,  0,  0,  0;
                 NaN,NaN,NaN,NaN,  0,  0;
                 NaN,NaN,NaN,NaN,NaN,  0];
di.dWdP=[  0,  0,  0,  0,  0,  0;
                 NaN,  0,  0,  0,  0,  0;
                 NaN,NaN,  0,  0,  0,  0;
                 NaN,NaN,NaN,  0,  0,  0;
                 NaN,NaN,NaN,NaN,  0,  0;
                 NaN,NaN,NaN,NaN,NaN,  0];
di.alpha=[0,0,0,0,0,0];
di.dalphadT=[0,0,0,0,0,0];
di.dalphadP=[0,0,0,0,0,0];
di.DQF=[0,0,0,0,-2.9,-1.5,-1.0];
di.dDQFdT=[0,0,0,0,0,0,0];
di.dDQFdP=[0,0,0,0,0,0,0];
          
o.minrl='o';
o.sys='NCFMASO';
o.endmember=[50,48,49,51,231,232,233];
o.endmass=[HP98(50).mass,HP98(48).mass,HP98(49).mass,HP98(51).mass,HP98(231).mass,HP98(232).mass,HP98(233).mass];
o.endnum=7;
o.ensym={'jd','di','hed','acm','om','cfm','jac'};
o.variables={'x(o)';'j(o)';'f(o)';'Q(o)';'Qa(o)';'Qf(o)'};
o.p={@(X)X(2)*(1-X(3))-X(4)-X(5);
              @(X)(1-X(2)-X(4))*(1-X(1)+X(6))-2*X(4)*X(1);
              @(X)(1-X(2)-X(4))*(  X(1)+X(6));
              @(X)X(3)*X(2)-X(5);
              @(X)2*X(4);
              @(X)2*X(4)*X(1)+2*X(6)*(X(2)-1+X(4));
              @(X)2*X(5)}';
o.a_ideal={ @(X)(                X(2)*(1-X(3))-X(4)+X(5))^(1/2) * (    X(2)*(1-X(3))+X(4)-X(5))^(1/2) * (  X(2)-X(4))^(1/2) * (  X(2)+X(4))^(1/2);
                    @(X)((1-X(2)+X(4))*(1-X(1)+X(6))-2*X(4)*X(6))^(1/2) * ((1-X(2)-X(4))*(1-X(1)-X(6)))^(1/2) * (1-X(2)+X(4))^(1/2) * (1-X(2)-X(4))^(1/2);
                    @(X)((1-X(2)+X(4))*(  X(1)-X(6))+2*X(4)*X(6))^(1/2) * ((1-X(2)-X(4))*(  X(1)+X(6)))^(1/2) * (1-X(2)+X(4))^(1/2) * (1-X(2)-X(4))^(1/2);
                    @(X)(                X(2)*   X(3)      -X(5))^(1/2) * (    X(2)*   X(3)      +X(5))^(1/2) * (  X(2)-X(4))^(1/2) * (  X(2)+X(4))^(1/2);
                    @(X)((1-X(2)+X(4))*(1-X(1)+X(6))-2*X(4)*X(6))^(1/2) * (    X(2)*(1-X(3))+X(4)-X(5))^(1/2) * (1-X(2)+X(4))^(1/2) * (  X(2)+X(4))^(1/2);
                    @(X)((1-X(2)+X(4))*(  X(1)-X(6))+2*X(4)*X(6))^(1/2) * ((1-X(2)-X(4))*(1-X(1)-X(6)))^(1/2) * (1-X(2)+X(4))^(1/2) * (1-X(2)-X(4))^(1/2);
                    @(X)(                X(2)*(1-X(3))-X(4)+X(5))^(1/2) * (    X(2)*   X(3)      +X(5))^(1/2) * (  X(2)-X(4))^(1/2) * (  X(2)+X(4))^(1/2)}';
o.model='sf';
o.W=[   26,   24,    5, 15.5, 25.2,    3;
                NaN,    4,   15,15.75,    2,21.05;
                NaN,  NaN,   14, 17.2,    2, 20.1;
                NaN,  NaN,  NaN, 12.8, 15.5,    3;
                NaN,  NaN,  NaN,  NaN,18.45, 19.3;
                NaN,  NaN,  NaN,  NaN,  NaN,21.05];
o.dWdT=[  0,  0,  0,  0,  0,  0;
                 NaN,  0,  0,  0,  0,  0;
                 NaN,NaN,  0,  0,  0,  0;
                 NaN,NaN,NaN,  0,  0,  0;
                 NaN,NaN,NaN,NaN,  0,  0;
                 NaN,NaN,NaN,NaN,NaN,  0];
o.dWdP=[  0,  0,  0,  0,  0,  0;
                 NaN,  0,  0,  0,  0,  0;
                 NaN,NaN,  0,  0,  0,  0;
                 NaN,NaN,NaN,  0,  0,  0;
                 NaN,NaN,NaN,NaN,  0,  0;
                 NaN,NaN,NaN,NaN,NaN,  0];
o.alpha=[0,0,0,0,0,0];
o.dalphadT=[0,0,0,0,0,0];
o.dalphadP=[0,0,0,0,0,0];
o.DQF=[0,0,0,0,-2.9,-1.5,-1.0];
o.dDQFdT=[0,0,0,0,0,0,0];
o.dDQFdP=[0,0,0,0,0,0,0];
          
% -----------------------------------------------------------------------

% 87:afchl;85:clin;88:daph;86:ames
chl.minrl='chl';
chl.sys='FMASH';
chl.endmember=[87,85,88,86];
chl.endmass=[HP98(87).mass,HP98(85).mass,HP98(88).mass,HP98(86).mass];
chl.endnum=4;
chl.ensym={'afchl','clin','daph','ames'};
chl.variables={'x(chl)','y(chl)','Q(chl)'};
chl.p={@(X)1-X(2)-X(3);
             @(X)2*X(3)-2/5*X(1)*(3-X(2));
             @(X)2/5*X(1)*(3-X(2));
             @(X)X(2)-X(3)}';
chl.a_ideal={@(X)  (1-X(1))^4 *((1-X(1))*(1-X(2)+X(3))) *((1-X(1))*(1-X(2)-X(3))) *(1-X(2))^2;
                   @(X)4*(1-X(1))^4 *((1-X(1))*(1-X(2)+X(3))) *(X(2)+X(3))              *X(2)      *(1-X(2));
                   @(X)4*X(1)^4     *(X(1)    *(1-X(2)+X(3))) *(X(2)+X(3))              *X(2)      *(1-X(2));
                   @(X)  (1-X(1))^4 *(X(2)-X(3))              *(X(2)+X(3))              *X(2)^2}';
chl.check=[1,0,1,0;0,1/2,1/2,1;0,1/2,1/2,0];
chl.model='sf';
chl.W=[18,14.5,20;
             NaN,2.5,18;
             NaN,NaN,13.5];
chl.dWdT=[0,0,0;NaN,0,0;NaN,NaN,0];
chl.dWdP=[0,0,0;NaN,0,0;NaN,NaN,0];
chl.alpha=[];
chl.dalphadT=[];
chl.dalphadP=[];
chl.DQF=[0,0,0,0];
chl.dDQFdT=[0,0,0,0];
chl.dDQFdP=[0,0,0,0];

% -----------------------------------------------------------------------

% 8:alm;7:py;10:gr;11:andr
g.minrl='g';
g.sys='CFMASO';
g.endmember=[8,7,10,11];
g.endmass=[HP98(8).mass,HP98(7).mass,HP98(10).mass,HP98(11).mass];
g.endnum=4;
g.ensym={'alm','py','gr','andr'};
g.variables={'x(g)';'z(g)';'f(g)'};
g.p={@(X)(1-X(2))*X(1);
            @(X)(1-X(2))*(1-X(1));
            @(X)X(2)-X(3);
            @(X)X(3)}';
g.a_ideal={@(X)((1-X(2))*X(1)    )^3 * (1-X(3))^2;
                  @(X)((1-X(2))*(1-X(1)))^3 * (1-X(3))^2;
                  @(X)X(2)               ^3 * (1-X(3))^2;
                  @(X)X(2)               ^3 * X(3)    ^2}';
g.model='asf';
g.W=[2.5, 10, 75
            NaN, 45, 90
            NaN,NaN,  0];
g.dWdT=[  0,  0,  0;
               NaN,  0,  0;
               NaN,NaN,  0];
g.dWdP=[  0,  0,  0;
               NaN,  0,  0;
               NaN,NaN,  0];
g.alpha=[1,1,3,3];
g.dalphadT=[0,0,0,0];
g.dalphadP=[0,0,0,0];
g.DQF=[0,0,0,-0.0044];
g.dDQFdT=[0,0,0,0.00014];
g.dDQFdP=[0,0,0,0];

% -----------------------------------------------------------------------

% 29:cz;31:ep;30:fep
ep.minrl='ep';
ep.sys='FASHO';
ep.endmember=[29,31,30];
ep.endmass=[HP98(29).mass,HP98(31).mass,HP98(30).mass];
ep.endnum=3;
ep.ensym={'cz','ep','fep'};
ep.variables={'f(ep)','Q(ep)'};
ep.p={@(X)1-X(1)-X(2);
            @(X)2*X(2);
            @(X)X(1)-X(2)}';
ep.a_ideal={@(X)(1 - X(1) + X(2)) * (1 - X(1) - X(2));
                  @(X)(1 - X(1) + X(2)) * (    X(1) + X(2));
                  @(X)(    X(1) - X(2)) * (    X(1) + X(2))}';
ep.check=[0,0;0.5,0.5;1,0];
ep.model='ideal';
ep.W=[];
ep.dWdT=[];
ep.dWdP=[];
ep.alpha=[];
ep.dalphadT=[];
ep.dalphadP=[];
ep.DQF=[0,0,0];
ep.dDQFdT=[0,0,0];
ep.dDQFdP=[0,0,0];

% -----------------------------------------------------------------------

% 101: abh; 104: an;
pl.minrl='pl';
pl.sys='NCAS';
pl.endmember=[101,104];
pl.endmass=[HP98(101).mass,HP98(104).mass];
pl.endnum=2;
pl.ensym={'abh','an'};
pl.variables={'ca(pl)'};
pl.p={@(X)1-X(1),@(X)X(1)};
pl.a_ideal={@(X)1-X(1),@(X)X(1)};
pl.model='asf';
pl.W=[3.1];
pl.dWdT=[0];
pl.dWdP=[0];
pl.alpha=[0.643,1];
pl.dalphadT=[0,0];
pl.dalphadP=[0,0];
pl.DQF=[0,7.03];
pl.dDQFdT=[0,-0.00466];
pl.dDQFdP=[0,0];

% -----------------------------------------------------------------------

% 75:mu;78:pa;76:cel;77:fcel
mu.minrl='mus';
mu.sys='NKFMASH';
mu.endmember=[75,78,76,77];
mu.endmass=[HP98(75).mass,HP98(78).mass,HP98(76).mass,HP98(77).mass];
mu.endnum=4;
mu.ensym={'mu','pa','cel','fcel'};
mu.variables={'x(mu)','y(mu)','na(mu)'};
mu.p={@(X)X(2)-X(3);
              @(X)X(3);
              @(X)(1-X(1))*(1-X(2));
              @(X)X(1)*(1-X(2))}';
mu.a_ideal={@(X)4* (1-X(3))* X(2)^2/2         *(1-X(2)/2);
                    @(X)4*    X(3) * X(2)^2/2         *(1-X(2)/2);
                    @(X)   (1-X(3))* (1-X(1))*(1-X(2))*(1-X(2)/2)^2;
                    @(X)   (1-X(3))*    X(1) *(1-X(2))*(1-X(2)/2)^2}';
mu.check=[0,1,0;0,1,1;0,0,0;1,0,0]';
mu.model='asf';
mu.W=[10.12,  0, 0;
                NaN, 52,52;
                NaN,NaN, 0];
mu.dWdT=[0.0034,0,0;
                      0,0,0;
                      0,0,0];
mu.dWdP=[0.353,0.2,0.2;
                      0,0,0;
                      0,0,0];
mu.alpha=[0.63,0.37,0.63,0.63];
mu.dalphadT=[0,0,0,0];
mu.dalphadP=[0,0,0,0];
mu.DQF=[0,0,0,0];
mu.dDQFdT=[0,0,0,0];
mu.dDQFdP=[0,0,0,0];

% -----------------------------------------------------------------------

% 80 : phl 81 : ann 201: obi 83 : east 202: tbi 203: fbi
bi.minrl='bi';
bi.sys='KFMASHTO';
bi.endmember=[80,81,201,83,202,203];
bi.endmass=[HP98(80).mass,HP98(81).mass,HP98(201).mass,HP98(83).mass,HP98(202).mass,HP98(203).mass];
bi.endnum=6;
bi.ensym={'phl','ann','obi','east','tbi','fbi'};
bi.variables={'x(bi)','y(bi)','f(bi)','t(bi)','Q(bi)'};
bi.p={@(X)1-X(1)-X(2)-X(3)-X(4)-2/3*X(5)+X(1)*(X(2)+X(3)+X(4));
               @(X)X(1) - 1/3*X(5);
               @(X)X(5) - X(1) * (X(2)+X(3)+X(4));
               @(X)X(2);
               @(X)X(4);
               @(X)X(3)}';
bi.a_ideal={@(X)4* (1-X(2)-X(3)-X(4)-2/3*X(5)-X(1)*(1-X(2)-X(3)-X(4))) * (1+1/3*X(5)-X(1))^2 *(1/2+1/2*(X(2)+X(3))) * (1/2-1/2*(X(2)+X(3))) * (1-X(4))^2;
                     @(X)4* (                 2/3*X(5)+X(1)*(1-X(2)-X(3)-X(4))) * ( -1/3*X(5)+X(1))^2 *(1/2+1/2*(X(2)+X(3))) * (1/2-1/2*(X(2)+X(3))) * (1-X(4))^2;
                     @(X)4* (                 2/3*X(5)+X(1)*(1-X(2)-X(3)-X(4))) * (1+1/3*X(5)-X(1))^2 *(1/2+1/2*(X(2)+X(3))) * (1/2-1/2*(X(2)+X(3))) * (1-X(4))^2;
                     @(X)   X(2)                                                * (1+1/3*X(5)-X(1))^2 *(1/2+1/2*(X(2)+X(3)))^2                       * (1-X(4))^2;
                     @(X)4* X(4)                                                * (1+1/3*X(5)-X(1))^2 *(1/2+1/2*(X(2)+X(3))) * (1/2-1/2*(X(2)+X(3))) * X(4)    ^2;
                     @(X)   X(3)                                                * (1+1/3*X(5)-X(1))^2 *(1/2+1/2*(X(2)+X(3)))^2                       * (1-X(4))^2;}';
bi.check=[0,0,0,0,0;
                   1,0,0,0,0;
                   1/3,0,0,0,1;
                   0,1,0,0,0;
                   0,0,0,1,0;
                   0,0,1,0,0]';
bi.model='sf';
bi.W=[  9,  3, 10,  0,  0;
               NaN,  6, -1, 10,  8;
               NaN,NaN, 10,  0,  0;
               NaN,NaN,NaN,  0,  0;
               NaN,NaN,NaN,NaN,  0];
bi.dWdT=[  0,  0,  0,  0,  0;
                  NaN,  0,  0,  0,  0;
                  NaN,NaN,  0,  0,  0;
                  NaN,NaN,NaN,  0,  0;
                  NaN,NaN,NaN,NaN,  0];
bi.dWdP=[  0,  0,  0,  0,  0;
                  NaN,  0,  0,  0,  0;
                  NaN,NaN,  0,  0,  0;
                  NaN,NaN,NaN,  0,  0;
                  NaN,NaN,NaN,NaN,  0];
bi.alpha=[];
bi.dalphadT=[];
bi.dalphadP=[];
bi.DQF=[0,-3,-10.73,0,78,0];
bi.dDQFdT=[0,0,0.0001,0,0.00013,0];
bi.dDQFdP=[0,0,0,0,0,0];
           
% -----------------------------------------------------------------------

% 93:ta;94:fta;95:tats
ta.minrl='ta';
ta.sys='FMASH';
ta.endmember=[93,94,95];
ta.endmass=[HP98(93).mass,HP98(94).mass,HP98(95).mass];
ta.endnum=3;
ta.ensym={'ta','fta','tats'};
ta.variables={'x(ta)','y(ta)'};
ta.p={@(X)1-X(2)-X(1)*(1-1/3*X(2)),@(X)X(1)*(1-1/3*X(2)),@(X)X(2)};
ta.a_ideal={@(X)(1-X(1))^3*(1-X(2))*(1-X(2)/2)^2,@(X)X(1)^3*(1-X(2))*(1-X(2)/2)^2,@(X)4*(1-X(1))^2*X(2)*(1-X(2)/2)*X(2)/2};
ta.check=[0,1,0;0,0,1];
ta.model='ideal';
ta.W=[];
ta.dWdT=[];
ta.dWdP=[];
ta.alpha=[];
ta.dalphadT=[0,0,0];
ta.dalphadP=[0,0,0];
ta.DQF=[0,0,0];
ta.dDQFdT=[0,0,0];
ta.dDQFdP=[0,0,0];

% -----------------------------------------------------------------------

% 78: pa;
pa.minrl='pa';
pa.sys='single';
pa.endmember=[78];
pa.endmass=[HP98(78).mass];
pa.endnum=1;
pa.ensym=[];
pa.variables=[];
pa.p=1;

% 32: law;
law.minrl='law';
law.sys='single';
law.endmember=[32];
law.endmass=[HP98(32).mass];
law.endnum=1;
law.ensym=[];
law.variables=[];
law.p=1;

% 100: ab;
ab.minrl='ab';
ab.sys='single';
ab.endmember=[100];
ab.endmass=[HP98(100).mass];
ab.endnum=1;
ab.ensym=[];
ab.variables=[];
ab.p=1;

% 120: ru;
ru.minrl='ru';
ru.sys='single';
ru.endmember=[120];
ru.endmass=[HP98(120).mass];
ru.endnum=1;
ru.ensym=[];
ru.variables=[];
ru.p=1;

% 43: sph;
sph.minrl='sph';
sph.sys='single';
sph.endmember=[43];
sph.endmass=[HP98(43).mass];
sph.endnum=1;
sph.ensym=[];
sph.variables=[];
sph.p=1;

% 17: ky;
ky.minrl='ky';
ky.sys='single';
ky.endmember=[17];
ky.endmass=[HP98(17).mass];
ky.endnum=1;
ky.ensym=[];
ky.variables=[];
ky.p=1;

% 28: zo;
zo.minrl='zo';
zo.sys='single';
zo.endmember=[28];
zo.endmass=[HP98(28).mass];
zo.endnum=1;
zo.ensym=[];
zo.variables=[];
zo.p=1;

% 105: q;
q.minrl='q';
q.sys='single';
q.endmember=[105];
q.endmass=[HP98(105).mass];
q.endnum=1;
q.ensym=[];
q.variables=[];
q.p=1;

