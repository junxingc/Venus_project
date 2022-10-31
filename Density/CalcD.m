function result = CalcD(T,P,file1,file2)
%P pressure in kbar
%T temperature in C
%file1 composition parameter file
%file2 mode file

load('Venus_NCKFMASHTO.mat')
filename = fullfile(file1);
fileID = fopen(filename);
C_data_1 = textscan(fileID,'%s %s %f','CommentStyle','//');
filename = fullfile(file2);
fileID = fopen(filename);
C_data_2 = textscan(fileID,'%s');
comp = C_data_1{3};
para = C_data_1{2};
mode = C_data_2{1};
N = length(mode)/2;
Nc = length(para);
for i=1:N
    mine{1,i} = mode(i,1);
    port(1,i) = str2num(mode{i+N,1});
end
%add or remove minerals if you changed ac
for i = 1:N
     if mine{1,i} == "act"
         mine{1,i} = act;
     elseif mine{1,i} == "bi"
         mine{1,i} = bi;
     elseif mine{1,i} == "chl"
         mine{1,i} = chl;
     elseif mine{1,i} == "di"
         mine{1,i} = di;
     elseif mine{1,i} == "ep"
         mine{1,i} = ep;
     elseif mine{1,i} == "g"
         mine{1,i} = g;
     elseif mine{1,i} == "gl"
         mine{1,i} = gl;
     elseif mine{1,i} == "hb"
         mine{1,i} = hb;
     elseif mine{1,i} == "mu"
         mine{1,i} = mu;
     elseif mine{1,i} == "o"
         mine{1,i} = o;
     elseif mine{1,i} == "pl"
         mine{1,i} = pl;
     elseif mine{1,i} == "ta"
         mine{1,i} = ta;
     elseif mine{1,i} == "pa"
         mine{1,i} = pa;
     elseif mine{1,i} == "law"
         mine{1,i} = law;
     elseif mine{1,i} == "ab"
         mine{1,i} = ab;
     elseif mine{1,i} == "ru"
         mine{1,i} = ru;
     elseif mine{1,i} == "sph"
         mine{1,i} = sph;
     elseif mine{1,i} == "ky"
         mine{1,i} = ky;
     elseif mine{1,i} == "zo"
         mine{1,i} = zo;
     elseif mine{1,i} == "q"
         mine{1,i} = q;         
     end
end
for i=1:N
    endnum(i)=mine{i}.endnum;
end
endmember=zeros(Nc,1);
k=1;
for i=1:N
    for j=1:endnum(i)
        endmember(k)=mine{i}.endmember(j);k=k+1;      
    end
end
k=0;
for i=1:N
   if mine{1,i}.sys == "single"
       M = mine{1,i}.endmass; 
       V = CalcV(P,T,HP98(mine{1,i}.endmember))
       D(1,i) = M/V; 
   else
       n = mine{1,i}.endnum-1;
       X= zeros(1,n);
       for j= 1:n
           X(1,j) = comp(k+j,1); 
       end
       V = CalcV_ss(P,T,mine{1,i},X)        
       p=ones(1,n+1);
       m=zeros(1,n+1)
       for j=1:n+1
            f=mine{1,i}.p{j};
            p(j)=f(X);
            m(j)=p(j)*mine{1,i}.endmass(j);
       end
       M = sum(m);
       D(1,i) = M/V;
       k = k + n;
   end
end
    %   single mineral density
        %   num = mineral(j).paranumber;
        %   para(1 to num) = comp(t+1 to t+num);
        %   t = t+num;
        %   density(1,i) = mineral(j).f(para(1 to num)))*p(1,i);
    %single density *mode & preserve
Density=0;
for i=1:N
    Density = Density + D(1,i)*port(1,i);
end
result = {Density,mine,port,D}
end