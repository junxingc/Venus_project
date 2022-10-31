%add or remove minerals if you changed ac
SP = ["Temperature","pressure","act","bi","chl","di","ep","g","gl","hb","mu","o","pl","ta","pa","law","ab","ru","sph","ky","zo","q"];
for s=1:101
SP = [SP;V1_Na_cold{s,1},V1_Na_cold{s,2},Calcmode(V1_Na_cold{s,4})];
end


function SP = Calcmode(file1)
filename = fullfile(file1);
fileID = fopen(filename);
C_data_1 = textscan(fileID,'%s');
mode = C_data_1{1};
N = length(mode)/2;
for i=1:N
    mine{1,i} = mode(i,1);
    port(1,i) = str2num(mode{i+N,1});
end
%add or remove minerals if you changed ac
SP=zeros(1,20);
for i = 1:N
     if mine{1,i} == "act"
         SP(1,1) = port(1,i);
     elseif mine{1,i} == "bi"
         SP(1,2) = port(1,i);
     elseif mine{1,i} == "chl"
         SP(1,3) = port(1,i);
     elseif mine{1,i} == "di"
         SP(1,4) = port(1,i);
     elseif mine{1,i} == "ep"
         SP(1,5) = port(1,i);
     elseif mine{1,i} == "g"
         SP(1,6) = port(1,i);
     elseif mine{1,i} == "gl"
         SP(1,7) = port(1,i);
     elseif mine{1,i} == "hb"
         SP(1,8) = port(1,i);
     elseif mine{1,i} == "mu"
         SP(1,9) = port(1,i);
     elseif mine{1,i} == "o"
         SP(1,10) = port(1,i);
     elseif mine{1,i} == "pl"
         SP(1,11) = port(1,i);
     elseif mine{1,i} == "ta"
         SP(1,12) = port(1,i);
     elseif mine{1,i} == "pa"
         SP(1,13) = port(1,i);
     elseif mine{1,i} == "law"
         SP(1,14) = port(1,i);
     elseif mine{1,i} == "ab"
         SP(1,15) = port(1,i);
     elseif mine{1,i} == "ru"
         SP(1,16) = port(1,i);
     elseif mine{1,i} == "sph"
         SP(1,17) = port(1,i);
     elseif mine{1,i} == "ky"
         SP(1,18) = port(1,i);
     elseif mine{1,i} == "zo"
         SP(1,19) = port(1,i);
     elseif mine{1,i} == "q"
         SP(1,20) = port(1,i);         
     end
end
end