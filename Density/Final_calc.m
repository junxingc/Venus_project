for i=1:101
D_result{i,1}=CalcD(V1_Na_cold{i,1},V1_Na_cold{i,2},V1_Na_cold{i,3},V1_Na_cold{i,4});
V1_Na_cold{i,5}=D_result{i,1}{1,1};
V1_Na_cold{i,6}=D_result{i,1}{1,2};
V1_Na_cold{i,7}=D_result{i,1}{1,3};
V1_Na_cold{i,8}=D_result{i,1}{1,4};
end