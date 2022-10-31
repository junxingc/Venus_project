for j=1:7
for i=1:8
    A(1,i)=V1Nacold(54*j-54+16,i+2);
    A(2,i)=V1Nacold(54*j-54+17,i+2);
end
fileID = fopen(V1_Na_cold{j,4},'w');
for i=1:8
fprintf(fileID,'%s ',A(1,i));
end
fprintf(fileID,'\n');
for i=1:8
fprintf(fileID,'%s ',A(2,i));
end
fclose(fileID);
for i=1:31
    B(i,1)=V1Nacold(54*j-54+23+i,1);
    B(i,2)=V1Nacold(54*j-54+23+i,2);
    B(i,3)=V1Nacold(54*j-54+23+i,3);
end
fileID = fopen(V1_Na_cold{j,3},'w');
for i=1:31
fprintf(fileID,'%s %s %s\n',B(i,1),B(i,2),B(i,3));
end
fclose(fileID);
end