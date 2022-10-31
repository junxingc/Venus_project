%a:file amount, 
% b:mineral numbers, 
% c:imported excel file name, 
% d:line amount for one cycle, 
% e:mode first line, 
% f:comp line amount, 
% g:comp first line,
% h:data file(with file name) name,
% example: for V3 Na warm: (98,8,V3Nawarm,47,13,26,21,V3_Na_warm)
%%% Remember to give a space to the numeber in the file
%%% Remember to import the Excel file as string
function []=txt_divid_func(a,b,c,d,e,f,g,h) 
for j=1:a
for i=1:b
    A(1,i)=c(d*j-(d-e),i+2);
    A(2,i)=c(d*j-(d-e-1),i+2);
end
fileID = fopen(h{j,4},'w');
for i=1:b
fprintf(fileID,'%s ',A(1,i));
end
fprintf(fileID,'\n');
for i=1:b
fprintf(fileID,'%s ',A(2,i));
end
fclose(fileID);
for i=1:f
    B(i,1)=c(d*j-(d-g+1)+i,1);
    B(i,2)=c(d*j-(d-g+1)+i,2);
    B(i,3)=c(d*j-(d-g+1)+i,3);
end
fileID = fopen(h{j,3},'w');
for i=1:f
fprintf(fileID,'%s %s %s\n',B(i,1),B(i,2),B(i,3));
end
fclose(fileID);
end