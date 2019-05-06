clear all;close all;
filepath='D:\MATLAB\GPR\test.txt';
fid = fopen(filepath,'r');
if fid == -1
    disp('can not open file!');
end

nextline = fgetl(fid);
DATA.FreqFrom = fscanf(fid,'%f',1);
DATA.FreqTo  = fscanf(fid,'%f',1);
DATA.FreqNum = fscanf(fid,'%d\n',1);
DATA.FreqList= fscanf(fid,'%f\n',DATA.FreqNum );
