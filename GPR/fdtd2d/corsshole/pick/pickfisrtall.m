function d_obs=pickfisrtall(gather,n1)
%
%   first arrival pick
%
[m,n,no]=size(gather);
data=[];
ne=round(n1*1.5);
MCM=zeros(1,n);
d_obs=zeros(no*n,1);
hw=waitbar(0,'处理中 ... 请稍后 ...');
for i=1:no
    data(:,:)=gather(:,:,i);
    ER=zeros(1,m-n1+1);
    for j=1:n   %按道循环
        SS=data(:,j)./max(data(:,j));
        for t=n1:m
            ER(t)=sum(SS(t-n1+1:t).^2)/(sum(SS(1:t).^2)+0.02);%beta
        end
        index=REPS(ER,ne);
        [~,bb]=max(diff(index));
        MCM(j)=bb;
    end
    d_obs((i-1)*no+1:i*no)=MCM(:);
	waitbar(i/no);
end
close(hw)
return

function out = REPS(s,ne)
% s为单道数据
% ne
% ne=27;
ns=length(s);
tianjia1=ones(1,ne-1)*s(1);
tianjia2=ones(1,ne-1)*s(ns);
data=[tianjia1,s,tianjia2];
%---------------------------------------eps----------------------------
nss=length(data);
c=zeros(ne,ne);
ss=zeros(ne,1);
out=zeros(ns,1);
for i=ne:nss-ne+1
    for j=1:ne
        c(j,:)=data(i+j-ne:i+j-1);
        ss(j)=std(c(j,:));
    end
    [~,mindex]=min(ss);
    ave=mean(c(mindex,:));
    out(i-ne+1)=ave;
end
return