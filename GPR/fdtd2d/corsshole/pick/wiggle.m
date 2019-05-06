d_obs=[];
for no=1:21
%     no=1;
    data(:,:)=gather(:,no,:);
    [m,n]=size(data);
    wiggledisplay(data,1:n,tout.*1e9,'vararea',100,1)
    ylabel('\itt\rm (ns)');
    xlabel('Trace No.');
    hold on
    MCM=pickfisrt(data,66).*(tout(2)-tout(1)).*1e9;
    plot(1:n,MCM,'r')
    d_obs=[d_obs;MCM];
end