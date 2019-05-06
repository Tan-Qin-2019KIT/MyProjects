close all; clear all;

filename='test_TQ';
zoom=1;

seis_file=char(['../par/su/' filename '_curl.su'],...
    ['../par/su/' filename '_div.su'],...
    ['../par/su/' filename '_p.su'],...
    ['../par/su/' filename '_vx.su'],...
    ['../par/su/' filename '_vy.su']);

for i=1:length(seis_file(:,1))
    seis_type=seis_file(i,:);
    seis_type(find(isspace(seis_type)))=[];
    [data]=su2matlab(seis_type);
    delta=data(1).dt*1e-6;
    nt=data(1).ns;
    timescale=delta:delta:delta*nt;
    amp_scale_fact=0.5;             %only for display, factor for scaling the amplitudes
    
    %tracewise normalizing
    for ii=1:1:length(data)
        source_order(ii)=data(ii).tracl;
        offsets(ii)=data(ii).offset*1e-3;
        data_trace(ii,1:nt)=data(ii).trace;
        data_trace(ii,1:nt)=data_trace(ii,:)./max(data_trace(ii,:)).*amp_scale_fact;
    end
    
    figure(i);
    hold on
    for ii=1:zoom:length(offsets)
        plot(offsets(ii)+data_trace(ii,:)*(offsets(2)-offsets(1))*zoom,timescale,'LineWidth',2);
    end
    hold off
    
    if i==1
        title('seismic gather (curl)');
    elseif i==2
        title('seismic gather (div)');
    elseif i==3
        title('seismic gather (p)');
    elseif i==4
        title('seismic gather (vx)');
    else
        title('seismic gather (vy)');
    end
    
    ylabel('Time in s');
    % xlabel('Trace number');
    xlabel('Distance in m');
    % legend('vert vy','hor vz','Diff');
    set(get(gca,'title'),'FontSize',16);
    set(get(gca,'title'),'FontWeight','bold');
    set(get(gca,'Ylabel'),'FontSize',16);
    set(get(gca,'Ylabel'),'FontWeight','bold');
    set(get(gca,'Xlabel'),'FontSize',16);
    set(get(gca,'Xlabel'),'FontWeight','bold');
    set(gca, 'YDir', 'reverse')
    set(gca,'FontSize',14);
    set(gca,'FontWeight','bold');
    set(gca,'Linewidth',1.0);
    set(gca,'Box','on');
    axis([offsets(1)-zoom*(offsets(2)-offsets(1)) offsets(ii)+zoom*(offsets(2)-offsets(1)) 0 delta*nt]);
    pause(0.1)
end