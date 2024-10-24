function [draws]=lineplot(fig,fTitle,data,time,peakcontract,peakrelax,thrsh,indval,lowpeaks,iv,ht,frames,peaklimit)
    %Makes the line plot of the data with the threshold and the areas.
    time=time(1:frames);
    data=data(1:frames);
    peakcontract=peakcontract(1:frames);
    peakrelax=peakrelax(1:frames);
    
    [~,~,~,t1,t2] = peakarea(data,frames,1,thrsh);
    t1=max(t1(1)-1,1); t2=min(t2(end)+1,frames);
    data=data(t1:t2);
    time=time(1:t2-t1+1);
    peakcontract=peakcontract(t1:t2);
    peakrelax=peakrelax(t1:t2);
    
    if(isa(fig,'matlab.ui.Figure'))
        hold on;
        draws(1)=plot(time,data,'Color','r','LineWidth',3);
        draws(2:3)=area(time,[peakcontract peakrelax]);
        draws(2).FaceColor = [0 0.5 0.5]; draws(3).FaceColor = [0 0.75 0.75];
        legend('Pressure change','Area contraction','Area relaxation','Threshold','Location','northeastoutside');
        %draws(4)=plot(indval(:,1)-time(2)*t1,indval(:,2),'bo');
        %draws(5)=plot(lowpeaks(:,1)-time(2)*t1,lowpeaks(:,2),'go');
        %draws(6)=plot(iv(:,1)-time(2)*t1,iv(:,2),'ko');
%         if(~isempty(ht(:,1)))
%                 draws(7)=plot(ht(:,1)-time(2)*t1,ht(:,2),'mo');
%                 draws(8)=plot(time,thrsh*ones(size(time)),':','Color','k','LineWidth',1);
%                 legend('Pa/s','Area contraction','Area relaxation','High maxima','Low maxima', ...
%                         'Low minima','High minima','Threshold','Location','northeastoutside');
%         else
%             draws(7)=plot(time,thrsh*ones(size(time)),':','Color','k','LineWidth',1);
%             legend('Pa/s','Area contraction','Area relaxation','High maxima','Low maxima', ...
%                     'Low minima','Threshold','Location','northeastoutside');
%         end
        
        
        xlabel('time (sec)'); xlim([time(1) time(end)]);
        ylabel('\DeltaPa'); ylim([0 peaklimit]);
        title(fTitle);
        grid on;
        hold off;
    else
        hold(fig,'on');
        draws(1)=plot(fig,time,data,'Color','r','LineWidth',3);
        draws(2:3)=area(fig,time,[peakcontract peakrelax]);
        draws(2).FaceColor=[0 0.5 0.5]; draws(3).FaceColor=[0 0.75 0.75];
        draws(4)=plot(fig,time,thrsh*ones(size(time)),':','Color','k','LineWidth',1);
        
%        Useful for debugging.
%        draws(4)=plot(fig,indval(:,1)-time(2)*t1,indval(:,2),'bo');
%         if(~isempty(lowpeaks(:,1)))
%             draws(5)=plot(fig,lowpeaks(:,1)-time(2)*t1,lowpeaks(:,2),'go');
%         end
%         draws(6)=plot(fig,iv(:,1)-time(2)*t1,iv(:,2),'ko');
%         if(~isempty(ht(:,1)))
%             draws(7)=plot(fig,ht(:,1)-time(2)*t1,ht(:,2),'mo');
%             draws(8)=plot(fig,time,thrsh*ones(size(time)),':','Color','k','LineWidth',1);
%             legend(fig,'Pressure change','Area contraction','Area relaxation','High maxima','Low maxima', ...
%                         'Low minima','High minima','Threshold','Location','northeastoutside');
%         else
%             draws(7)=plot(fig,time,thrsh*ones(size(time)),':','Color','k','LineWidth',1);
%             legend(fig,'Pressure change','Area contraction','Area relaxation','High maxima','Low maxima', ...
%                         'Low minima','Threshold','Location','northeastoutside');
%         end
        
        legend(fig,'Pressure change','Area contraction','Area relaxation','Threshold','Location','northeastoutside');
        xlabel(fig,'time (sec)'); xlim(fig,[time(1) time(end)]);
        ylabel(fig,'\DeltaPa'); ylim(fig,[0 peaklimit]);
        title(fig,fTitle);
        grid(fig,'on');
        hold(fig,'off');
    end
end
