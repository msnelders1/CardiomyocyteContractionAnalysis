function [allvalues,peakcontract,peakrelax,thrsh,newindval,lowmax,iv,highmin,indval2] = calculating(data,frames,maxframes,fps,dist,dt,manualthresh,time,peakPriority)
    data=data(1:frames);
    dt=dt(1:frames-1);
    
    allvalues = [0; 0; ...
                0; 0; 0; 0; ...
                0; 0; 0; ...
                0; 0; 0; ...
                0; 0; 0; ...
                0; 0; ...
                0; 0; 0];
    
    meanpresorspeed = mean(data); 
    %if there is no manual threshold, the mean will be used
    if isempty(manualthresh)
        [indval,lowmax,abrt] = findextremes(data,fps,dist,0,meanpresorspeed); 
        %indval contains time indices,y-values and frame-indices of the
        %maxima in data, lowpeaks contains the same but for the maxima below the
        %threshold
        if abrt ~= 0 
            fprintf('There are no peaks.\n')
            %figure()
            %plot(time,data)
            %title('Pressure vs time')
            %ylabel('Pa/s')
            %xlabel('Time (s)')
            return  
        end
    else
        [indval,lowmax,abrt] = findextremes(data,fps,dist,0,manualthresh);
        if abrt ~= 0 
            fprintf('There are no peaks.\n')
            %figure()
            %plot(time,data)
            %title('Pressure vs time')
            %ylabel('Pa/s')
            %xlabel('Time (s)')
            return
        end
    end
    if isempty(manualthresh)
        if ~isempty(lowmax)
            [iv,highmin] = findextremes(data,fps,1,1,min(indval(:,2))*0.85); %threshold for minima is 85% of the y-value of the lowest contraction-/relaxation-peak
            %iv contains the same information as indval (without the frame indices), but for the minima of
            %data, ht has it for the minima higher than the threshold
            thrsh1 = max(lowmax(:,2));
            thrsh2 = max(iv(:,2));
            if ~isempty(thrsh2)
                thrsh2=0;
            end
            %new threshold to determine all peak variables is either
            %the highest maximum below the mean, or the highest minimum
            %below the minimum-threshold
            if thrsh1>thrsh2
                thrsh = thrsh1;
            else
                thrsh = thrsh2;
            end
        else
            thrsh = meanpresorspeed;    
            iv = zeros(1,2);
            highmin = zeros(1,2);
        end
    else
        if ~isempty(lowmax)
            [iv,highmin] = findextremes(data,fps,1,1,manualthresh);
            thrsh = manualthresh;
        else
            thrsh = manualthresh;    
            iv = zeros(1,2);
            highmin = zeros(1,2);
        end
    end
    
    if isempty(thrsh)
        thrsh=0;
    end
       
    thrsh=ceil(thrsh*10000)/10000;
    
    [peakvalues,peaks,newindval,peakbegin,peakend] = peakarea(data,frames,fps,thrsh); 

    %cntrl=0 if the first peak is a contraction, cntrl=1 if the first peak is a relaxation
    if size(newindval,1) <= 2    %when you have two peaks you have to assume the most left peak is the contraction
        cntrl = peakPriority;
    %elseif newindval(2,3)-newindval(1,3) < newindval(3,3)-newindval(2,3)    %determine what is contraction and relaxation based on the peak distance in x
    %    cntrl = peakPriority;
    else
        cntrl = 1-peakPriority;
    end
    
%         FIXED: This old part does not remove non-peak values
%         peakcontract = peakvalues; %only keep the contraction peaks
%         for k = (2-cntrl):2:numel(peakend)
%             peakcontract(peakbegin(k):peakend(k)) = 0;
%         end
%         peakrelax = peakvalues; %only keep the relaxation peaks
%         for k = (1+cntrl):2:numel(peakend)
%             peakrelax(peakbegin(k):peakend(k)) = 0;
%         end
    
    peakbgr = peakvalues; % Only keep the background noise.
    for k = 1:1:numel(peakend)
        peakbgr(peakbegin(k):peakend(k))=0;
    end
    % Subtract the background value from the peaks.
    meanbgr = mean(nonzeros(peakbgr));
    
    peakcontract = zeros(numel(peakvalues),1); % Only keep the contraction peaks.
    for k = (1+cntrl):2:numel(peakend)
        peakcontract(peakbegin(k):peakend(k)) = peakvalues(peakbegin(k):peakend(k))-meanbgr;
    end
    peakrelax = zeros(numel(peakvalues),1); % Only keep the relaxation peaks.
    for k = (2-cntrl):2:numel(peakend)
        peakrelax(peakbegin(k):peakend(k)) = peakvalues(peakbegin(k):peakend(k))-meanbgr;
    end
    
    meanrelax = mean(nonzeros(peakrelax)); %average data of a moving pixel over all relaxation frames
    meancont = mean(nonzeros(peakcontract)); %average data of a moving pixel over all contraction frames
    
    [bpm,avgbeatduration,betweentimeaverage,betweentimemax,betweentimedev,betweencontrelaverage,avgmaxcont,avgmaxrelax,maxratio] = peakcalc(newindval,frames,fps,cntrl,peakbegin,peakend);
    [avgconttime,avgrelaxtime,timeratio,avgcontarea,avgrelaxarea,arearatio] = timecalc(peaks,cntrl);
    
    slopethresh = max(dt)/20;
    [indval2,avgmaxcontslope,avgmaxrelaxslope] = findpeakslope(dt,fps,dist,slopethresh,frames); % assumptation that contraction slope is always bigger than relaxation slope.
    %indval2 are the maxima of the derivative of b in N/m^2 per second^2 or in m/s^2, the 2
    %averages display the max slope of the contraction and relaxation in the same units as indval2
    
    tempM=zeros(maxframes,1);
    tempM(1:frames)=peakcontract;
    peakcontract=tempM;
    
    tempM=zeros(maxframes,1);
    tempM(1:frames)=peakrelax;
    peakrelax=tempM;
    
    allvalues = [bpm; avgbeatduration; ...
                avgconttime; avgrelaxtime; timeratio; betweencontrelaverage; ...
                betweentimeaverage; betweentimemax; betweentimedev; ...
                meancont; meanrelax; meanpresorspeed; ...
                avgmaxcont; avgmaxrelax; maxratio; ...
                avgmaxcontslope; avgmaxrelaxslope; ...
                avgcontarea; avgrelaxarea; arearatio];
end