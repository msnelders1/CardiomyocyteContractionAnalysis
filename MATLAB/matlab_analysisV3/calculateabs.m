function absresults = calculateabs(data,time)
    figure; plot(time,data);
    
    % Detect start of peaks.
    [in,~]=find(data~=0);
    peakStart=[];
    peakStart(1)=in(1);
    for ii=2:length(in)
        if(in(ii)-in(ii-1)~=1), peakStart(end+1)=in(ii); end
    end
     % Detect end of peaks.
    in=flip(in,1);
    peakEnd=[];
    peakEnd(1)=in(1);
    for ii=2:length(in)
        if(in(ii-1)-in(ii)~=1), peakEnd=cat(2,in(ii),peakEnd); end
    end
    
    % Calculate the average maximum for each contraction.
    peakMax=0;
    for ii=1:length(peakStart)
        peakMax=peakMax+max(data(peakStart(ii):peakEnd(ii)));
    end
    peakMax=peakMax/length(peakStart);









    absresults=peakMax;
end
