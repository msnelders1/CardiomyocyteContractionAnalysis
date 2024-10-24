function [bpm,avgbeatduration,betweenbeatavg,betweenbeatdif,betweenbeatstd,betweencontrelaverage,avgmaxcont,avgmaxrelax,maxratio] = peakcalc(indval,frames,fps,cntrl,peakbegin,peakend)
    %Measure frequency
    if size(indval,1) < 3
        spb = frames/fps;
        bpm = 0; %don't display this maximum frequency in allvalues
        %fprintf('There is only 1 contraction in this timeframe. The frequency is smaller than %6.4f beats per minute.\n',60/spb)
        avgbeatduration=0;
        betweenbeatavg=0;
        betweenbeatdif=0;
        betweenbeatstd=0;
        betweencontrelaverage=0;
        avgmaxcont=0;
        avgmaxrelax=0;
        maxratio=0;
        return;
    else
        spb = (max(indval(1:2:end,1))-min(indval(1:2:end,1))) ./ ceil(size(indval,1)/2-1);
        bpm = 60/spb;
        %fprintf('The frequency is equal to %6.4f beats per minute.\n',bpm)
    end
    
    %Measure time of the beat duration
    beattime = peakend(2+cntrl:2:end)-peakbegin(1+cntrl:2:end-1); %end-1 for the case the final peak is a contraction        
    avgbeatduration = mean(beattime)/fps; %in seconds
    %fprintf('The average time of one contraction and relaxation is %6.4f seconds.\n',avgbeatduration)
    
    %Measure time between the beat's end of contraction and start of relaxation
    betweentime=peakbegin(2+cntrl:2:end)-peakend(1+cntrl:2:end-1);
    betweencontrelaverage=mean(betweentime)/fps; %seconds
    %fprintf('The average time between a contraction and relaxation is %6.4f seconds.\n',betweentimeaverage)
    
    %Measure resting time between beats
    %betweenbeataverage=spb-avgbeatduration; (see below)
    %FIXME: needs verification that the new one is better.
    betweenbeatdur=(peakbegin(3+cntrl:2:end)-peakend(2+cntrl:2:end-1))/fps;
    betweenbeatavg=mean(betweenbeatdur);
    betweenbeatdif=max(betweenbeatdur)-min(betweenbeatdur);
    if isnan(betweenbeatavg)
        betweenbeatavg=0;
        betweenbeatdif=0;
    end
        
    
    %FIXME: Just an idea, but gives a misguided value with short times.
    %Measures the consistency of beating frequency (0-1)
    %1: highly consistent, 0: inconsistent beating.
    %if(betweenbeatdif>0)
    %    betweenbeatstd=1-(std(betweenbeatdur)/betweenbeatdif);
    %else
    %    betweenbeatstd=0;
    %end
    betweenbeatstd=std(betweenbeatdur);
    
    %Measure average maximum of contraction and relaxation
    contforce = indval(1+cntrl:2:end,2); 
    relaforce = indval(2-cntrl:2:end,2);
    avgmaxcont = mean(contforce); %average maximum contraction pressure (Pa/s) or speed (m/s) 
    avgmaxrelax = mean(relaforce); %average maximum relaxation pressure (Pa/s) or speed (m/s) 
    
    %Measure ratio between contraction and relaxation pressure
    maxratio = avgmaxcont/avgmaxrelax;
    %fprintf('The average maximum peak ratio is %6.4f.\n\n',maxratio)
end
