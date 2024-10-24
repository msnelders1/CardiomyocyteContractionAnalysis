function [avgconttime,avgrelaxtime,timeratio,avgcontarea,avgrelaxarea,arearatio] = timecalc(peaks,cntrl) 
    %Measure average contraction and relaxation time
    conttimes = peaks(1+cntrl:2:end,3); 
    relatimes = peaks(2-cntrl:2:end,3); 

    avgconttime = mean(conttimes); %average contraction time in seconds
    avgrelaxtime = mean(relatimes); %average relaxation time in seconds
    %fprintf('The average contraction time is %5.4f seconds.\nThe average relaxation time is %5.4f seconds.\n',avgconttime,avgrelaxtime)
    
    %Measure ratio between contraction and relaxation time
    timeratio = avgconttime/avgrelaxtime;
    %fprintf('The average peak time ratio is %5.4f.\n\n',timeratio)
    
    %Measure the average contraction and relaxation area
    avgcontarea = mean(peaks(1+cntrl:2:end,4)); %average contraction area in Pa or m
    avgrelaxarea = mean(peaks(2-cntrl:2:end,4)); %average relaxation area in Pa or m
    
    %calculate the area ratio
    arearatio = avgcontarea/avgrelaxarea;
end
