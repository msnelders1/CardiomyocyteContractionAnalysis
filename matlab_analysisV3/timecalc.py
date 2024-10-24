def timecalc(peaks,cntrl):
    import statistics as s
    # Measure average contraction and relaxation time
    conttimes = peaks[cntrl:-1:2,2];
    relatimes = peaks[1-cntrl:-1:2,2];

    avgconttime = s.mean(conttimes); # Average contraction time in seconds
    avgrelaxtime = s.mean(relatimes); # Average relaxation time in seconds

    # Measure ratio between contraction and relaxation time
    timeratio = avgconttime/avgrelaxtime;

    # Measure the average contraction and relaxation area
    avgcontarea = s.mean(peaks[cntrl:-1:2,3]); # Average contraction area in Pa or m
    avgrelaxarea = s.mean(peaks[1-cntrl:-1:2,3]); # Average relaxation area in Pa or m

    # Calculate the area ratio
    arearatio = avgcontarea/avgrelaxarea;
    return avgconttime,avgrelaxtime,timeratio,avgcontarea,avgrelaxarea,arearatio
