def peakcalc(indval,frames,fps,cntrl,peakbegin,peakend):
    import statistics as s
    import numpy as np
    # Measure frequency
    if len(indval,1) < 3:
        spb = frames/fps
        bpm = 0 # Don't display this maximum frequency in allvalues
        avgbeatduration=0
        betweenbeatavg=0
        betweenbeatdif=0
        betweenbeatstd=0
        betweencontrelaverage=0
        avgmaxcont=0
        avgmaxrelax=0
        maxratio=0
        return
    else:
        spb = np.divide(max(indval[0:-1:2,0])-min(indval[0:-1:2,0])),np.ceil(len(indval,0)/2-1);
        bpm = 60/spb;

    # Measure time of the beat duration
    beattime = peakend[1+cntrl:-1:2]-peakbegin[cntrl:-2:2] # -2 in case the final peak is a contraction
    avgbeatduration = s.mean(beattime)/fps; # In seconds

    # Measure time between the beat's end of contraction and start of relaxation
    betweentime=peakbegin[1+cntrl:-1:2]-peakend[cntrl:-2:2];
    betweencontrelaverage=s.mean(betweentime)/fps; # In seconds

    # Measure resting time between beats
    #betweenbeataverage=spb-avgbeatduration; (see below)
    #FIXME: needs verification that the new one is better.
    betweenbeatdur=(peakbegin[2+cntrl:-1:2]-peakend[1+cntrl:-2:2])/fps;
    betweenbeatavg=s.mean(betweenbeatdur);
    betweenbeatdif=max(betweenbeatdur)-min(betweenbeatdur);
    if np.nan(betweenbeatavg):
        betweenbeatavg=0
        betweenbeatdif=0

    #Just an idea, but gives a misguided value with short times.
    #Measures the consistency of beating frequency (0-1)
    #1: highly consistent, 0: inconsistent beating.
    #if(betweenbeatdif>0)
    #    betweenbeatstd=1-(std(betweenbeatdur)/betweenbeatdif);
    #else
    #    betweenbeatstd=0;
    #end
    betweenbeatstd=s.stdev(betweenbeatdur)

    # Measure average maximum of contraction and relaxation
    contforce = indval[cntrl:-1:2,1];
    relaforce = indval[1-cntrl:-1:2,1];
    avgmaxcont = s.mean(contforce); # Average maximum contraction pressure (Pa/s) or speed (m/s)
    avgmaxrelax = s.mean(relaforce); # Average maximum relaxation pressure (Pa/s) or speed (m/s)

    # Measure ratio between contraction and relaxation pressure
    maxratio = avgmaxcont/avgmaxrelax;
    return bpm,avgbeatduration,betweenbeatavg,betweenbeatdif,\
        betweenbeatstd,betweencontrelaverage,avgmaxcont,avgmaxrelax,maxratio
