def calculateabs(data,time):
    import numpy as np

    # Detect start of peaks.
    inn=np.argwhere(data)
    peakStart=[];
    peakStart[0]=inn[0];
    for ii in range(1,len(inn)):
        if(inn[ii]-inn[ii-1]!=1):
            peakStart.append(inn[ii])

    # Detect end of peaks.
    inn=np.flip(inn,0);
    peakEnd=[];
    peakEnd[0]=inn[0];
    for ii in range(1,len(inn)):
        if(inn[ii-1]-inn[ii]!=1):
            peakEnd.append(inn[ii])

    # Calculate the average maximum for each contraction.
    peakMax=0
    for ii in range(0,len(peakStart)):
        peakMax=peakMax+max(data[peakStart[ii]:peakEnd[ii]])
    peakMax=peakMax/len(peakStart);

    return peakMax;
