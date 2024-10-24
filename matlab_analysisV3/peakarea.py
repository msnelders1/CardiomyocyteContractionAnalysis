def peakarea(data,frames,fps,thrsh,indval):
    import numpy as np
    import random as rand
    # Distinguishes peaks from noise and determines the peak positions
    peakvalues = np.zeros(frames,1)
    peakindices = np.argwhere(data>thrsh)
    peakvalues[peakindices] = data[peakindices] # Extract all values above threshold

    jj = 0
    peaks=[]
    peakbegin=[]
    peakend=[]
    for ii in range(0,frames): # Assign each frame to its corresponding peak
        if peakvalues[ii]:
            if ii==0 or not peakvalues[ii-1]: # New peak = new element
                peaks[jj,0] = 0
                peaks[jj,1] = 0
                jj+=1
            peaks[jj,0] = peaks[jj,0]+1 # Peak time (frames)
            peaks[jj,1] = peaks[jj,1]+peakvalues[ii] # Peak area ((Pa/s)*frames) or ((m/s)*frames)

    if peakvalues[-1]: # Remove peak at end of graph if it doesn't end with zero
        right = peaks[-1,0]
        peaks[-1,:] = [] # Remove this peak from peaks
        peakvalues[-1-right:-1] = 0  # Remove this entire peak from peakvalues

    if peakvalues[0]: # Remove peak at the start of the graph if it does not start with 0
        left = peaks[0,0]
        peakvalues[0:left] = 0 # Remove this peak from peakvalues
        peaks[0,:] = [] # Remove this peak from peaks

    peaks[:,2] = peaks[:,0]/fps # Peak time (seconds)
    peaks[:,3] = peaks[:,1]/fps # Peak area (Pa) or (m)

    peakindex = np.argwhere(peakvalues>0) # Determine where peaks in peakvalues start and end
    peakbegin[0] = peakindex[0] # The first is not found with the for loop
    peakend = [] # Predetermine for the case there is only one peak, the for loop will then never create matrix peakend
    jj = 0
    for n in range(1,len(peakindex)):
        if peakindex[n] != peakindex(n-1)+1:
            jj = jj+1
            peakend[jj] = peakindex[n-1]
            peakbegin[jj+1] = peakindex[n]

    peakend.append(peakindex[-1]) # The last peak is not found by the for loop

    # Prevent similar peaks from being called twice
    peakvalues=peakvalues+0.00001*rand.randint(0,len(peakvalues))

    # Make sure each peak has only one maximum in indval, and that all peaks in peakvalues are also in indval
    indval = np.zeros(len(peakend),3)
    for m in range(1,len(peakend)):
        indval[m,1] = max(peakvalues[peakbegin[m]:peakend[m]])
        indval[m,2] = np.argwhere(peakvalues == indval[m,1])
    indval[:,0] = indval[:,2]/fps
    indval = np.sortrows(indval)

    return peakvalues,peaks,indval,peakbegin,peakend
