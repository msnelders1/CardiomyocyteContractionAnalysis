def findpeakslope(X,fps,minp,thrsh,frames):
    import statistics as s
    import numpy as np

    def calc(indval,frames,fps,cntrl):
        # Measure average maximum contraction and relaxation slope
        contslope = indval[cntrl:-1:2,1]
        relaslope = indval[1-cntrl:-1:2,1]
        avgmaxcontslope = s.mean(contslope)
        avgmaxrelaxslope = s.mean(relaslope)
        return avgmaxcontslope,avgmaxrelaxslope

    def allpeaks(X):
        # Looks for derivative sign change from positive to negative for maxima
        V = np.full(len(X),False)
        V[1:-2,:] = np.sign(np.diff(X[0:-2,:]))-np.sign(np.diff(X[1:-1,:])) > 1
        i = np.argwhere(V>0)
        v = X[i] # Makes array of all maximum values and their indices
        return v,i

    def funct(I,V,minp,thrsh):
        val = []
        ind = []
        if I:
            val = [val, V]
            ind = [ind, I]
            # Removes peaks that are smaller than the set threshold
            x = np.argwhere(val<thrsh)
            val[x] = []
            ind[x] = []
            # Removes peaks that are closer than a minimal value to a previous peak
            [ind,val] = removepeaks(ind,val,minp)
        return val,ind

    def removepeaks(Is,Vs,minp):
        counter=0
        rightpeak = Is[1:-1]
        leftpeak = Is[0:-2]
        difference = rightpeak-leftpeak
        tooclose = np.argwhere(difference<=minp)
        heightleft = Vs(tooclose)
        heightright = Vs(tooclose+1)
        for ii in range(0,len(heightleft)):
            if heightleft(ii)>=heightright(ii):
                Vs[tooclose(ii)+1+counter] = []
                Is[tooclose(ii)+1+counter] = []
            else:
                Vs[tooclose(ii)+counter] = []
                Is[tooclose(ii)+counter] = []
            counter=counter-1
        return Is,Vs

    [v,i] = allpeaks(X);
    [val,ind] = funct(i,v,minp,thrsh);
    indval = [ind/fps,val,ind];
    indval = np.sortrows(indval);

    if indval(1,2) > indval(2,2):
        cntrl = 0;
    else:
        cntrl = 1;

    [avgmaxcontslope,avgmaxrelaxslope] = calc(indval,frames,fps,cntrl);
    return indval,avgmaxcontslope,avgmaxrelaxslope
