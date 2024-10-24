function [indval,avgmaxcontslope,avgmaxrelaxslope] = findpeakslope(X,fps,minp,thrsh,frames)
    [v,i] = allpeaks(X); 
    
    [val,ind] = funct(i,v,minp,thrsh);
    
    indval = [ind/fps,val,ind];     
    indval = sortrows(indval);
    
    if indval(1,2) > indval(2,2)
        cntrl = 0;
    else
        cntrl = 1;
    end
    
    [avgmaxcontslope,avgmaxrelaxslope] = calc(indval,frames,fps,cntrl);

    %%
    function [v,i] = allpeaks(X)
        V = false(size(X));
        V(2:end-1,:) = sign(diff(X(1:end-1,:)))-sign(diff(X(2:end,:))) > 1;     %Looks for derivative sign change from positive to negative for maxima
        i = find(V>0);
        v = X(i);                                                               %Makes array of all maximum values and their indices
    end

    %%
    function [val,ind] = funct(I,V,minp,thrsh)
        val = [];
        ind = [];
        if ~isempty(I)                                                    
            val = [val; V];
            ind = [ind; I];          
            x = find(val<thrsh);                                                %Removes peaks that are smaller than the set threshold
            val(x) = [];
            ind(x) = [];

            [ind,val] = removepeaks(ind,val,minp);                              %Removes peaks that are closer than a minimal value to a previous peak  
        end
    end

    %%
    function [Is,Vs] = removepeaks(Is,Vs,minp)
        counter=0;
        rightpeak = Is(2:end);
        leftpeak = Is(1:end-1);
        difference = rightpeak-leftpeak;    
        tooclose = find(difference<=minp);
        heightleft = Vs(tooclose);
        heightright = Vs(tooclose+1);
        for ii=1:numel(heightleft)
            if heightleft(ii)>=heightright(ii) 
                Vs(tooclose(ii)+1+counter) = [];
                Is(tooclose(ii)+1+counter) = [];
            else
                Vs(tooclose(ii)+counter) = [];
                Is(tooclose(ii)+counter) = []; 
            end
            counter=counter-1;
        end
    end

    %%
    function [avgmaxcontslope,avgmaxrelaxslope] = calc(indval,frames,fps,cntrl)
        %Measure average maximum contraction and relaxation slope
        contslope = indval(1+cntrl:2:end,2);
        relaslope = indval(2-cntrl:2:end,2);
        avgmaxcontslope = mean(contslope);
        avgmaxrelaxslope = mean(relaslope);
    end
end