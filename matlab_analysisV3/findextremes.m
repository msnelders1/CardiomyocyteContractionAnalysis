function [indval,removedextr,abrt] = findextremes(data,fps,minp,extr,thrsh)
%Finds the appropriate maxima/minima of the data
    if extr == 0 %find maxima
        [v,i] = allmax(data);
        [val,ind,vallowmax,indlowmax,abrt] = funcmax(i,v,minp,thrsh);  
        indval = [ind/fps,val,ind];     
        indval = sortrows(indval);
        removedextr = [indlowmax/fps,vallowmax,indlowmax];
    elseif extr == 1 %find minima
        [v,i] = allmin(-data);
        [indval,removedextr] = funcmin(i,v,fps,thrsh);
        if isempty(indval)
            [indval,removedextr] = funcmin(i,v,fps,mean(data)); %if there're no minima on 0.85*mean, search again using the mean
        end
    end
           
end

%%
function [v,i] = allmax(data)
    V = false(size(data));
    V(2:end-1,:) = sign(diff(data(1:end-1,:)))-sign(diff(data(2:end,:))) > 1;     %Looks for derivative sign change from positive to negative for maxima
    i = find(V>0);
    v = data(i);                                                               %Makes array of all maxima values and their indices
end

%%
 function [v,i] = allmin(data)
    V = false(size(data));
    V(2:end-1,:) = sign(diff(data(1:end-1,:)))-sign(diff(data(2:end,:))) > 1;     %Looks for derivative sign change from negative to positive for minima
    i = find(V>0);
    v = -data(i);                                                              %Makes array of all minima values and their indices
 end
%%
function [val,ind,vallowmax,indlowmax,abrt] = funcmax(I,V,minp,thrsh)
    val = [];
    ind = [];
    vallowmax = [];
    indlowmax = [];
    abrt=0;
    if ~isempty(I)            
        val = [val; V];
        ind = [ind; I];          
        i = find(val<thrsh);                                                %Removes peaks that are smaller than the set threshold
        val(i) = [];
        ind(i) = [];
            
        vallowmax = V;
        indlowmax = I;           
        k = find(vallowmax>thrsh);
        vallowmax(k) = []; 
        indlowmax(k) = [];
        
        [ind,val] = removeextremes(ind,val,minp);                             %Removes peaks that are closer than a minimal value to a previous peak 
    end
    if isempty(val)
        abrt = 1;
        ind = 0;
        val = 0;
    end
end

%%
function [indval,highmin] = funcmin(I,V,fps,thrsh)
    val = [];
    ind = [];
    if ~isempty(I)                                         
        val = [val; V];
        ind = [ind; I];
        indval = [ind,val];     
        indval = sortrows(indval);
        i = indval(:,2)>thrsh;                                              %Removes throughs that are larger than the set threshold
        indval(i,:) = [];
        ind = indval(:,1);
        val = indval(:,2); 
        valhighmin = V;
        indhighmin = I;
        highmin = [indhighmin/fps,valhighmin];
        highmin = sortrows(highmin);
        k = highmin(:,2)<thrsh;
        highmin(k,:) = [];        
        indval = [ind/fps,val];
    end
end

%% if there are 2 maxima close together, only save the highest
function [Is,Vs] = removeextremes(Is,Vs,minp)
    counter=0;
    rightpeak = Is(2:end);
    leftpeak = Is(1:end-1);
    difference = rightpeak-leftpeak;    
    tooclose = find(difference<=minp);
    heightleft = Vs(tooclose);
    heightright = Vs(tooclose+1);
    for i=1:numel(heightleft)
        if heightleft(i)>=heightright(i) 
            Vs(tooclose(i)+1+counter) = [];
            Is(tooclose(i)+1+counter) = [];
        else
            Vs(tooclose(i)+counter) = [];
            Is(tooclose(i)+counter) = []; 
        end
        counter=counter-1;
    end
end



