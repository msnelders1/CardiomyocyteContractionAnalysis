function [peakvalues,peaks,indval,peakbegin,peakend] = peakarea(data,frames,fps,thrsh,indval)
%distinguishes peaks from noise and determines the peak positions  
    peakvalues = zeros(frames,1);
    peakindices = find(data>thrsh);
    peakvalues(peakindices) = data(peakindices); %set all values in data below the threshold to zero
    
    jj = 0;
    for ii = 1:frames                                                       %Assign each frame to its corresponding peak
        if peakvalues(ii) ~= 0   
            if ii>1 && peakvalues(ii-1)==0 || ii==1                       %New peak = new element
                jj = jj+1;
                peaks(jj,1) = 0;
                peaks(jj,2) = 0;
            end
            peaks(jj,1) = peaks(jj,1)+1;                                    %Peak time (frames)
            peaks(jj,2) = peaks(jj,2)+peakvalues(ii);                       %Peak area ((Pa/s)*frames) or ((m/s)*frames)
        end
    end

    if peakvalues(end) ~= 0 % remove peak at end of graph if it doesn't end with zero
        right = peaks(end,1);
        peaks(end,:) = []; % remove this peak from peaks
        peakvalues(end-right:end) = 0; %remove this entire peak from peakvalues
    end
    
    if peakvalues(1) ~= 0  %Remove peak at the start of the graph if it does not start with 0
        left = peaks(1,1);
        peakvalues(1:left) = 0; %remove this peak from peakvalues
        peaks(1,:) = []; %remove this peak from peaks
    end
      
    peaks(:,3) = peaks(:,1)/fps;                                            %Peak time (seconds)
    peaks(:,4) = peaks(:,2)/fps;                                            %Peak area (Pa) or (m) 
    
    peakindex = find(peakvalues ~= 0);  %determine where peaks in peakvalues start and end
    peakbegin(1) = peakindex(1); %the first is not found with the for loop
    peakend = []; %predetermine for the case there is only one peak, the for loop will then never create matrix peakend
    jj = 0;
        for n = 2:numel(peakindex)
            if peakindex(n) ~= peakindex(n-1)+1
                jj = jj+1;
                peakend(jj) = peakindex(n-1);
                peakbegin(jj+1) = peakindex(n);
            end
        end
    peakend(end+1) = peakindex(end); %the last peak is not found by the for loop
     
    %Prevent similar peaks from being called twice
    peakvalues=peakvalues+0.00001*rand(size(peakvalues));
    
    %make sure each peak has only one maximum in indval, and that all peaks in peakvalues are also in indval
    indval = zeros(numel(peakend),3);
    for m = 1:numel(peakend)
        indval(m,2) = max(peakvalues(peakbegin(m):peakend(m)));
        indval(m,3) = find(peakvalues == indval(m,2));
    end
    indval(:,1) = indval(:,3)/fps;
    indval = sortrows(indval);
end