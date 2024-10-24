%% Cardiomyocyte Contraction Analysis
%{
Version 3.4
TODO: Center all gui windows on the screen.
TODO: Exclude videos without enough peaks.
TODO: Higher framerate is higher Pa? Is this correct?
TODO: Combine the waitbars into a single window.
TODO: Cumulative plots; normalize and start at the first contraction peak.
TODO: Stiffness can be set per graph, like threshold.
TODO: Load previously used threshold and peak priority by opening
results.csv
TODO: List the number of peaks/contractions (need at least 5 for a valid
TODO: Load thresholds from previous run if the results.csv file is in the
dir.
deltaPa at different recording speeds are not comparable right now.
FIXME1: Threshold and switch peaks buttons are put in queue (execution should interrupt).
FIXME1: Changing FPS is unreliable. Find out how to best solve this.

17jul2023:
    - Fixed unable to load movies in subfolders.
15jun2023:
    - Added saving threshold in DataVal file.
    - Fixed issues the normalized contraction pressure graph.
22nov2022:
    - Added a checkbox to exclude data from analysis.
21oct2022:
    - Analysis now stores peak priorities,exclusion and stiffness for subsequent runs.
18oct2022:
    - Fixed issue where thresholds were rounded before final calculations.
    - Added correction for relaxation peaks exceeding contraction peaks.
01sep2022:
    - Removed unnecessary graph components for faster graph drawing.
18aug2022:
    - Updated to MATLAB version R2022a.
04jan2022:
    - Added normalized accumulative peak graph exporting.
    - Added manually setting Shear Modulus per movie.
27dec2021:
    - Fixed peakcalc.m -> error caused by return; out-variables were not set.
22nov2021:
    - Fixed graphmaking waitbar into a single window.
28feb2021:
    - Fixed peakcalc.m -> beating variation variables showing NaN.
    - Added Rcsv2mat.m -> background subtraction.
19feb2021:
    - Fixed lineplot.m -> crashes caused by empty values.
27aug2020:
    - Fixed duplicate row names in the final result table.
30jul2020:
    - Fixed X-axis (time) limits.
    - Fixed Y-axis (Pa/s) limits to scale according to highest peak in all movies.
29jul2020:
    - Fixed uncalculated peaks showing up on graphs (ends are cropped to analysed region).
28jul2020:
    - Fixed time starting at t=1 instead of t=0.
02jul2020:
	- Added loading bars.
	- Fixed missing high minima values causing crashes.
29jun2020:
	- Redesigned user input windows.
27jun2020:
    - Added results are now exported to selected folder as results.csv.
    - Added results table is now shown after running the analysis.
    - Added graph exporting is now optional.
26jun2020:
    - Added ImageProcessing.ijm -> pixel resolution added to .csv file.
    - Added RCsv2mat.m -> gets pixel resolution from .csv file.
    - Fixed RCsv2mat.m -> No longer relies on image to get resolution.
09jun2020:
    - Added calculating.m -> peakbgr; background signal is now substracted from the peak values.
    - Added differential and accumulative force graphs.
    - Fixed calculating.m -> peakrelax/peakcontract contained non-peak values.
    - Fixed peakPriority not applying on final calculation.
    - Fixed manualThresh not applying on final calculation.
01jun2020:
    - Fixed peakcalc.m -> betweenbeatdur returning 1x0 empty vector.
05may2020:
    - Fixed usage of movies with different framecounts.
14apr2020:
    - Fixed trailing/leading spaces in r_file char vector.
06apr2020:
    - Initial version.
%}


%% Initialization
clearvars; close all force; format compact;
warning('off','MATLAB:legend:IgnoringExtraEntries');
[filepath,~] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath([filepath '\matlab_analysisV3\']);
%cd([filepath '\data\']);
%cd('V:\ddata\GENE\essers\Matthijs Snelders\Data\Data_Leiden\MuscleMotion');
%cd('\\store\department\bmw\GENE\Jeroen Essers\Heartchip\HeartCHIP_I\Example - Isoproterenol');
%cd('\\store\department\bmw\GENE\Jeroen Essers\Heartchip\MachineLearning');
%cd('V:\ddata\GENE\essers\Matthijs Snelders\Data\Data_Leiden\Raw\Processed');
%cd('\\store\department\bmw\GENE\Jeroen Essers\Heartchip\Data_Leiden\FrameComparison\Samples');
%cd('\\store\department\bmw\GENE\Jeroen Essers\Heartchip\Synthetic data contraction\synth_recordings');
cd('\\store\department\bmw\GENE\Jeroen Essers\Heartchip\Microscope Sessions\Contraction recordings\12jul2023-comptox-50ms');
clear filepath;


%% User input
%Let user input the settings
prompt={'Pixel size (um)','Shear modulus (kPa)','Frames per second'};
defaultans={'0.33','15','20'};
answer = str2double(inputdlg(prompt,'Input',1,defaultans));
if isempty(answer), fprintf('Analysis aborted. No parameters filled in.\n'); return;
end
pix_size = answer(1)*(1E-6); %size of one pixel in (m).
ShearM = answer(2)*(1E3); %shear modulus.
fps = answer(3); %frames per second.
frames=9999; %maximum framecount.
clear prompt defaultans answer;


%% Movie selection
%selecting files
folder=uigetdir(pwd,'Select a folder');
if(folder==0), fprintf('Analysis aborted. No folder selected.\n'); return;
else, files=parseFiles(folder,{});
end

%get the list of movies.
numMov=length(files); %number of movies.
fName={numMov}; %file names.
fTitle=strings(numMov,1); %file titles.
fPath=strings(numMov,1); %file path+name without extension.
numCond=0; %number of conditions.
if numMov==0, fprintf('Analysis aborted. No movies found.\n'); return;
end
for ii = 1:numMov
  %parsing movie name
  fName{ii}=split(extractBefore(files(ii).name,strlength(files(ii).name)-13),'_');
  fTitle(ii)=strjoin(string(fName{ii}));
  fPath(ii)=[files(ii).folder,filesep,files(ii).name];
  fPath(ii)=extractBefore(fPath(ii),strlength(fPath(ii))-13);
  numCond=max(numCond,length(fName{ii}));
  if(ii==1), r_file=[files(ii).folder,filesep,files(ii).name];
  else, r_file=char(r_file,[files(ii).folder,filesep,files(ii).name]);
  end
end


%% Initialize variables and constants
%preallocating variables
sumforce = zeros(frames,numMov);
sumpressure = zeros(frames,numMov);
dt = zeros(frames-1,numMov);
data = zeros(frames,numMov);
allvalues = zeros(numMov,20);
nz2=zeros(1,numMov);
framez=zeros(1,numMov);
stiffness(1:numMov)=ShearM;

t = 1:frames; %frame indices
dist = 2; %measure for maximal x-distance between peaks
time = (t-1)/fps; %time indices (s)
A = pix_size^2; %Area of one pixel (m^2)

%Ibidi micro-slide angiogenesis inner well
v_m = 10E-6*1E-3; % Volume of the well of the experiment, default=10ul (m^3)
rhe = 2E-3; % radius of the wells of the experiment, default=2mm (m)
l = v_m/(rhe^2*pi)/pix_size; %length in number of pixels (pixels)


%% Data processing
%fprintf('Processing movies.\n');
lb=waitbar(0,'Processing movies','CloseRequestFcn',[]);
movegui(lb,'center');
tic %start timing
for ii = 1:numMov
    waitbar(ii/(2*numMov),lb,['Reading out movie ' num2str(ii) ' of ' num2str(numMov)]);
    [r(ii),dX(:,ii),framez(ii),pixs] = RCsv2mat(strtrim(r_file(ii,:)),frames); %dX is the sum of the displacement over all pixels/frame
    dX(:,ii) = dX(:,ii)*fps; %sum of the movement in pixels/sec
    s = (r(ii)/framez(ii)); %average number of moving pixels/sec
    dX(:,ii) = dX(:,ii)/s; %average movement of a single moving pixel
    sumpressure(:,ii) = (stiffness(ii).*dX(:,ii))/l; % average N/m^2 of a moving pixel/sec
    data(:,ii) = sumpressure(:,ii); % average N/m^2 or Pa's of pressure/sec
    dt(:,ii) = diff(data(:,ii)); %derivative of (average N/m^2 of a moving pixel/sec) over time
end

manualthresh=[]; %manual threshold has to be empty for the if-statements in calculating
thrsh=zeros(1,numMov);
indval=cell(1,numMov);
lowmax=cell(1,numMov);
iv=cell(1,numMov);
highmin=cell(1,numMov);
peakcontract=zeros(frames,numMov);
peakrelax=zeros(frames,numMov);
peakPriority=zeros(1,numMov);
bInvalid=zeros(1,numMov);
if(isfile([folder,filesep,'DataVal.csv']))
    temp=readmatrix([folder,filesep,'DataVal.csv']);
    if(size(temp,1)==numMov)
        peakPriority=temp(1:numMov,1)';
        bInvalid=temp(1:numMov,2)';
        stiffness=temp(1:numMov,3)';
        manualthresh=temp(1:numMov,4)';
    end
end

bAbort=false;
for ii = 1:numMov
    waitbar((numMov+ii)/(2*numMov),lb,['Processing movie ' num2str(ii) ' of ' num2str(numMov)]);
    try
        if length(manualthresh)>=ii
            [allvalues(ii,:),peakcontract(:,ii),peakrelax(:,ii),thrsh(ii),indval{ii},lowmax{ii},iv{ii},highmin{ii}] = calculating(data(:,ii),framez(ii),frames,fps,dist,dt(:,ii),manualthresh(ii),time(1:framez(ii)),peakPriority(ii));
        else
            [allvalues(ii,:),peakcontract(:,ii),peakrelax(:,ii),thrsh(ii),indval{ii},lowmax{ii},iv{ii},highmin{ii}] = calculating(data(:,ii),framez(ii),frames,fps,dist,dt(:,ii),manualthresh,time(1:framez(ii)),peakPriority(ii));
        end
    catch
            fprintf('Movie "%s" produced invalid data.\n',fPath(ii)); bAbort=true;
    end
end
manualthresh=thrsh;
peaklimit=max(data);
close(lb,'force');

if bAbort, fprintf('Analysis aborted. invalid data.\n'); return;
end


%% Plotting graphs
%fprintf('Preparing graphs.\n');
[result,graphs]=graphmaking(fTitle,time,data,peakcontract,peakrelax,thrsh,indval,lowmax,iv,highmin,fps,dist,framez,frames,numMov,dt,manualthresh,allvalues,peaklimit,dX,stiffness,A,l,peakPriority,bInvalid);
endtime = toc; %end timing
waitfor(result);
doExport=questdlg('Export additional graphs?','Export graphs','Yes','No','Yes');


%% Finalizing
%fprintf('Finalizing movies.\n');
lb=waitbar(0,'Finalizing movies','CloseRequestFcn',[]);
movegui(lb,'center');
tic %start timing
for ii=1:numMov
    % Final calculations.
    %fprintf('Finalizing movie %1.0f of %1.0f.\n',ii,numMov);
    waitbar(ii/numMov,lb,['Finalizing movie ' num2str(ii) ' of ' num2str(numMov)]);
    
    sumforce(:,ii) = (stiffness(ii)*A/l).*dX(:,ii); % average N per second on a single moving pixel
    sumpressure(:,ii) = sumforce(:,ii)/A; % average N/m^2 of a moving pixel per second  
    data(:,ii) = sumpressure(:,ii); % average N/m^2 of a moving pixel per second 
    dt(:,ii) = diff(data(:,ii)); %derivative of (average N/m^2 of a moving pixel per second) over time
    [allvalues(ii,:),peakcontract(:,ii),peakrelax(:,ii),thrsh(ii),indval{ii},lowmax{ii},iv{ii},highmin{ii}] = calculating(data(:,ii),framez(ii),frames,fps,dist,dt(:,ii),manualthresh(ii),time(1:framez(ii)),peakPriority(ii));

    % Make a contraction/relaxation peak pattern.
    peaks=peakcontract(1:framez(ii),ii)-peakrelax(1:framez(ii),ii);

    % Crop to the first contraction and last relaxation.
    t1=find(peakcontract(1:framez(ii),ii)~=0);
    t2=find(peakrelax(1:framez(ii),ii)~=0);
    peaks=[0; peaks(t1(1):t2(end)); 0];
    
    % Normalize the relaxation peaks, cannot exceed contraction!
    %peaks(peaks<0)=peaks(peaks<0)*abs(max(peaks(peaks>0))/min(peaks(peaks<0)));
    arC=trapz(peaks(peaks>0));
    arR=trapz(peaks(peaks<0));
    peaks(peaks<0)=peaks(peaks<0)*abs(arC/arR);
    
    % Make the peak pattern over time.
    peak=zeros(length(peaks),1);
    for nn=1:length(peaks), peak(nn+1)=max(peak(nn)+peaks(nn),0); end
    peak=peak(2:end-1)/fps; % Divide by fps to get back to delta per timepoint.
    
    % Baseline correction.
    peaknorm=peak;
    id=find(peaks==0);
    for nn=2:length(id)-1
        if(peaks(id(nn)-1)<0)
            peaknorm(id(nn):end)=max(peaknorm(id(nn):end)-peaknorm(id(nn)),0);
        end
    end
    peaknorm(end+1)=0;
    
    if(string(doExport)=="No"), continue; end
    % Export the result graph.
    fig=figure('visible','off');
    lineplot(fig,fTitle(ii),data(:,ii),time(1:framez(ii)),peakcontract(:,ii),peakrelax(:,ii),thrsh(ii),indval{ii},lowmax{ii},iv{ii},highmin{ii},framez(ii),peaklimit(ii));
    saveas(fig,strjoin([fPath(ii),'_graph.tif'],''));
    if(bInvalid(ii)==1), continue; end
    % Export the differential pressure graph.
    fig=figure('visible','off'); grid on;
    plot(time(1:length(peaks)),peaks);
    xlabel('time (sec)'); ylabel('\DeltaPa'); grid on;
    title(strcat(fTitle(ii),' - Differential'));
    saveas(fig,strjoin([fPath(ii),'_graphDP.tif'],''));
    % Export the normalized accumulative pressure graph.
    fig=figure('visible','off');
    grid on;
    plot(time(1:length(peaknorm)),peaknorm);
    xlabel('time (sec)'); ylabel('Pa'); grid on;
    title(strcat(fTitle(ii),' - Normalized Cumulative'));
    saveas(fig,strjoin([fPath(ii),'_graphNCP.tif'],''));
end
endtime=endtime+toc; %end timing
fprintf('Elapsed time per movie is %4.1f seconds.\n',endtime/numMov);

%Write up the Thresholds we used and put them in allvalues.
allvaluestext=cat(2,stiffness'/(1E3),allvalues);
%Write up the Thresholds we used and put them in allvalues.
for ii=1:numMov
   if(manualthresh(ii)==0)
       manualthresh(ii)=thrsh(ii);
   end
end
allvaluestext=cat(2,manualthresh',allvaluestext);

%Add the Path,Movie names and Exclusions to allvalues.
fPath=fPath + "_processed.csv";
allvaluestext=cat(2,bInvalid',allvaluestext);
allvaluestext=cat(2,fTitle,allvaluestext);
allvaluestext=cat(2,fPath,allvaluestext);

%Add table headers and copy to clipboard.
tHeader=[...
    "File",...
    "Movie",...
    "Excluded",...
    "Threshold", ...
    "ShearMod (kPa)", ...
    "Beats per minute (bpm)", ...
    "Beat duration (s)", ...
    "Contraction time (s)", ...
    "Relaxation time (s)", ...
    "Time ratio", ...
    "Time between contraction and relaxation (s)", ...
    "Time between beats (s)", ...
    "Max between beats (s)", ...
    "StDev between beats (s)", ...
    "Mean contraction pressure (\DeltaPa)", ...
    "Mean relaxation pressure (\DeltaPa)", ...
    "Mean pressure per frame (\DeltaPa)", ...
    "Max contraction pressure (\DeltaPa)", ...
    "Max relaxation pressure (\DeltaPa)", ...
    "Max pressure ratio", ...
    "Max contraction slope (dP/dt)", ...
    "Max relaxation slope (dP/dt)", ...
    "Average total contraction pressure (Pa)", ...
    "Average total relaxation pressure (Pa)", ...
    "Area ratio"...
    ];
allvaluestext=cat(1,tHeader,allvaluestext);

%export the results to a .csv file.
waitbar(1,lb,'Exporting data');
writematrix(allvaluestext,[folder,filesep,'Results.csv']);
writematrix(cat(2,peakPriority',bInvalid',stiffness',manualthresh'),[folder,filesep,'DataVal.csv']);
fprintf('Results exported to %s.\n',[folder,filesep,'Results.csv']);


%% Results
%show the results in a table.
waitbar(1,lb,'Summarizing data');
res=showResults(allvaluestext,fName,bInvalid);
drawnow; close(lb,'force');
fprintf('Analysis finished.\n');
return;

%% Absolute pressure calculation
% for ii=1:numMov
%     absResults(:,ii)=calculateabs(peaksnormalized(1:framez(ii),ii),time(1:framez(ii)));
% end
% export the absolute pressure results to a .csv file.
% writematrix(absresults,[folder,filesep,'AbsResults.csv']);
% fprintf('Results exported to %s.\n',[folder,filesep,'AbsResults.csv']);














































%% Functions
function f=showResults(d,fName,bInvalid)
    %get the units per variable.
    tUnits=["Grp","bpm","s","s","s","AU","s","s","s","s","Pa/s","Pa/s","Pa/s","Pa/s","Pa/s","Ratio","dP/dt","dP/dt","Pa","Pa","AU"];
    
    %ignore excluded data and remove unneeded columns.
    fName=fName(bInvalid'==0);
    bInvalid=[0,bInvalid];
    d=d(bInvalid'==0,:);
    d(:,3:5)=[];
    d(:,1)=[];
    
    %find and assign groups.
    groups=strings(length(d(:,1))-1,1);
    for ii=1:length(fName)
        groups(ii)=fName{ii}{max(end-1,1),1};
    end
    groups=['Group';groups];
    d=[d(:,1),groups,d(:,2:end)];
    vName=d(1,2:end); rName=d(2:end,2); d=d(2:end,2:end);
    
    %avoid unique values in row names
    [~,idxu,idxc] = unique(rName);
    %count unique values
    [count, ~, idxcount] = histcounts(idxc,numel(idxu));
    %where is greater than one occurence
    idxkeep = count(idxcount)>1;
    %add extra numbers
    for ii=1:length(idxkeep)
        if idxkeep(ii), rName(ii)=join([rName(ii) num2str(ii)]);
        end
    end
    d=array2table(d,'RowNames',rName,'VariableNames',vName);
    
    %show the new figure.
    f=uifigure('Name','Results','Resize','off','Visible','off');
    f.Position(3:4)=[800 400];
    figure_size=get(f,'outerposition');
    movegui(f,'center');
    
    %add the figure to the right.
    g=uiaxes(f);
    g.Position=[figure_size(3)/2 0 figure_size(3)/2 figure_size(4)-44];
    
    %add the dropdown box to the upper right.
    c=uidropdown(f,'Items',vName(2:end),...
        'ValueChangedFcn',@(c,event) changeGroups());
    c.Position(1:3)=[figure_size(3)/2 figure_size(4)-22 figure_size(3)/2];
    
    %add the table to the left.
    uit=uitable(f,'Data',d(1:end,:),...
            'DisplayDataChangedFcn',@(uit,event) changeGroups());
    uit.ColumnSortable=true;
    edit=zeros(1,width(d)); edit(1)=1;
    uit.ColumnEditable=logical(edit);
    uit.Position=[0 0 figure_size(3)/2 figure_size(4)];
    
    %set up the initial graph.
    drw=[]; resetGraph(d);
    drawnow; figure(f);
    
    function changeGroups()
        d=uit.DisplayData;
        resetGraph(d);
    end
    
    function resetGraph(d)
        %get the selected variable.
        val=c.Value;
        d=d(:,{'Group' val});
        groupLabels=d{:,1};
        unitval=string(uit.Data.Properties.VariableNames)==val;
        
        %seperate the groups.
        [~,~,Index]=unique(d(:,1));
        d(:,1)=table(Index);
        groupLabels=unique(groupLabels);
        
        %convert to double values and calculate.
        for jj=1:2, d.(jj)=str2double(d{:,jj});
        end
        test=grpstats(d,'Group',{'mean' 'std'});
        
        %draw the graph.
        if(~isempty(drw)), delete(drw);
        end
        hold(g,'on');
        drw(1)=bar(g,test.(3));
        drw(2)=errorbar(g,test.(3),test.(4),'k','CapSize',25,'LineStyle','none');
        xticks(g,linspace(1,length(groupLabels),length(groupLabels)));
        xticklabels(g,groupLabels); xtickangle(g,45);
        ylabel(g,tUnits(unitval)); ylim(g,'auto');
        hold(g,'off');
    end
end

function [f,graphs]=graphmaking(fTitle,time,data,peakcontract,peakrelax,thrsh,indval,lowmax,iv,highmin,fps,dist,frames,maxframes,numMov,dt,manualthresh,allvalues,peaklimit,dX,stiffness,A,l,peakPriority,bInvalid)   
    lb=waitbar(0,'Preparing graphs','CloseRequestFcn',[]);
    f=uifigure('Name','Result graphs','Resize','off','Visible','off');
    figure_size=get(f,'outerposition');
    tg=uitabgroup(f,'Position',[0 0 figure_size(3:4)]);
    graphs(numMov)=struct();
    movegui(f,'center');
    
    %add tab per movie.
    for ii=1:numMov
        waitbar(ii/(numMov*2),lb,['Preparing graph ' num2str(ii) ' of ' num2str(numMov)]);
        tabs(ii)=uitab(tg,'Title',fTitle(ii));
        figure_size=get(tabs(ii),'outerposition');
        graphs(ii).fig=uiaxes(tabs(ii),'Position',[0 0 figure_size(3) figure_size(4)-28]);
    end
    figure_size=get(tabs(1),'outerposition');
    
    % Draw the graphs.
    for ii = 1:numMov
        waitbar((numMov+ii)/(numMov*2),lb,['Drawing graph ' num2str(ii) ' of ' num2str(numMov)]);
        graphs(ii).drw=lineplot(graphs(ii).fig,fTitle(ii),data(:,ii),time(1:frames(ii)),peakcontract(:,ii),peakrelax(:,ii),thrsh(ii),indval{ii},lowmax{ii},iv{ii},highmin{ii},frames(ii),peaklimit(ii));
        title(graphs(ii).fig,fTitle(ii));

        % Add an input box for switching peaks.
        swp=uicheckbox(tabs(ii),'Text','Switch peaks',...
        'Value',peakPriority(ii),...
        'Position',[figure_size(3)-120 124 100 22],...
        'Interruptible','on',...%FIXME1: does not work!
        'BusyAction','cancel',...%FIXME1: does not work!
        'ValueChangedFcn',@(swp,event) switchPeaks(ii,swp.Value));
        graphs(ii).swp=swp;
        
        % Add an input box for setting threshold.
        trh=uieditfield(tabs(ii),'Tooltip','Manual threshold',...l
        'Value',num2str(thrsh(ii)),...
        'Position',[figure_size(3)-120 100 60 22],...
        'Interruptible','on',...%FIXME1: does not work!
        'BusyAction','cancel',...%FIXME1: does not work!
        'ValueChangedFcn',@(trh,event) switchThreshold(ii,trh.Value));
        uilabel(tabs(ii),'Position',[figure_size(3)-60 100 60 22],'Text',"Threshold");
        graphs(ii).trh=trh;
        
        % Add an input box for setting stiffness.
        stf=uieditfield(tabs(ii),'Tooltip','Shear Modulus (kPa)',...
        'Value',num2str(stiffness(ii)/(1E3),3),...
        'Position',[figure_size(3)-120 76 60 22],...
        'Interruptible','on',...%FIXME1: does not work!
        'BusyAction','cancel',...%FIXME1: does not work!
        'ValueChangedFcn',@(stf,event) switchStiffness(ii,str2double(stf.Value)*(1E3)));
        uilabel(tabs(ii),'Position',[figure_size(3)-60 76 60 22],'Text',"kPa");
        
        % Add an input box to exclude data.
        inv=uicheckbox(tabs(ii),'Text','Exclude data',...
        'Value',bInvalid(ii),...
        'Position',[figure_size(3)-120 52 100 22],...
        'Interruptible','on',...%FIXME1: does not work!
        'BusyAction','cancel',...%FIXME1: does not work!
        'ValueChangedFcn',@(inv,event) markInvalid(ii,inv.Value));
        graphs(ii).inv=inv;

        % Force draw.
        drawnow limitrate nocallbacks
    end
    drawnow
    figure(f);
    close(lb,'force');
    
    function markInvalid(movie,bInvalid)
        evalin('base',char(['bInvalid(' string(movie) ') = ' string(bInvalid)]));
    end
    
    function switchPeaks(movie,peakPriority)
        evalin('base',char(['peakPriority(' string(movie) ') = ' string(peakPriority)]));
        
        %calculate all values again with the new peaks
        [allvalues(movie,:),peakcontract(:,movie),peakrelax(:,movie),thrsh(movie),indval{movie},lowmax{movie},iv{movie},highmin{movie}] = calculating(data(:,movie),frames(movie),maxframes,fps,dist,dt(:,movie),manualthresh(movie),time(1:frames(movie)),peakPriority);
        
        %Overwrite allvalues in the Workspace
        assignin('base','allvalues',allvalues);
        
        % Create new plots
        delete(graphs(movie).drw);
        graphs(movie).drw=lineplot(graphs(movie).fig,fTitle(movie),data(:,movie),time(1:frames(movie)),peakcontract(:,movie),peakrelax(:,movie),thrsh(movie),indval{movie},lowmax{movie},iv{movie},highmin{movie},frames(movie),peaklimit(movie));
    end
    
    function switchThreshold(movie,mthresh)
        manualthresh(movie)=str2double(mthresh); %correctly put the manual threshold in a variable
        
        % Calculate all values again with the new threshold
        peakPriority=evalin('base',char(['peakPriority(' string(movie) ')']));
        [allvalues(movie,:),peakcontract(:,movie),peakrelax(:,movie),thrsh(movie),indval{movie},lowmax{movie},iv{movie},highmin{movie}] = calculating(data(:,movie),frames(movie),maxframes,fps,dist,dt(:,movie),manualthresh(movie),time(1:frames(movie)),peakPriority);
        % Add the manual threshold in the Workspace
        evalin('base',char(['manualthresh(' string(movie) ') = ' string(manualthresh(movie))]));
        % Overwrite allvalues in the Workspace
        assignin('base','allvalues',allvalues);
        
        % Create new plots
        delete(graphs(movie).drw);
        graphs(movie).drw=lineplot(graphs(movie).fig,fTitle(movie),data(:,movie),time(1:frames(movie)),peakcontract(:,movie),peakrelax(:,movie),thrsh(movie),indval{movie},lowmax{movie},iv{movie},highmin{movie},frames(movie),peaklimit(movie));
    end
    
    function switchStiffness(movie,mstiff)
        mstiff=max(mstiff,0.0001); %we do not allow 0 kPa Shear Modulus.
        stiffchange=evalin('base',char(['stiffness(' string(movie) ')']));
        evalin('base',char(['stiffness(' string(movie) ') = ' string(mstiff)]));
        stiffchange=mstiff/stiffchange;
        
        % Recalculate pressure.
        sumforce(:,movie) = (mstiff*A/l).*dX(:,movie); % average N per second on a single moving pixel
        sumpressure(:,movie) = sumforce(:,movie)/A; % average N/m^2 of a moving pixel per second  
        data(:,movie) = sumpressure(:,movie); % average N/m^2 of a moving pixel per second 
        dt(:,movie) = diff(data(:,movie)); %derivative of (average N/m^2 of a moving pixel per second) over time
        
        % Calculate all values again with the new peaks
        %peakPriority=evalin('base',char(['peakPriority(' string(movie) ')']));
        %[allvalues(movie,:),peakcontract(:,movie),peakrelax(:,movie),thrsh(movie),indval{movie},lowmax{movie},iv{movie},highmin{movie}] = calculating(data(:,movie),frames(movie),maxframes,fps,dist,dt(:,movie),manualthresh(movie),time(1:frames(movie)),peakPriority);
        peaklimit(movie)=max(data(:,movie));
        
        % Overwrite allvalues in the Workspace
        assignin('base','allvalues',allvalues);
        assignin('base','peaklimit',peaklimit);
        
        % Create new plots
        graphs(movie).trh.Value=string(str2double(graphs(movie).trh.Value)*stiffchange);
        switchThreshold(movie,graphs(movie).trh.Value);
        %delete(graphs(movie).drw);
        %graphs(movie).drw=lineplot(graphs(movie).fig,fTitle(movie),data(:,movie),time(1:frames(movie)),peakcontract(:,movie),peakrelax(:,movie),thrsh(movie),indval{movie},lowmax{movie},iv{movie},highmin{movie},frames(movie),peaklimit);
    end
end
