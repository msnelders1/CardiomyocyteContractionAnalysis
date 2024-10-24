function fileOut=parseFiles(folder,fileIn)
    fileOut=fileIn;
    %check subfolders first.
    subfolder=dir(convertCharsToStrings(folder));
    subfolder(~[subfolder.isdir])=[]; %get rid of non folders
    subfolder(ismember({subfolder.name},{'.', '..'})) = []; %get rid of . and ..
    subfolder=fullfile(folder,{subfolder.name}); %construct list of names
    %add subfolders to the list.
    for ii=1:length(subfolder)
        fileOut=parseFiles(subfolder(ii),fileOut);
    end
    %add the current folder's files to the list.
    if(~isempty(fileOut))
        fileOut=[fileOut;dir(fullfile(convertCharsToStrings(folder),'*_processed.csv'))];
    else, fileOut=dir(fullfile(convertCharsToStrings(folder),'*_processed.csv'));
    end
end
