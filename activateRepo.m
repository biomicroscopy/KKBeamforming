
currentDir = pwd;

subFolders = {"Datasets","inc","EigenMEX","SubFunctions"};

for i = 1:length(subFolders)
    addpath(genpath(fullfile(currentDir,subFolders{i})));
end

clearvars