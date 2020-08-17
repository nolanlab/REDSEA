%% Extract single cell data and do the REDSEA compensation
% Based on the original script from Leeat Keren
% This now reads in any folder of TIFFs with individual channels from the 
% same field of view, and uses that to recreate a countsNoNoise based on 
% the massDS order.
% Then performs REDSEA compensation as implemented by Yunhao Bai
% The outputs will be REDSEA compensated and non-compensated FCS files
% 4May2020, Yunhao Bai, Sizun Jiang


% Main path for the all the data
mainPath = 'sampleData_MIBI';

% This is a csv file for your channel labels within
massDS = dataset('File',[mainPath,'/sampleData.csv'],'Delimiter',',');

% This assumes the path points to a folder containing all the Points from 
% the run. Your segmentationParams.mat from each point should be in the 
% each Point's folder
pathTiff = [mainPath,'/originalTiff']; 

% This is where the FCS file output will go to
pathResults = [mainPath,'/FCS_output'];

% Select the channels that are expected to be expressed. Cells with minimal
% expression of at least one of these channels will be removed
clusterChannels = massDS.Label; % exclude elemental channels (here the 
% sample dataset does not contain 
[~, clusterChannelsInds] = ismember(clusterChannels,massDS.Label);

% boundaryMod determines the type of compensation done.
% 1:whole cell compensation
% 2:boundary compensation (default)
boundaryMod = 2;
% REDSEAChecker determines the type of compensation done.
% 0:only subtraction; 
% 1:subtraction and reinforcement (default)
REDSEAChecker = 1;

% for boundary compensation, needs to specify elementShape. 
% 1:Sudoku style, 2:Cross style
elementShape = 2;
% elementSize. How many pixels around the center to be considered for the
% elementShape, can be selected from 1-4.
% As a default, keep elementShape and elementSize as 2.
elementSize = 2;

% Select channels for REDSEA compensation. Surface markers are recommended
% boundary compensation codes
% selected channels to do the boundary compensation
normChannels = {'CD4';'CD56';'CD21 (CR2)';'CD163';'CD68';'CD3';'CD20';'CD8a'};
%normChannels = {'x7500y3500_1700_DAPI';'x7500y3500_1700_CD3';'x7500y3500_1700_CD4';'x7500y3500_1700_CD8a';'x7500y3500_1700_CD11b';'x7500y3500_1700_CD20';'x7500y3500_1700_CD45';'x7500y3500_1700_CD68'};
[~, normChannelsInds] = ismember(normChannels,massDS.Label);
channelNormIdentity = zeros(length(massDS.Label),1);
% Getting an array of flags for whether to compensate or not
for i = 1:length(normChannelsInds)
    channelNormIdentity(normChannelsInds(i)) = 1;
end

% Whether what to plot scatter to check the REDSEA result and effect,
% default=0 for not, 1 for plotting.
% Note that if multiple channels selected (in normChannels), to plot out all
% the sanity plots need long time.
plotSanityPlots = 0;

%%
mkdir(pathResults);

for p=1:1
    disp(['point',num2str(p)]);
    pointNumber = p;
    % load tiffs to recreate countsNoNoise
    for i=1:length(massDS.Label)
        t = imread([pathTiff, '/Point', num2str(pointNumber), '/', massDS.Label{i}, '.tiff']);
        d = double(t);
        % imshow(d)
        countsNoNoise(:,:,i) = d;
    end
        
    % load segmentation file
    load([pathTiff, '/Point', num2str(pointNumber), '/segmentationParams.mat']);
    labelNum = max(max(newLmod));
    channelNum = length(massDS);
    stats = regionprops(newLmod,'Area','PixelIdxList'); % Stats on cell size. Region props is DF with cell location by count
    countsReshape= reshape(countsNoNoise,size(countsNoNoise,1)*size(countsNoNoise,2),channelNum);
    
    % make a data matrix the size of the number of labels x the number of markers
    % Include one more marker for cell size
    data = zeros(labelNum,channelNum);
    dataScaleSize = zeros(labelNum,channelNum);
    cellSizes = zeros(labelNum,1);
    
    % for each label extract information
    for i=1:labelNum
        currData = countsReshape(stats(i).PixelIdxList,:);
        data(i,1:channelNum) = sum(currData,1);
        dataScaleSize(i,1:channelNum) = sum(currData,1) / stats(i).Area;
        cellSizes(i) = stats(i).Area;
    end
    
    %% do cell boundary compensation
    if boundaryMod == 1
        dataCompen = MIBIboundary_compensation_wholeCellSA(newLmod,data,channelNormIdentity,REDSEAChecker);
    elseif boundaryMod == 2
        dataCompen = MIBIboundary_compensation_boundarySA(newLmod,data,countsNoNoise,channelNormIdentity,elementShape,elementSize,REDSEAChecker);
    end
    dataCompenScaleSize = dataCompen./repmat(cellSizes,[1 channelNum]);

    % %% Add point number 
    % pointnum = double(repmat(p, 1, length(data))); 
    % pointnum = repelem(p, length(data),[1]);
    % data = [data, pointnum];
    % dataScaleSize = [dataScaleSize, pointnum];
    % dataCompen = [dataCompen, pointnum];
    % dataCompenScaleSize = [dataCompenScaleSize, pointnum];

    %%
    % get the final information only for the labels with 
    % 1.positive nuclear identity (cells)
    % 2.that have enough information in the clustering channels to be
    % clustered
    labelIdentityNew2 = labelIdentityNew([1:end-1]); % fix bug resulting from previous script
    sumDataScaleSizeInClusterChannels = sum(dataScaleSize(:,clusterChannelsInds),2);
    labelIdentityNew2(sumDataScaleSizeInClusterChannels<0.1) = 2;
    
    dataCells = data(labelIdentityNew2==1,:);
    dataScaleSizeCells = dataScaleSize(labelIdentityNew2==1,:);
    dataCompenCells = dataCompen(labelIdentityNew2==1,:);
    dataCompenScaleSizeCells = dataCompenScaleSize(labelIdentityNew2==1,:);

    labelVec=find(labelIdentityNew2==1);
    
    % get cell sizes only for cells
    cellSizesVec = cellSizes(labelIdentityNew2==1);

    dataL = [labelVec,cellSizesVec,dataCells,repmat(p,[length(labelVec) 1])];
    dataScaleSizeL = [labelVec,cellSizesVec,dataScaleSizeCells,repmat(p,[length(labelVec) 1])];
    dataCompenL = [labelVec,cellSizesVec,dataCompenCells,repmat(p,[length(labelVec) 1])];
    dataCompenScaleSizeL = [labelVec,cellSizesVec,dataCompenScaleSizeCells,repmat(p,[length(labelVec) 1])];

    % dataTransStdL = [labelVec,dataCellsTransStd];
    % dataScaleSizeTransStdL = [labelVec,dataScaleSizeCellsTransStd];
    % dataStdL = [labelVec,dataCellsStd];
    % dataScaleSizeStdL = [labelVec,dataScaleSizeCellsStd];

    % channelLabelsForFCS = ['cellLabelInImage';'cellSize';massDS.Label];
    channelLabelsForFCS = ['cellLabelInImage';'cellSize';massDS.Label;'PointNum'];
    
    %% output
    outputPath = [pathResults,'/Point',num2str(pointNumber),'/BM=',num2str(boundaryMod),'_RC=',num2str(REDSEAChecker),'_Shape=',num2str(elementShape),'_Size=',num2str(elementSize)];
    mkdir(outputPath);

    % plot sanity scatter images
    if plotSanityPlots == 1
        pathSanityPlots = [outputPath,'/sanityPlots/'];
        mkdir(pathSanityPlots);
        MIBIboundary_compensation_plotting(dataScaleSizeCells,dataCompenScaleSizeCells,normChannels,normChannelsInds,pathSanityPlots);
    end    
    
    % save fcs
    TEXT.PnS = channelLabelsForFCS;
    TEXT.PnN = channelLabelsForFCS;
    
    save([outputPath,'/cellData.mat'],'labelIdentityNew2','labelVec','cellSizesVec','dataCells','dataScaleSizeCells','dataCompenCells','dataCompenScaleSizeCells','channelLabelsForFCS');
    writeFCS([outputPath,'/dataFCS.fcs'],dataL,TEXT);
    writeFCS([outputPath,'/dataScaleSizeFCS.fcs'],dataScaleSizeL,TEXT);
    writeFCS([outputPath,'/dataRedSeaFCS.fcs'],dataCompenL,TEXT);
    writeFCS([outputPath,'/dataRedSeaScaleSizeFCS.fcs'],dataCompenScaleSizeL,TEXT);

    % writeFCS([outputPath,'/dataTransFCS.fcs'],dataTransL,TEXT);
    % writeFCS([outputPath,'/dataScaleSizeTransFCS.fcs'],dataScaleSizeTransL,TEXT);
    % writeFCS([outputPath,'/dataStdFCS.fcs'],dataStdL,TEXT);
    % writeFCS([outputPath,'/dataScaleSizeStdFCS.fcs'],dataScaleSizeStdL,TEXT);
    % writeFCS([outputPath,'/dataTransStdFCS.fcs'],dataTransStdL,TEXT);
    % writeFCS([outputPath,'/dataScaleSizeTransStdFCS.fcs'],dataScaleSizeTransStdL,TEXT);
end