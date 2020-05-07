% MibiExtractSingleCellDataFromSegmentation
% Based on the Original Script by Leeat Keren
% This now reads in any folder of TIFFs with individual channels from the 
% same field of view, and uses that to recreate a countsNoNoise based on 
% the massDS order.
% Then performs REDSEA compensation as implemented by Yunhao Bai
% The outputs will be REDSEA compensated and non-compensated FCS files
% 4May2020 Yunhao Bai, Sizun Jiang

% This is a csv file for your channels within
massDS = MibiReadMassData('example_channel_inforamtion.csv');
path = 'Inputs'; % This assumes the path points to a folder 
% containing all the Points from the run. Your segmentationParams.mat from 
% each point should be in the each Point's folder
% There should be a segmentation.mat in the same folder containing the
% segmentation of the image.

% This is where the FCS file output will go to
pathSegment = 'result/';

% Select the channels that are expected to be expressed. Cells with minimal
% expression of at least one of these channels will be removed
clusterChannels = massDS.Label(6:46); % exclude elemental channels
[~, clusterChannelsInds] = ismember(clusterChannels,massDS.Label);

% boundaryMod determines the type of compensation done for REDSEA.
% elementShape. 1:Sudoku style, 2: Cross style
% elementSize. How many pixels around the center to be considered for the
% elementShape
% As a default, keep elementShape and elementSize as 2.
elementShape = 2;
elementSize = 2;
% Select channels for REDSEA compensation. Surface markers are recommended
%boundary compensation codes
%selected channels to do the boundary compensation
normChannels = {'CD16';'CD209 (DC-SIGN)';'CD4';'CD11c';'CD56';'CD39';'CD21 (CR2)';'PD-1';'CCR7';'CD163';'CD68';'CD8';'CD3';'CD45-RA';'CD86';'CTLA-4';'CD20';'MPO';'HLA-DR';'CD169 (Sialoadhesin)';'CD8a';'CD11b';'CD36';'Digoxigenin (DIG)';'CD25';'CD45'};
[~, normChannelsInds] = ismember(normChannels,massDS.Label);
channelNormIdentity = zeros(length(massDS.Label),1);
% Getting an array of flags for whether to compensate or not
for i = 1:length(normChannelsInds)
    channelNormIdentity(normChannelsInds(i)) = 1;
end

%%
mkdir(pathSegment);

for p=1:1
    disp(['point',num2str(p)]);
    pointNumber=p;
    % Load tiffs to recreate countsNoNoise
    for i=1:length(massDS.Label)
        t = imread([path, '/Point', num2str(pointNumber), '/', massDS.Label{i}, '.tiff']);
        d = double(t);
    %     imshow(d)
        countsNoNoise(:,:,i) = d;
    end
        
    % Load segmentation file
    load([path, '/Point', num2str(pointNumber), '/watershed_result/Point1_0.01_0.35/segmentationParams.mat']);
    labelNum = max(max(newLmod));
    channelNum = length(massDS);
    stats = regionprops(newLmod,'Area','PixelIdxList'); % Stats on cell size. Region props is DF with cell location by count
    countsReshape= reshape(countsNoNoise,size(countsNoNoise,1)*size(countsNoNoise,2),channelNum);
%     % make a data matrix the size of the number of labels x the number of markers
%     data = zeros(labelNum,channelNum);
%     dataScaleSize = zeros(labelNum,channelNum);
%     cellSizes = zeros(labelNum,1);

    % make a data matrix the size of the number of labels x the number of markers
    % Include one more marker for cell size
    data = zeros(labelNum,channelNum);
    dataScaleSize = zeros(labelNum,channelNum);
    cellSizes = zeros(labelNum,1);
    
    % for each label extract information
    for i=1:labelNum
        %
        currData = countsReshape(stats(i).PixelIdxList,:);
        
        data(i,1:channelNum) = sum(currData,1);
        dataScaleSize(i,1:channelNum) = sum(currData,1) / stats(i).Area;
        cellSizes(i) = stats(i).Area;
    end
    
    %% do cell boundary compensation
    dataCompen = MIBIboundary_compensation_boundarySA(newLmod,data,countsNoNoise,channelNormIdentity,elementShape,elementSize);
    
    dataCompenScaleSize = dataCompen./repmat(cellSizes,[1 channelNum]);

%     %% Add point number 
% %     pointnum = double(repmat(p, 1, length(data))); 
%     pointnum = repelem(p, length(data),[1]);
%     data = [data, pointnum];
%     dataScaleSize = [dataScaleSize, pointnum];
%     dataCompen = [dataCompen, pointnum];
%     dataCompenScaleSize = [dataCompenScaleSize, pointnum];
% 
%     
    %%

    % get the final information only for the labels with 
    % 1.positive nuclear identity (cells)
    % 2. That have enough information in the clustering channels to be
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

%     dataTransStdL = [labelVec,dataCellsTransStd];
%     dataScaleSizeTransStdL = [labelVec,dataScaleSizeCellsTransStd];
%     dataStdL = [labelVec,dataCellsStd];
%     dataScaleSizeStdL = [labelVec,dataScaleSizeCellsStd];

%    channelLabelsForFCS = ['cellLabelInImage';'cellSize';massDS.Label];
    channelLabelsForFCS = ['cellLabelInImage';'cellSize';massDS.Label;'PointNum'];

    
    
    %% save fcs
    TEXT.PnS = channelLabelsForFCS;
    TEXT.PnN = channelLabelsForFCS;
    mkdir([pathSegment,'/Point',num2str(pointNumber)]);
    save([pathSegment,'/Point',num2str(pointNumber),'/cellData.mat'],'labelIdentityNew2','labelVec','cellSizesVec','dataCells','dataScaleSizeCells','dataCompenCells','dataCompenScaleSizeCells','channelLabelsForFCS');
    writeFCS([pathSegment,'/Point',num2str(pointNumber),'/dataFCS.fcs'],dataL,TEXT);
    writeFCS([pathSegment,'/Point',num2str(pointNumber),'/dataScaleSizeFCS.fcs'],dataScaleSizeL,TEXT);
    writeFCS([pathSegment,'/Point',num2str(pointNumber),'/dataRedSeaFCS.fcs'],dataCompenL,TEXT);
    writeFCS([pathSegment,'/Point',num2str(pointNumber),'/dataRedSeaScaleSizeFCS.fcs'],dataCompenScaleSizeL,TEXT);


    % writeFCS([resultsDir,'/dataFCS_p',num2str(pointNumber),'.fcs'],dataL,TEXT);
    % writeFCS([resultsDir,'/dataScaleSizeFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeL,TEXT);
    % writeFCS([resultsDir,'/dataTransFCS_p',num2str(pointNumber),'.fcs'],dataTransL,TEXT);
    % writeFCS([resultsDir,'/dataScaleSizeTransFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeTransL,TEXT);
    
%     writeFCS([resultsDir,'/dataStdFCS_p',num2str(pointNumber),'.fcs'],dataStdL,TEXT);
%     writeFCS([resultsDir,'/dataScaleSizeStdFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeStdL,TEXT);
%     writeFCS([resultsDir,'/dataTransStdFCS_p',num2str(pointNumber),'.fcs'],dataTransStdL,TEXT);
%     writeFCS([resultsDir,'/dataScaleSizeTransStdFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeTransStdL,TEXT);
    
end