% MibiExtractSingleCellDataFromSegmentation
% assign expression data to a cell based on segmentation
% I use regionprops to speed up. Gives same results as version2, but much
% faster!

massDS = MibiReadMassData(['original.csv']);
%path = '1FOV/';
%pathSegment = '1FOV';
path = 'Inputs/Point1/';
pathSegment = 'Inputs';
resultsDir = [path,'dataPerCell'];
mkdir(resultsDir);
clusterChannels = {'C';'Na';'P';'S';'Fe';'dsDNA';'Vimentin';'Histone H3';'CD16';'SMA';'CD209 (DC-SIGN)';'NFkB-p100 (pS865)';'CD4';'CD11c';'CD56';'FoxP3';'CD39';'Granzyme B';'Biotin';'CD21 (CR2)';'Ki-67';'PD-1';'Pax-5';'CCR7';'CD163';'CD68';'FoxO1';'CD8';'CD3';'CD45-RA';'Dinitrophenyl (DNP)';'CD86';'CTLA-4';'CD20';'Lamin AC';'MPO';'HLA-DR';'IL10';'CD169 (Sialoadhesin)';'CD8a';'Pan-Keratin';'CD11b';'CD36';'Digoxigenin (DIG)';'CD25';'CD45';'Ta';'Au'};
[~, clusterChannelsInds] = ismember(clusterChannels,massDS.Label);

%boundary compensation codes
%selected channels to do the boundary compensation
normChannels = {'dsDNA';'Histone H3';'CD16';'SMA';'CD209 (DC-SIGN)';'NFkB-p100 (pS865)';'CD4';'CD11c';'CD56';'FoxP3';'CD39';'Granzyme B';'Biotin';'CD21 (CR2)';'Ki-67';'PD-1';'Pax-5';'CCR7';'CD163';'CD68';'FoxO1';'CD8';'CD3';'CD45-RA';'Dinitrophenyl (DNP)';'CD86';'CTLA-4';'CD20';'Lamin AC';'MPO';'HLA-DR';'IL10';'CD169 (Sialoadhesin)';'CD8a';'Pan-Keratin';'CD11b';'CD36';'Digoxigenin (DIG)';'CD25';'CD45'};
[~, normChannelsInds] = ismember(normChannels,massDS.Label);
channelNormIdentity = zeros(length(clusterChannels),1);
for i = 1:length(normChannelsInds)
    channelNormIdentity(normChannelsInds(i)) = 1;
end

for p=1:1
    disp(['point',num2str(p)]);
    pointNumber=p;
    % load data and get nuclear markers
    %load(['segmentated_DNA_1FOV_resize_back/Point1_0.01_0.35/segmentationParams.mat']);
    load(['segmentationParams.mat']);
    countsNoNoiseRec = zeros(size(newLmod,1),size(newLmod,2),length(clusterChannels));
    
    for i=1:length(clusterChannels)
        tiffDir = strcat(path,clusterChannels(i),'.tiff');
        t = Tiff(tiffDir{1},'r');
        nucImTemp = read(t);
        countsNoNoiseRec(:,:,i) = nucImTemp;
    end
    save([path,'/dataDeNoiseCohort.mat'],'countsNoNoiseRec');
    
    % original 
    load([path,'/dataDeNoiseCohort.mat']);
    
    labelNum = max(max(newLmod));
    channelNum = length(massDS);
    stats = regionprops(newLmod,'Area','PixelIdxList','Centroid');
    %
    countsReshape= reshape(countsNoNoiseRec,size(countsNoNoiseRec,1)*size(countsNoNoiseRec,2),channelNum);
    
    % make a data matrix the size of the number of labels x the number of markers
    data = zeros(labelNum,channelNum);
    dataScaleSize = zeros(labelNum,channelNum);
    cellSizes = zeros(labelNum,1);

    % for each label extract information
    for i=1:labelNum
        %
        currData = countsReshape(stats(i).PixelIdxList,:);
        
        data(i,:) = sum(currData,1);
        dataScaleSize(i,:) = sum(currData,1) / stats(i).Area;
        cellSizes(i) = stats(i).Area;
    end

    % get the final information only for the labels with 
    % 1.positive nuclear identity (cells)
    % 2. That have enough information in the clustering channels to be
    % clustered
    labelIdentityNew2 = labelIdentityNew([1:end-1]); % fix bug resulting from previous script
    sumDataScaleSizeInClusterChannels = sum(dataScaleSize(:,clusterChannelsInds),2);
    labelIdentityNew2(sumDataScaleSizeInClusterChannels<0.1) = 2;
    dataCells = data(labelIdentityNew2==1,:);
    dataScaleSizeCells = dataScaleSize(labelIdentityNew2==1,:);
    labelVec=find(labelIdentityNew2==1);
    
    % get cell sizes only for cells
    cellSizesVec = cellSizes(labelIdentityNew2==1);

    % asinh transform
    dataScaleSizeCellsTrans = asinh(dataScaleSizeCells);
    dataCellsTrans = asinh(dataCells);
    
%     % standardize
%     dataScaleSizeCellsTransStd = zscore(dataScaleSizeCellsTrans);
%     dataCellsTransStd = zscore(dataCellsTrans);
%     dataCellsStd = zscore(dataCells);
%     dataScaleSizeCellsStd = zscore(dataScaleSizeCells);

    dataTransL = [labelVec,cellSizesVec,dataCellsTrans];
    dataScaleSizeTransL = [labelVec,cellSizesVec, dataScaleSizeCellsTrans];
    dataL = [labelVec,cellSizesVec, dataCells];
    dataScaleSizeL = [labelVec,cellSizesVec, dataScaleSizeCells];
%     dataTransStdL = [labelVec,dataCellsTransStd];
%     dataScaleSizeTransStdL = [labelVec,dataScaleSizeCellsTransStd];
%     dataStdL = [labelVec,dataCellsStd];
%     dataScaleSizeStdL = [labelVec,dataScaleSizeCellsStd];

    channelLabelsForFCS = ['cellLabelInImage';'cellSize';massDS.Label];

    %% save fcs
    TEXT.PnS = channelLabelsForFCS;
    TEXT.PnN = channelLabelsForFCS;
    save([pathSegment,'/cellData.mat'],'labelIdentityNew2','labelVec','cellSizesVec','dataCells','dataScaleSizeCells','dataScaleSizeCellsTrans','dataCellsTrans','channelLabelsForFCS');
    writeFCS([pathSegment,'/dataFCS.fcs'],dataL,TEXT);
    writeFCS([pathSegment,'/dataScaleSizeFCS.fcs'],dataScaleSizeL,TEXT);
    writeFCS([pathSegment,'/dataTransFCS.fcs'],dataTransL,TEXT);
    writeFCS([pathSegment,'/dataScaleSizeTransFCS.fcs'],dataScaleSizeTransL,TEXT);
%     writeFCS([path,'/Point',num2str(pointNumber),'/dataStdFCS.fcs'],dataStdL,TEXT);
%     writeFCS([path,'/Point',num2str(pointNumber),'/dataScaleSizeStdFCS.fcs'],dataScaleSizeStdL,TEXT);
%     writeFCS([path,'/Point',num2str(pointNumber),'/dataTransStdFCS.fcs'],dataTransStdL,TEXT);
%     writeFCS([path,'/Point',num2str(pointNumber),'/dataScaleSizeTransStdFCS.fcs'],dataScaleSizeTransStdL,TEXT);
    
    writeFCS([resultsDir,'/dataFCS_p',num2str(pointNumber),'.fcs'],dataL,TEXT);
    writeFCS([resultsDir,'/dataScaleSizeFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeL,TEXT);
    writeFCS([resultsDir,'/dataTransFCS_p',num2str(pointNumber),'.fcs'],dataTransL,TEXT);
    writeFCS([resultsDir,'/dataScaleSizeTransFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeTransL,TEXT);
%     writeFCS([resultsDir,'/dataStdFCS_p',num2str(pointNumber),'.fcs'],dataStdL,TEXT);
%     writeFCS([resultsDir,'/dataScaleSizeStdFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeStdL,TEXT);
%     writeFCS([resultsDir,'/dataTransStdFCS_p',num2str(pointNumber),'.fcs'],dataTransStdL,TEXT);
%     writeFCS([resultsDir,'/dataScaleSizeTransStdFCS_p',num2str(pointNumber),'.fcs'],dataScaleSizeTransStdL,TEXT);
    
end