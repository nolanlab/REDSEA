

<p align="left"><img width=20%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/redsea.jpg"></p>

# REDSEA
We present **RE**inforcement **D**ynamic **S**pillover **E**limin**A**tion (REDSEA) as a solution for spillover compensation without loss of signal. This is the step-by-step guidance of how to use REDSEA to produced a FCS file with compensated channel information.

<p align="center"><img width=40%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/overview.png"></p>

### How it works

In breif, the REDSEA algorithm identifies the boudary region for each cell based on the segmentation mask. Subsequently, for one cell's each channel, signals were subtracted based on the shared boundary with neighbooring cells ,and their corresponding signal.  Moreover, the removed signal from adjacent cells can be reinforced back to the cell by option. For detailed information look at the Material and Methods section of the [paper](www.facebook.com).

## Table of content

- [Required Inputs](#required-inputs)
    - [MATLAB scripts](#matlab-scripts)
    - [tiff images and channel information](#tiff-images-and-channel-information)
    - [Segmentation Mask](#segmentation-mask)

- [REDSEA](#redsea)
    - [Parameters](#parameters)
    - [Sanity plots](#sanity-plots)
    - [Output](#output)

## Required Inputs


### MATLAB scripts
The [MATLAB scripts](https://github.com/BokaiZhu/REDSEA/tree/master/code) include the REDSEA algorithm, along with other scripts created by [Leeat Keren](https://github.com/lkeren/MIBIAnalysis) for processing multiplexed images.

### tiff images and channel information
There should be a CSV file containing the channel name information.

<p align="left"><img width=20%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/csv_example.png"></p>

And also a folder containing .tiff images for each channel, where the names should be same as in the CSV file.

### Segmentation Mask
This method requires you to provide a cell segmentation mask. A cell nuclei probability matrix can be produced by your own chose (popular options include [ilastik](https://www.ilastik.org/) or [deepcell](https://github.com/vanvalenlab/deepcell-tf)). In our case for the MIBI data we implemented a trained-in-house deepcell CNN model.

<p align="center"><img width=50%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/probability_matrix.png"></p>

DeepCell is also easy to implement on different imaging modalities. Here is a prediction model we trained with ~ 1500 cells in the [cycIF](https://www.nature.com/articles/s41597-019-0332-y) dataset:

<p align="center"><img width=70%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/cycIF.png"></p>


After producing a nuclei probablity mask, we will then use the script ```MibiSegmentByDeepProbWithPerim3.m``` to implement a watershed algorithm for whole cell segmentation. This will produce something like this:

<p align="center"><img width=50%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/watershed-result.png"></p>

<details><summary>See MibiSegmentByDeepProbWithPerim3.m Script</summary>
<p>



```MATLAB
%% Pipeline for nuclear segmentation using pixel probabilities from deepCell
% Changed the pipeline to use the boundaries from the deepcell as boudaries
% instead of the raw nuclear intensities.

t1 = clock;
% maxs = imextendedmax(probNuc,0.015);
% probNuc>0.05

path = 'inputs/Point1/dsDNA.tiff';
%path = '1FOV/Histone H3.tiff';
resultsPath = 'watershed_result';
%resultsPath = 'segmentated_H3_1FOV_resize_back';
%deepPath = 'deepcell_1FOV';
deepPath = 'deepcell';

%segmentThres = 0.015; %change from 0.1 to 0.02
%probNucThres = 0.5; %change from 0.05 to 0.5
for segmentThres=0.01:0.005:0.02
for probNucThres=0.05:0.15:0.5
for p=1:1
    disp(['point',num2str(p)]);
    pointNumber=p;
    
    %{
    %% Get perimiters of nuclei for plotting
    % load data and get nuclear markers
    load([path,'/Point',num2str(pointNumber),'/dataDeNoiseCohort.mat']);
    imSize = (size(countsNoNoise,1));
    % sum nuclear markers to increase contrast
    nucleiChannels = {'Histone H3'};
    [tf loc] = ismember(nucleiChannels,massDS.Label);
    nucIm = sum(countsNoNoise(:,:,loc),3);
    %}
    
    % read tiff image back to matrix for prediction
    t = Tiff(path,'r');
    nucIm = read(t);
    maxv=50;
    rgb_image = MibiGetRGBimageFromMat(nucIm,maxv);

    % %% Get maxima from deep learning probabilities
    % read nuclear segmentation from Ilastik/deepcell
    probNuc = double(imread([deepPath,'/nuclei-probability-matrix.tif']));
    figure;
    imagesc(probNuc);

    % find local maxima in probability map
    maxs = imextendedmax(probNuc,segmentThres); % change from 0.1 to 0.02
    rgb_image_perim_extMax = imoverlay(rgb_image , maxs, [1 0 0]);   
    figure;
    imagesc(rgb_image_perim_extMax);
    
    %% watershed over the deep results
    bw1 = zeros(size(probNuc));
    bw1(probNuc>probNucThres) = 1; %change from 0.05 to 0.5
    figure; imagesc(bw1);
    bw2 = bwareaopen(bw1,40);
    SE = strel('disk',4);
    bw=imdilate(bw2,SE);
    
    [B,L] = bwboundaries(maxs,4,'noholes');
    maxsFix = bw & maxs;

    % modify the image so that the background pixels and the extended maxima pixels are forced to be the only local minima in the image.
    Jc = imcomplement(probNuc);
    I_mod = imimposemin(Jc, ~bw | maxsFix);
    L = watershed(I_mod);
    labeledImage = label2rgb(L);
    
    cellPerimNewMod= L;
    cellPerimNewMod(L>0) = 100;
    cellPerimNewMod(cellPerimNewMod==0)=1;
    cellPerimNewMod(cellPerimNewMod==100)=0;
    
    rgb_image_cellPerim = imoverlay(rgb_image, cellPerimNewMod, [1 0 0]);
    figure;
    imagesc(rgb_image_cellPerim);

    %% 1. For each label, decide whether it is of a nucleus/ background
    t=40;
    labelNum = length(unique(L(:)));
    labelIdentity = zeros (labelNum,1);
    labelPixelsPercentInNucleiMask = zeros(labelNum,1);
    labelSize = zeros(labelNum,1);
    
    %spend a lot of time during this loop 
    for i=1:labelNum
        [r c] = find(L==i);
        labelSize(i) = length(r); 
        labelMask = (L==i);
        labelPixelsNumInNucleiMask = sum(bw(labelMask));
        labelPixelsPercentInNucleiMask(i) = labelPixelsNumInNucleiMask / labelSize(i);
        if (labelPixelsPercentInNucleiMask(i) > 0.7)
            labelIdentity(i) = 1;
        end
    end

    % 2. Merge small regions within the nuclei mask with their neighbours
    keepVec = ones(labelNum,1);
    newL = L;
    for i=1:labelNum
        if (labelIdentity(i) == 1) && (labelSize(i) < t)
            disp(['Removing label ',num2str(i),'. Size: ',num2str(labelSize(i))]);
            % get neighbour with largest border that is also in nuclear region
            [neighbourLabels , neighbouringRegionSize] = MibiGetNeighbourLabels (newL, i);
            found = 0;
            [neighbouringRegionSizeS , neighbouringRegionSizeSInd] = sort(neighbouringRegionSize,'descend');
            neighbourLabelsS = neighbourLabels(neighbouringRegionSizeSInd);
            maxInd = 1;
            while ~found
                mergeLabelId = neighbourLabelsS(maxInd);
                if (~(mergeLabelId == 0) && (labelIdentity(mergeLabelId) == 1))
                    found = 1;
                else
                    maxInd = maxInd+1;
                end
                if (maxInd >length(neighbourLabelsS)) % reached end of neighbours with no good merging candidate
                    disp (['Warning: no good merging target found for label', num2str(i), '. Keeping it.']);
                    break;
                end
            end
            % update
            if (maxInd <= length(neighbourLabelsS))
                [newL] = MibiMergeLabels (newL, i, mergeLabelId);
                keepVec(i) = 0;
            end
        end
    end

    % Update label numbers to account for deleted labels
    allLabels = [1:labelNum];
    currLabels =  allLabels(keepVec == 1);
    labelIdentityNew = zeros(length(currLabels),1); 
    newLmod = newL;
    for i = 1:length(currLabels)
        newLmod(newLmod == currLabels(i)) = i;
        labelIdentityNew(i) = labelIdentity(currLabels(i));
    end

    cellPerimNewMod= bwperim(newLmod);
    cellPerimNewMod= newLmod;
    cellPerimNewMod(newL>0) = 100;
    cellPerimNewMod(cellPerimNewMod==0)=1;
    cellPerimNewMod(cellPerimNewMod==100)=0;
    rgb_image_cellPerimNewMod = imoverlay(rgb_image, cellPerimNewMod, [1 0 0]);
    % add ilastik segmentation
    % ilsegParams = load([path,'/SegmentPerim/Point',num2str(pointNumber),'/segmentationParams.mat']);
    % rgb_image_cellPerimNewModIl = imoverlay(rgb_image, ilsegParams.cellPerimNewMod, [0 0 0]);
    % rgb_image_cellPerimNewModBoth = imoverlay(rgb_image_cellPerimNewModIl, cellPerimNewMod, [1 0 0]);
    
    figure;
    imagesc(rgb_image_cellPerimNewMod);
    mkdir([resultsPath,'/Point',num2str(pointNumber),'_',num2str(segmentThres),'_',num2str(probNucThres)]);
    imwrite(rgb_image_cellPerimNewMod,[resultsPath,'/Point',num2str(pointNumber),'_',num2str(segmentThres),'_',num2str(probNucThres),'/compareIlastikDeep.tif'],'tif');
    imwrite(cellPerimNewMod,[resultsPath,'/Point',num2str(pointNumber),'_',num2str(segmentThres),'_',num2str(probNucThres),'/PureSegmentation.tif'],'tif');
    save([resultsPath,'/Point',num2str(pointNumber),'_',num2str(segmentThres),'_',num2str(probNucThres),'/segmentationParams.mat'],'newLmod','cellPerimNewMod','labelIdentityNew');
    %close all;
end
end
end

t2=clock;
t = etime(t2,t1);
disp(['elapsed time: ',num2str(t)]);

```

</p>
</details>

## REDSEA

### Parameters
With the .tiff images, .mat segmentation file and .csv channel information, we are now ready to implement REDSEA for boundary compensation.

There are **two methods** for boundary compensation in REDSEA: ```Sudoku``` and ```Cross```. The algorithm walks through the boundaries of each cell, and decides the area to extract signal. You would need to choose one of the two methods and deside how many pixels to expand from the boundary pixel: 

<p align="center"><img width=40%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/method_show.png"></p>

The pixel number for expansion should be **proportional** to the cell size: in our MIBI data, the average cell size is 107 pixels, and we used 2 pixels for expansion; In the [cycIF tonsil dataset](https://www.synapse.org/#!Synapse:syn17796423) the average cell size is 325 pixels, so 3-4 pixels is recommended. 

Also, you need to supply a list of channel names to perform the compensation process: normally you should only compensate for the **surface markers**, like in our case:
```'CD16';'CD209 (DC-SIGN)';'CD4';'CD11c';'CD56';'CD39';'CD21 (CR2)';'PD-1';'CCR7';'CD163';'CD68';'CD8';'CD3';'CD45-RA';'CD86';'CTLA-4';'CD20';'MPO';'HLA-DR';'CD169 (Sialoadhesin)';'CD8a';'CD11b';'CD36';'CD25';'CD45'```

Take a look at the annotated code in the block under:

<details><summary>See MibiExtractSingleCellDataFromSegmentationAndTiff_REDSEA.m Script</summary>
<p>



```MATLAB
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
% elementShape, can be selected from 1-4.
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

% Whether what to plot scatter to check the REDSEA result and effect,
% default=0 for not, 1 for plotting.
plotSanityPlots = 1;
pathSanityPlots = strcat('result/sanityPlots/', 'Shape',num2str(elementShape), 'elementSize', num2str(elementSize), '/');

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
%     load([path,'/Point',num2str(pointNumber),'/segmentationParams.mat']);
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

    %% plot sanity scatter images
    if plotSanityPlots == 1
        mkdir(pathSanityPlots);
        MIBIboundary_compensation_plotting(dataScaleSizeCells,dataCompenScaleSizeCells,normChannels,normChannelsInds,pathSanityPlots);
    end    
    
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
```

</p>
</details>

### Sanity plots

To visually inspect if the parameters discribed in the previous section are optimal for your data, we provided with a flag in the script ```plotSanityPlots = 1```. Once set to 1, it will produce the pairwise combination scatter plots of all the compensated channels. User can use these plots to evaluate if the compensation is optimal (for example looking at CD4-CD8, CD3-CD20 etc.) It is suggested to run this first will **less** channels, find the optimal parameters, then run the full list.  




### Output

REDSEA will produce the 4 FCS files for downstream analysis:

```dataFCS.fcs``` = raw counts for each single cell without compensation

```dataRedSeaFCS``` = raw counts for each single cell with REDSEA compensation

```dataScaleSizeFCS.fcs``` = counts for each single cell as scaled by cell size (counts/cellSize) without compensation

```dataRedSeaScaleSizeFCS.fcs``` = counts for each single cell as scaled by cell size (counts/cellSize) with REDSEA compensation

It is recommended to use the ```dataRedSeaScaleSizeFCS.fcs``` file. We can see by using REDSEA compensation, the boundary signal spillover is dynamically eliminated and reinforced.

For the MIBI dataset:
<p align="center"><img width=55%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/fcs_result.png"></p>

And REDSEA also works comparably on the cycIF data:
<p align="center"><img width=57%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/cycIF_fcs.png"></p>


