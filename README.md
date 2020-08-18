

<p align="left"><img width=20%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/redsea.jpg"></p>

# REDSEA
We present **RE**inforcement **D**ynamic **S**pillover **E**limin**A**tion (REDSEA) as a solution for spillover compensation without loss of signal. This is the step-by-step guidance of how to use REDSEA to produced a FCS file with compensated channel information.

<p align="center"><img width=40%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/overview.png"></p>

### How it works

In brief, the REDSEA algorithm identifies the boundary region for each cell, based on the segmentation mask. Subsequently, for one cell's each channel, signals were subtracted based on the shared boundary with neighboring cells, and their corresponding signal. Moreover, the removed signal from adjacent cells can be reinforced back to the cell. For detailed information look at the Materials and Methods section of the [paper](https://www.biorxiv.org/).

## Table of content

- [Example Data](#example-data)
    - [MIBI data](#mibi-data)
    - [CyCIF data](#cycif-data)

- [Required Inputs](#required-inputs)
    - [MATLAB scripts](#matlab-scripts)
    - [tiff images and channel information](#tiff-images-and-channel-information)
    - [Segmentation Mask](#segmentation-mask)
    
- [REDSEA](#redsea)
    - [Parameters](#parameters)
    - [Sanity plots](#sanity-plots)
    - [Output](#output)

## Example Data

Here we have provided [example datasets](https://github.com/BokaiZhu/REDSEA/tree/master/code) from two different multiplexed imaging modalities, [MIBI](https://www.nature.com/articles/nm.3488) and [CyCIF](https://www.pnas.org/content/110/29/11982).

### MIBI data
Multiplexed Ion Beam Imaging (MIBI) is a method that uses secondary ion mass spectrometry to image antibodies tagged with isotopically pure elemental metal reporters. Here we will use data of non-human primate tissue samples acquired by MIBI (Unpublished).

### CyCIF data
CyCIF is a method for highly multiplexed immuno-fluorescence imaging of formalin-fixed, paraffin-embedded (FFPE) specimens mounted on glass slides. Here we will use human tonsil data acquired with t-CyCIF method. The acquisition of this data is described in this [paper](https://www.nature.com/articles/s41597-019-0332-y) (Rashid, Rumana, et al, Scientific Data 2019)

## Required Inputs


### MATLAB scripts
The [MATLAB scripts](https://github.com/BokaiZhu/REDSEA/tree/master/code) include the REDSEA algorithm, along with other scripts created by [Leeat Keren](https://github.com/lkeren/MIBIAnalysis) for processing multiplexed images, and a fcs processing script ```writeFCS``` created by Jakub Nedbal (https://www.mathworks.com/matlabcentral/fileexchange/42603-writefcs-fname-data-text-other).

### tiff images and channel information
There should be a CSV file containing the channel name information.

<p align="left"><img width=20%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/csv_example.png"></p>

And also a folder containing .tiff images for each channel, where the names should be same as in the CSV file.

### Segmentation Mask
This method requires you to provide a cell segmentation mask. A cell nuclei probability matrix can be produced by your own chose (popular options include [ilastik](https://www.ilastik.org/) or [deepcell](https://github.com/vanvalenlab/deepcell-tf)). In our case for the MIBI data we implemented a trained-in-house deepcell CNN model.

<p align="center"><img width=40%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/feature_1_frame_1_p1_dsDNA.png"></p>

DeepCell is also easy to implement on different imaging modalities. Here is a prediction model we trained with ~ 1500 cells in the [CyCIF](https://www.nature.com/articles/s41597-019-0332-y) dataset:

<p align="center"><img width=100%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/cycIF.png"></p>

In the [example data folder](https://github.com/BokaiZhu/REDSEA/tree/master/code), we have provided a nucleus prediction matrix for the MIBI and CyCIF data, with the name ```feature_1_frame_1_p1_.tif```, in the ```deepCell``` sub folder.

After producing a nuclei probablity mask, we will then use the script ```MibiSegmentByDeepProbWithPerim3.m``` to implement a watershed algorithm for whole cell segmentation. This will produce something like this:

<p align="center"><img width=100%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/ws-segmentation3.png"></p>

<details><summary>See MibiSegmentByDeepProbWithPerim3.m Script</summary>
<p>



```MATLAB
%% Nuclear segmentation using pixel probabilities from deepCell
% Modified from original script from Leeat Keren
% Changed the pipeline to use the probablities from the deepcell instead of 
% the raw nuclear intensities.
% 25July2020, by Yunhao Bai and Sizun Jiang

% Main path for the all the data
mainPath = 'sampleData_MIBI'; %for MIBI
%mainPath = 'sampleData_cycIF'; %for CyCIF
resultsPath = [mainPath,'/segmentResults/']; % for the next step, put the 
% output segmentationParams.mat and PureSegmentation.tif to the
% originalTiff/Point? folder

for p=1:1
for segmentThres=0.01 % change from 0.01 to 0.1, segmentThres defines the 
% local maximum, higher segmentThres leads to fewer cells and more merged
% regions
for probNucThres=0.35 % change from 0.05 to 0.5, probNucThres defines how 
% board the cells will expand
    pointNumber=p;disp(['point',num2str(p)]);
    
    % read .tiff image of nucleus marker to matrix for image  
    pathNucleusMarker = [mainPath,'/originalTIFF/Point',num2str(p),'/dsDNA.tiff']; %for MIBI
    %pathNucleusMarker = [mainPath,'/originalTIFF/Point',num2str(p),'/x7500y3500_1700_DAPI.tif']; %for CyCIF
    t = Tiff(pathNucleusMarker,'r');
    nucIm = read(t);
    %the max value of nucleus channel .tiff
    maxv=25; %for MIBI
    %maxv=3000; %for CyCIF
    rgb_image = MibiGetRGBimageFromMat(nucIm,maxv);

    % %% Get maxima from deepCell probabilities
    % read possibility map from deepCell/ilastik/other segmentation methods
    probNuc = double(imread([mainPath,'/deepCell/feature_1_frame_1_p',num2str(p),'_dsDNA.tif'])); %for MIBI
    %probNuc = double(imread([mainPath,'/deepCell/feature_1_frame_1_p',num2str(p),'_DAPI.tif'])); %for CyCIF
    figure;imagesc(probNuc);

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

    % modify the image so that the background pixels and the extended maxima 
    % pixels are forced to be the only local minima in the image.
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

    cellPerimNewMod = bwperim(newLmod);
    cellPerimNewMod = newLmod;
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
end
end
end

disp('Please put the PureSegmentation.tif and segmentationParams.mat to the corresponding originalTiff folder for next step.');

```

</p>
</details>

This process produced a file called ```segmentationParams.mat``` and stored in a subfolder. This file will be used automatically by the REDSEA algorithm.


## REDSEA

### Parameters
With the ```.tiff``` images, ```.mat``` segmentation file and ```.csv``` channel information, we are now ready to implement REDSEA for boundary compensation.

There are **two methods** for boundary compensation in REDSEA: ```Sudoku``` and ```Cross```. The algorithm walks through the boundaries of each cell, and decides the area to extract signal. You would need to choose one of the two methods and deside how many pixels to expand from the boundary pixel: 

<p align="center"><img width=50%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/method_show.png"></p>

The pixel number for expansion should be **proportional** to the cell size: in our MIBI data, the average cell size is 107 pixels, and we used 2 pixels for expansion; In the [CyCIF tonsil dataset](https://www.synapse.org/#!Synapse:syn17796423) the average cell size is 325 pixels, so 3-4 pixels is recommended. 

Also, you need to supply a list of channel names to perform the compensation process: normally you should only compensate for the **surface markers**, like in our case:
```'CD4';'CD56';'CD21 (CR2)';'CD163';'CD68';'CD3';'CD20';'CD8a'```

Take a look at the annotated code in the block under:

<details><summary>See MibiExtractSingleCellDataFromSegmentationAndTiff_REDSEA.m Script</summary>
<p>



```MATLAB
%% Extract single cell data and do the REDSEA compensation
% Based on the original script from Leeat Keren
% This now reads in any folder of TIFFs with individual channels from the 
% same field of view, and uses that to recreate a countsNoNoise based on 
% the massDS order.
% Then performs REDSEA compensation as implemented by Yunhao Bai
% The outputs will be REDSEA compensated and non-compensated FCS files
% 4May2020, Yunhao Bai, Sizun Jiang


% Main path for the all the data
mainPath = 'sampleData_MIBI'; %for MIBI
%mainPath = 'sampleData_cycIF'; %for CyCIF

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
normChannels = {'CD4';'CD56';'CD21 (CR2)';'CD163';'CD68';'CD3';'CD20';'CD8a'}; %for MIBI
%normChannels = {'x7500y3500_1700_DAPI';'x7500y3500_1700_CD3';'x7500y3500_1700_CD4';'x7500y3500_1700_CD8a';'x7500y3500_1700_CD11b';'x7500y3500_1700_CD20';'x7500y3500_1700_CD45';'x7500y3500_1700_CD68'}; %for CyCIF
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
        t = imread([pathTiff, '/Point', num2str(pointNumber), '/', massDS.Label{i}, '.tiff']); %mind the .tif and .tiff
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
```

</p>
</details>

### Sanity plots

To visually inspect if the parameters described in the previous section are optimal for your data, we provided with a flag in the script ```plotSanityPlots = 1```. Once set to 1, it will produce the pairwise combination scatter plots of all the compensated channels. User can use these plots to evaluate if the compensation is optimal (for example looking at CD4-CD8, CD3-CD20 etc.) It is suggested to run this first with **less** channels, find the optimal parameters, then run the full list.  


Here we showed a sainity plot of CD3-CD20 for the CyCIF data:
<p align="center"><img width=55%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/Sanity_check_x7500y3500_1700_CD3_x7500y3500_1700_CD20.png"></p>
The x and y axis are raw intensity/counts of the channels, scaled by each individual cell's size. In the plot we can see a reduction of CD3-CD20 double positive cells (less orange dots than blue dots in the middle of the plot), and at the same time retention of single-positive signals.


The sanity plots will be produced in the ```sanityPlots``` sub-folder, along with the output of fcs files.

### Output

REDSEA will produce the 4 FCS files for downstream analysis:

```dataFCS.fcs``` = raw counts for each single cell without compensation

```dataRedSeaFCS``` = raw counts for each single cell with REDSEA compensation

```dataScaleSizeFCS.fcs``` = counts for each single cell as scaled by cell size (counts/cellSize) without compensation

```dataRedSeaScaleSizeFCS.fcs``` = counts for each single cell as scaled by cell size (counts/cellSize) with REDSEA compensation

It is recommended to use the ```dataRedSeaScaleSizeFCS.fcs``` file. We can see by using REDSEA compensation, the boundary signal spillover is dynamically eliminated and reinforced.

For the MIBI dataset:
<p align="center"><img width=70%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/mibi_redsea_fcsPlot_new.png"></p>

And REDSEA also works comparably on the CyCIF data:
<p align="center"><img width=40%% src="https://github.com/BokaiZhu/REDSEA/blob/master/media/cycif_redsea_fcs_plot.png"></p>


