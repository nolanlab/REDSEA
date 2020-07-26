%% Yunhao Bai for MIBIsubgroup, 20200508
% accessory function for
% MibiExtractSingleCellDataFromSegmentationAndTiff_REDSEA.m to plot scatter
% of every combination of normchannels

function MIBIboundary_compensation_plotting(dataScaleSizeCells,dataCompenScaleSizeCells,normChannels,normChannelsInds,pathSanityPlots)
sz = 7; % size of dots for scatter
perentile = 99; % remove a small percentage of extreme values
for channelSele1 = 1:length(normChannelsInds)
    for channelSele2 = channelSele1+1:length(normChannelsInds)
        f1=figure('Position', [10 10 910 910],'visible','off');
        channelSele1_per = prctile(cat(1,dataScaleSizeCells(:,normChannelsInds(channelSele1)),dataCompenScaleSizeCells(:,normChannelsInds(channelSele1))),perentile);
        channelSele2_per = prctile(cat(1,dataScaleSizeCells(:,normChannelsInds(channelSele2)),dataCompenScaleSizeCells(:,normChannelsInds(channelSele2))),perentile);
        scatter(dataScaleSizeCells(:,normChannelsInds(channelSele1)),dataScaleSizeCells(:,normChannelsInds(channelSele2)),sz,'filled'); hold on;
        scatter(dataCompenScaleSizeCells(:,normChannelsInds(channelSele1)),dataCompenScaleSizeCells(:,normChannelsInds(channelSele2)),sz,'filled');
        legend([char(normChannels(channelSele1)),' to ',char(normChannels(channelSele2)),' pre-compensation, scaled'],...
            [char(normChannels(channelSele1)),' to ',char(normChannels(channelSele2)),' post-compensation, scaled']);
        % in case of extreme situation like no counts in specific channel
        if channelSele1_per > 0
            xlim([0 channelSele1_per]);
        end
        if channelSele2_per > 0
            ylim([0 channelSele2_per]);
        end
        xlabel(char(normChannels(channelSele1)));
        ylabel(char(normChannels(channelSele2)));
        title([char(normChannels(channelSele1)),' to ',char(normChannels(channelSele2)),', scaled by cell area']);
        saveas(f1,[pathSanityPlots,'Sanity_check_',char(normChannels(channelSele1)),'_',char(normChannels(channelSele2)),'.png']);
    end
end
end

