function [medhalftime,stehalftime] = bootstrapwedges(allcells,Rvals,allcellind,allhalftimes,bootstrapcells,whichcells)

if (~exist('bootstrapcells'))
    % bootstrap at the individual cell level
    bootstrapcells = 1;
end

% if (exist('whichcells'))
%     % only use data from particular cells
%     
%     for rc = 1:length(Rvals)
%         keepind = ismember(allcellind{rc},whichcells);
%         allhalftimes{rc} = allhalftimes{rc}(keepind);
%         allcellind{rc} = allh
%     
% end

%% medians and bootstrapping of individual cells
nboot = 100;
stehalftime = zeros(size(Rvals));
clear boothalftimes

for rc = 1:length(Rvals)
    %medhalftimectrl(rc) = median(allhalftimesctrl{rc});
    %medhalftimeRTN4OE(rc) = median(allhalftimesRTN4OE{rc});
    
     %medhalftimectrl(rc) = mean(allhalftimesctrl{rc});
    %medhalftimeRTN4OE(rc) = mean(allhalftimesRTN4OE{rc});
    
    if (rc > length(allcellind))
        medhalftime(rc) = NaN;
        stdhalftime(rc) = NaN;
        continue
    else
        cellind = allcellind{rc};
        % ignore datapoints coming from only one cell
        if (length(unique(cellind))<2)
            medhalftime(rc) = NaN;
            stdhalftime(rc) = NaN;
            continue
        end
    end
    
    halftimes = allhalftimes{rc};
    cellind = allcellind{rc};
    
    %% bootstrapping            
    for bc = 1:nboot    
        
        if (bootstrapcells) % bootstrap whole cells            
            if (exist('whichcells','var'))
                % only use data from specific cells
                keepcells = allcells(whichcells);
            else
                keepcells = allcells;
            end            
            cellsamp = randsample(length(keepcells),length(keepcells),true);      
            
            alltimes = [];
            for cc = 1:length(cellsamp)
                alltimes = [alltimes halftimes(cellind==cellsamp(cc))];
            end
        else % bootstrap individual wedges
            if (exist('whichcells','var'))
                % only use data from specific cells
                keeptimes = halftimes(ismember(cellind,whichcells));
            else
                keeptimes = halftimes;
            end                     
            wedgesamp = randsample(length(keeptimes),length(keeptimes),true);
            alltimes = keeptimes(wedgesamp);
        end
        boothalftimes(bc) = nanmedian(alltimes);
        %boothalftimesctrl(bc) = mean(allhalftimesctrl{rc}(cellsamp));
    end    
    
   medhalftime(rc) = nanmean(boothalftimes);  
   stehalftime(rc) = nanstd(boothalftimes);
end