function [medhalftime,stehalftime] = bootstrapwedges(allcells,Rvals,allcellind,allhalftimes,options)

% default parameters
opt = struct();
% bootstrap at the individual cell level
opt.bootstrapcells = 1;

% ignore R values that come from less than some min number of cells
opt.mincells = 2;

% which specific cells to use in analysis, if not provided, use all
opt.whichcells = NaN;

if (exist('options','var'))
    opt = copyStruct(options,opt);
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
        if (length(unique(cellind))<opt.mincells)
            medhalftime(rc) = NaN;
            stdhalftime(rc) = NaN;
            continue
        end
    end
    
    halftimes = allhalftimes{rc};
    cellind = allcellind{rc};
    
    %% bootstrapping            
    for bc = 1:nboot    
        if (opt.bootstrapcells) % bootstrap whole cells        
            if (~isnan(opt.whichcells))
                % only use data from specific cells
                keepcells = allcells(opt.whichcells);
            else
                keepcells = allcells;
            end            
            cellsamp = randsample(length(keepcells),length(keepcells),true);      
            
            alltimes = [];
            for cc = 1:length(cellsamp)
                alltimes = [alltimes halftimes(cellind==cellsamp(cc))];
            end
        else % bootstrap individual wedges
            if (~isnan(opt.whichcells))
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