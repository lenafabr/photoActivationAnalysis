function [halftimesbycell, meanbycell, medianbycell] = resortHalfTimesByCell(allhalftimes,allcellind,rvalind)
% resort half times for several different R values by which cell they
% belong to

halftimesbycell = {};
for rc = rvalind 
    for hc = 1:length(allhalftimes{rc})
        cc = allcellind{rc}(hc);
        if (length(halftimesbycell)<cc)
            halftimesbycell{cc} = [];
        end
        halftimesbycell{cc}(end+1) = allhalftimes{rc}(hc);
    end
end

for cc = 1:length(halftimesbycell)
    meanbycell(cc) = mean(halftimesbycell{cc});
    medianbycell(cc) = median(halftimesbycell{cc});
end