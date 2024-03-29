function plotTraces(ROIs,PlotAxes,textlabels)

% plot some (normalized) traces from ROIs, with Gaussian peaks identified
% PlotAxes = which axes to use
% textlabels = include number labels for the traces? default is off

if (~exist('textlabels','var'))
    textlabels = false;
end

regionTraces = vertcat(ROIs.avgsignal);

% normalize traces
minvals = min(regionTraces')';
tracerange = range(regionTraces,2);
regionTracesNorm = (regionTraces-minvals)./tracerange;

% shift up traces to view many at once
regionTraceshift = regionTracesNorm + (1:size(regionTraces,1))';
tlist = 1:size(regionTraces,2);

xvals0= linspace(-5,5,100);
gauss0 = normpdf(xvals0)*sqrt(2*pi);

roicmap = jet(length(ROIs));
textlabels = true;

for rc = 1:length(ROIs)
    plot(PlotAxes, tlist,regionTraceshift(rc,:),'Color',roicmap(rc,:))
    hold(PlotAxes,'all')
    if (isfield(ROIs(1),'puffind'))
        puffind = ROIs(rc).puffind;
        plot(PlotAxes, tlist(puffind),regionTraceshift(rc,puffind),'o','Color',roicmap(rc,:))
    end
    if (isfield(ROIs(rc),'puffwidths'))
        % plot Gaussian signal overlay
        puffind = ROIs(rc).puffind;
        signal = ROIs(rc).avgsignal;
        tlist = 1:length(ROIs(rc).avgsignal);
        basesignal = ROIs(rc).basesignal;

        for cc= 1:length(puffind)
            base = ROIs(rc).basesignal(puffind(cc));
            w = ROIs(rc).puffwidths(cc);
            gauss = gauss0*(signal(puffind(cc))-base) + base;
            xvals = xvals0*w/2 + tlist(puffind(cc));

            gaussshift  = (gauss-minvals(rc))/tracerange(rc)+rc;
            baseshift = (basesignal-minvals(rc))/tracerange(rc)+rc;
            
            plot(PlotAxes,tlist,baseshift,'Color',roicmap(rc,:))
            plot(PlotAxes,xvals,gaussshift,'--','Color',roicmap(rc,:))
        end
    end
    if (textlabels)
        text(PlotAxes,tlist(end),regionTraceshift(rc,end),sprintf('%d',rc),'Color',roicmap(rc,:))
    end
end
hold(PlotAxes,'off')
end