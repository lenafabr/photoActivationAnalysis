% ---------------
%% load in data
% This load cell objects that were created with processManyCells.m
% ---------------

load('../celldata/allcells_ctrl.mat')
load('../celldata/allcells_RTN4OE.mat')

%%  calculate halftimes for all regions

options = struct('dodisplay',1,'maxtscl',5);
[avghalftime_ctrl,allhalftimes_ctrl] = getAvgHalfTimes(allcellsctrl,options);
%%
[avghalftime_RTN4OE,allhalftimes_RTN4OE] = getAvgHalfTimes(allcellsRTN4OE,options);

%% avg size of activated region
actsizectrl = []; actsizeRTN4OE = [];
for cc = 1:length(allcellsctrl)
    CL = allcellsctrl(cc);
    actsizectrl(end+1) = CL.actROI.rad/CL.pxperum;
end
for cc = 1:length(allcellsRTN4OE)
    CL = allcellsRTN4OE(cc);
    actsizeRTN4OE(end+1) = CL.actROI.rad/CL.pxperum;
end

%% plot half-time versus radius
figure(1)
plot(Rvals,avghalftime_ctrl,'b.-','LineWidth',2,'MarkerSize',15)
hold all
plot(Rvals,avghalftime_RTN4OE,'r.-','LineWidth',2,'MarkerSize',15)
hold off

xlabel('radius (um)')
ylabel('avg half-time (s)')

legend('control','RTN4OE')

%% plot scaling on log-log axes
figure(2)
% shift R axis by typical size of activated region
loglog(Rvals-2,avghalftime_ctrl,'b.-','LineWidth',2,'MarkerSize',15)
hold all
loglog(Rvals-2,avghalftime_RTN4OE,'r.-','LineWidth',2,'MarkerSize',15)

loglog(Rvals-2,0.5*(Rvals-2).^2,'k--','LineWidth',2)
loglog(Rvals-2,2*(Rvals-2),'m--','LineWidth',2)
hold off

xlabel('distance from center - 2 um')
ylabel('avg half-time (s)')

legend('control','RTN4OE', '~R^2','~R')

%% view some of the exponential fits
figure(3)
CL = allcellsctrl(2);
wedgeind = find(strcmp({CL.ROIs.type},'wedge'));
ringind = find(strcmp({CL.ROIs.type},'ring'));
tvals= (1:CL.NFrame)*CL.dt;
cmat = jet(length(wedgeind));

subplot(1,2,1)
imshow(CL.ERimg,[])
drawrectangle('Position',CL.actROI.rect,'Color','magenta')
hold all
for rc = 1:10:length(wedgeind)
    bound = CL.ROIs(wedgeind(rc)).bound;
    plot(bound(:,1),bound(:,2),'LineWidth',2,'Color',cmat(rc,:))
end
hold off
title(sprintf('%s',CL.Name),'Interpreter','none')

subplot(1,2,2)
for rc = 1:10:length(wedgeind)
    plot(tvals,CL.ROIs(wedgeind(rc)).avgsignal,'Color',cmat(rc,:))
    hold all
    tfit = tvals(CL.startPA+1:end);
    plot(tfit,CL.fitfunc(CL.ROIs(wedgeind(rc)).cfit,tfit),'Color',cmat(rc,:),'LineWidth',2)
end
hold off
xlabel('time (sec)')
ylabel('signal')