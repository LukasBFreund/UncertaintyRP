%% Decomposition
%--------------------------------------------------------------------------

NumModels = 3;

load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_baseline')); mIRF0_1=mIRFProp_zUncertainty_EMAS; 
uPos = strmatch('u',vNames,'exact');
PiiPos = strmatch('Pii',vNames,'exact');

mIRF0_1(:,uPos,:) = mIRF0_1(:,uPos,:)*vEMAS(uPos);
mIRF0_1(:,PiiPos,:) = mIRF0_1(:,PiiPos,:)*vEMAS(PiiPos); % don't annualize for cumulation

load(fullfile('.', 'Output/IRFs/', 'IRFs_FLR_noRP')); mIRF0_2=mIRFProp_zUncertainty_EMAS; 
mIRF0_2(:,uPos,:) = mIRF0_2(:,uPos,:)*vEMAS(uPos);
mIRF0_2(:,PiiPos,:) = mIRF0_2(:,PiiPos,:)*vEMAS(PiiPos);

% Simple decomposition by turning off covariance
mIRF_1 = mIRF0_1;
mIRF_2 = mIRF0_1 - mIRF0_2;     % Risk premium
mIRF_3 = mIRF0_2;               % Demand/other 

% Design
OptionLegend = 1;
LegendPosition = 'southeast';
OptionPlotExpectations = 0; 
vVariables = {'u','Pii'};
vVNames = {'Unemployment rate','Inflation rate (ann.)'};
NumSubplotV = 1;
NumSubplotH = 1;

Labels{1} = 'Total';
Labels{2} = 'Risk premium';
Labels{3} = 'Demand';

%% Plot 
%--------------------------------------------------------------------------

uPos = strmatch('u',vNames,'exact');
PiiPos = strmatch('Pii',vNames,'exact');
vVNames = {'Unemployment rate','Inflation rate'};

x = categorical(vVNames);
x = reordercats(x,vVNames);

y = 100*[
    sum(mIRF_1(:,uPos)) sum(mIRF_2(:,uPos)) sum(mIRF_3(:,uPos)) ;
    sum(mIRF_1(:,PiiPos)) sum(mIRF_2(:,PiiPos)) sum(mIRF_3(:,PiiPos)) ;
];

% contribution of RP
%y(1,2)/y(1,1)

a = [y(1,:)' nan(NumModels,1)]; % unemployment
b = [nan(NumModels,1) y(2,:)']; % inflation
myC= [vColors{1};vColors{2};vColors{3};vColors{4};[120,12,31]/255];  

hfig=figure(2);
yyaxis left
H=bar(x,a');
ylabel('Deviation (ppt)','Color','k','FontSize',fontsizeAxis,'fontname',fonttype, 'fontweight',fontweight);

for k=1:NumModels
  set(H(k),'facecolor',myC(k,:))
  set(H(k),'linestyle',vLinestyle{k})
end

boundaryY = round(max(y,[],'all')*1.2,1);
hold on
yyaxis right
set(gca,'Ylim',[-boundaryY +boundaryY])
set(gca,'ytick',[-boundaryY:0.1:boundaryY],'FontSize',fontsizeAxisticks,'fontname',fonttype,'fontweight',fontweight);
H=bar(x,b');
ylabel('Deviation (ppt)','Color','k','FontSize',fontsizeAxis,'fontname',fonttype,'fontweight',fontweight);
for k=1:NumModels
  set(H(k),'facecolor',myC(k,:))
  set(H(k),'linestyle',vLinestyle{k})
end
legend1 = legend(Labels(1:NumModels));
set(legend1,'fontname',fonttype,'Location','northeast','FontSize',fontsizeLegend,'fontweight',fontweight);%,'interpreter','latex')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%align zero for left and right
yyaxis left; 
set(gca,'Ylim',[-boundaryY +boundaryY])
set(gca,'ytick',[-boundaryY:0.2:boundaryY],'FontSize',fontsizeAxisticks,'fontname',fonttype);

ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
yyaxis right; yliml = get(gca,'Ylim');
if yliml(2)*ratio<yliml(1)
    set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
else
    set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
end
set(gca,'ytick',[-boundaryY:0.2:boundaryY],'FontSize',fontsizeAxisticks,'fontname',fonttype);

set(gca,'fontname',fonttype);


%% Print 
%--------------------------------------------------------------------------

xSize = 10; 
ySize = 6.25;
xCut = 0;
yCut = 0;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-xCut ySize-yCut],'PaperPositionMode','auto')
  
if OptionPrint == 1      
    if optionGreycolor == 0
         FigureNamepdf =horzcat(horzcat(TargetPath,horzcat(name,'_CS')),'.pdf');       
    else
         FigureNamepdf =horzcat(horzcat(TargetPath,name),'.pdf');
    end
        exportgraphics(hfig,FigureNamepdf,'ContentType', 'vector'); 
end
