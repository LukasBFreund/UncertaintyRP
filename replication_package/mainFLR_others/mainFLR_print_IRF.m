%%=========================================================================
% 
% This file prints the IRFs using a standardized setup
% Called by: mainFLR_plot
%
%%=========================================================================


%%  Combine IRF files just loaded
%--------------------------------------------------------------------------

% Put together
if NumModels == 1
aIRFProp = reshape([mIRFProp_1],[size(mIRFProp_1,1),size(mIRFProp_1,2),NumModels]);
elseif NumModels == 2
aIRFProp = reshape([mIRFProp_1,mIRFProp_2],[size(mIRFProp_1,1),size(mIRFProp_1,2),NumModels]);
elseif NumModels == 3
aIRFProp = reshape([mIRFProp_1,mIRFProp_2,mIRFProp_3],[size(mIRFProp_1,1),size(mIRFProp_1,2),NumModels]);
elseif NumModels == 4
aIRFProp = reshape([mIRFProp_1,mIRFProp_2,mIRFProp_3,mIRFProp_4],[size(mIRFProp_1,1),size(mIRFProp_1,2),NumModels]);
end

% Adjust values that are virtually zero (not exactly due to use of a solver) to exactly zero to avoid confusing graphs
aIRFProp(abs(aIRFProp)<1e-10)=0;

%% Plot 
%--------------------------------------------------------------------------

hfig=figure;
for iV = 1:numel(vVariables)
    subplot(NumSubplotV,NumSubplotH,iV)
for iN = 1:NumModels
box on 
hold on
    p1=plot(0:IRFPeriodsNum,100*aIRFProp(1:IRFPeriodsNum+1,strmatch(vVariables{iV},vNames,'exact'),iN),vLinestyle{iN},'LineWidth',linewidthDefault,'Color',vColors{iN});
    hold on
    plot(0:IRFPeriodsNum,zeros(IRFPeriodsNum+1,1),styleZeros,'HandleVisibility','off','Color',colorZeros,'Linewidth',0.5);
    xlim([0 IRFPeriodsNum]);
    set(gca,'XTick',[0:6:IRFPeriodsNum],'FontSize',fontsizeAxisticks,'fontname',fonttype, 'fontweight','bold');
    ax = gca;
    ax.YAxis.Exponent = 0; 
end
    title(vVNames{iV},'FontSize',fontsizeDefault,'fontname',fonttype,'FontWeight',fontweight);
    
    if ismember(vVariables(iV),{'Pii','RFR','R','RReal','u','f','h'}) == 1
        ylabel('Deviation (ppt)','FontSize',fontsizeAxis,'fontname',fonttype, 'fontweight',fontweight);
    elseif ismember(vVariables(iV),{'RP'}) == 1
        ylabel('Deviation (bps)','FontSize',fontsizeAxis,'fontname',fonttype, 'fontweight',fontweight);
    else
       ylabel('Deviation (pct)','FontSize',fontsizeAxis,'fontname',fonttype, 'fontweight',fontweight);
    end    
if iV >= (NumSubplotV*NumSubplotH)-1
    xlabel('Time (months)','FontSize',fontsizeAxis,'fontname',fonttype,'fontweight',fontweight);
end
end

%% Legend 
if OptionLegend == 1
    subplot(NumSubplotV,NumSubplotH,1)
    legend1 = legend(legend_labels);
    set(legend1,'fontname',fonttype,'Location',LegendPosition,'FontSize',fontsizeLegend,'interpreter','latex')
end

%% Print 
%--------------------------------------------------------------------------
if NumSubplotV*NumSubplotH==1
xSize = 17.5/2; 
ySize = 1* 6.25;
xCut = 0;
yCut = 0;
elseif NumSubplotV*NumSubplotH==2
xSize = 17.5; 
ySize = 6.25;
xCut = 1;
yCut = 0;    
elseif NumSubplotV*NumSubplotH==4
xSize = 17.5; 
ySize = 12.5;
xCut = 1;
yCut = 0;  
elseif NumSubplotV*NumSubplotH==6
xSize = 17.5; 
ySize = 3*6.25; 
xCut = 2;
yCut = 0.5;
end

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

