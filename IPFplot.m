
function [GrainColor]=IPFplot(o,cs,ss,IPFflag)

%%%%%%%% calculate color key and plot color bar
odf = calcODF(o); %% Calculate an ODF
h = [Miller(0,0,1,cs) Miller(1,0,1,cs) Miller(1,1,1,cs)];
cK = ipfHSVKey(cs);
cK.inversePoleFigureDirection=zvector;
GrainColor = cK.orientation2color(o);

% plot grain IPF
figure;
% plotPDF(o,GrainColor,h,'points','all');
plotIPDF(o,GrainColor,vector3d(0,0,1),'points','all','MarkerSize',12,'resolution',0.5*degree);
% LabelContents=arrayfun(@num2str,[1:length(o)],'uniformOutput',false);
% annotate(o,'label',LabelContents{i},'textaboveMarker');
% annotate(rotation('Euler',o(i).phi1,o(i).Phi,o(i).phi2),...
% 'marker','s','MarkerSize',6,'MarkerFaceColor','r',...
% 'label','A','color','w');
% plotIPDF(o,GrainColor,[Miller(0,0,1,cs) Miller(1,0,1,cs) Miller(1,1,1,cs)],'points','all');
% % print([FileName '_IPF001'],'-dtiff','-r300');

% IPF color bar
oM = ipfHSVKey(o);
if IPFflag==1
    figure;
    plot(oM,vector3d(0,0,1),'resolution',0.5*degree,'FontSize',36);
    % print('ColorMap','-dtiff','-r300');
end

