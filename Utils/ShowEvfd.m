function  ShowEvfd(EigenValuesMat,ColorsMat,tVec,FigsFolder)

% Plot eigenvalues flow diagram
GeodesicEigenValuesRawStack=log(abs(EigenValuesMat(:)));
GeodesicCorrEigenValuesRawStack=ColorsMat(:);
tRawStack=kron(tVec(1:size(EigenValuesMat,2)),ones(1,size(EigenValuesMat,1)))';
EigenValuesTupple=[GeodesicEigenValuesRawStack,tRawStack];

% figure();
scatter((EigenValuesTupple(:,1)),EigenValuesTupple(:,2),50,(GeodesicCorrEigenValuesRawStack),'filled','o');
colormap jet;
% caxis([0 1]);
% colorbar;
% ylabel('$t$','FontSize',25);xlabel('$-\log(\mu^t_i)$','FontSize',25);title('Eigenvalues flow diagram','FontSize',25);
xlim([min(EigenValuesTupple(:,1)),max(EigenValuesTupple(:,1))]);
if nargin>3
    set(gcf, 'Position', get(0, 'Screensize'));
    MySaveFig(gcf,FigsFolder,"EigenvaluesFlowDiagram");close(gcf);
end

end

