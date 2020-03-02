function [EigenValuesMat,ColorsMat,tVec]=GetEvfd(EvfdParams,Interpulator,ColorVariable)
if nargin==2
    ColorVariable=1;
v2struct(EvfdParams);
tVec=linspace(0,1,NumberOfPointAlongTheGeodesicPath);
EigenValuesMat=[];
% ColorsMat=[];
% f=waitbar(0,sprintf('Calculating eigenvalues flow diagram'));
tmp=fprintf(sprintf('Calculating eigenvalues flow diagram'));
for t_ind=1:length(tVec)
    
    
    K=Interpulator(tVec(t_ind));
    
    [V,D]=eigs(K,NumberOfEigenVals);
    EigenVals=real(diag(D));[~,i]=sort(EigenVals,'descend');
    EigenVals=EigenVals(i)';V=real(V(:,i));
    
    EigenValuesMat=[EigenValuesMat;EigenVals(2:NumberOfEigenVals)];
%     ColorsMat=[ColorsMat;abs(corr(ColorVariable',V(:,2:NumberOfEigenVals)))];
%     fprintf(repmat('\b',numel(tmp)));
%     waitbar(t_ind/length(tVec),f,sprintf('Calculating eigenvalues flow diagram'));
    fprintf(repmat('\b',1,(tmp)));
    tmp=fprintf(sprintf('Calculating eigenvalues flow diagram %g',t_ind/length(tVec)));
end
% close(f)
    fprintf(repmat('\b',1,(tmp)));
EigenValuesMat=EigenValuesMat';
ColorsMat=ones(size(EigenValuesMat));

end

