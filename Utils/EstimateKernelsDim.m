function [Dim] = EstimateKernelsDim(K1,K2,debug,TolFac)
if nargin<4
    TolFac=1e-6;
end
[~,DD1]=eig(K1);[~,DD2]=eig(K2);
d1=sort(diag(real(DD1)),'descend'); d2=sort(diag(real(DD2)),'descend');
Dim=min(find(d1<TolFac,1,'first'),find(d2<TolFac,1,'first'));
if isempty(Dim), Dim=-1;,end;
if debug
figure(); plot(d1); hold on; plot(d2);legend('1','2');
end
end

