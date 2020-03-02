function [D1,D2,A1,A2,K1,K2,Scale1,Scale2]  = GetKernels(S1,S2,KernelsParams)
v2struct(KernelsParams);

D1= pdist2(S1,S1,'euclidean');
D2= pdist2(S2,S2,'euclidean');

if isscalar(NormFac)
    Scale1=median(D1(:))*NormFac;
    Scale2=median(D2(:))*NormFac;
    switch UseMutualScaleFac
        case 'Min'
            Scale1=min([Scale1,Scale2]);Scale2=Scale1;
        case 'Mean'
            Scale1=(Scale1+Scale2)/2;Scale2=Scale1;
            % eps1=sqrt(eps1*eps2);
            % eps1=sqrt(eps1^2+eps2^2)'
    end
else
    Scale1=median(D1(:))*NormFac(1);
    Scale2=median(D2(:))*NormFac(2);
end



if any(NormFac>0)
    A1=exp(-D1.^2/(Scale1^2));%tmp=sum(A1,2); K1=diag(1./tmp)*A1;
    K1 = sinkhornKnopp(A1);
    A2=exp(-D2.^2/(Scale2^2));%tmp=sum(A2,2); K2=diag(1./tmp)*A2;
    K2 = sinkhornKnopp(A2);
else
    if isscalar(NormFac)
        [A1,~]=AffinityFromDistance( D1, abs(NormFac));% [A1] = SuppressDiagonal(A1);
        K1=sinkhornKnopp(A1);K1=K1^1;
        
        [A2,~]=AffinityFromDistance( D2,abs(NormFac));% [A2] = SuppressDiagonal(A2);
        K2=sinkhornKnopp(A2);K2=K2^1;
    else
        [A1,~]=AffinityFromDistance( D1, abs(NormFac(1)));% [A1] = SuppressDiagonal(A1);
        K1=sinkhornKnopp(A1);K1=K1^1;
        
        [A2,~]=AffinityFromDistance( D2,abs(NormFac(2)));% [A2] = SuppressDiagonal(A2);
        K2=sinkhornKnopp(A2);K2=K2^1;
    end
end


