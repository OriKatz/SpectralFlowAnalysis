function [A, r, c] = sinkhornKnopp(A, varargin)
%SINKHORNKNOPP normalises a matrix to be doubly stochastic
%   M = SINKHORNKNOPP(A) takes a nonnegative NxN matrix A and normalises it
%   so that the sum of each row and the sum of each column is unity. M is
%   equal to DAE where D and E are diagonal matrices with positive values
%   on the diagonal.
% 
%   M = SINKHORNKNOPP(A, NAME1, VALUE1, ...) allows other parameters to be
%   set. These are:
% 
%       'Tolerance' - a positive scalar, default EPS(N). The value is the
%       maximum error in the column sums of M; the row sums will be correct
%       to within rounding error.
% 
%       'MaxIter' - a positive integer, default Inf. The maximum number of
%       iterations to carry out.
% 
%   [M, R, C] = SINKHORNKNOPP(...) also returns normalising vectors R and C
%   such that M = diag(R) * A * diag(C) to within a small tolerance.
% 
%   Example
%   -------
% 
%       a = toeplitz(1:6);
%       m = sinkhornKnopp(a);
%       disp('a'); disp(a);
%       disp('m'); disp(m);
%       disp('Row and column sums'); disp(sum(m,1)); disp(sum(m,2));
% 
%   Convergence
%   -----------
% 
%   The algorithm will converge for positive matrices, but may not converge
%   if there are too many zeros in A, depending on their distribution. In
%   such cases it may be necessary to set 'MaxIter' and to check the column
%   sums of M. For some applications adding a small constant to A is
%   recommended.
%
%   Algorithm
%   ---------
% 
%   The Sinkhorn-Knopp algorithm, also known as the RAS method and
%   Bregman's balancing method, is used. The code is modified from Knight
%   (2008), avoiding the matrix transpose.
%
%   Reference
%   ---------
%
%   Philip A. Knight (2008) The Sinkhorn?Knopp Algorithm: Convergence and
%   Applications. SIAM Journal on Matrix Analysis and Applications 30(1),
%   261-275. doi: 10.1137/060659624

% Copyright David Young 2015

% Input parameter parsing and checking
validateattributes(A, {'numeric'}, {'nonnegative' 'square'});
N = size(A, 1);

inp = inputParser;
inp.addParameter('Tolerance', eps(N), ...
    @(x) validateattributes(x, {'numeric'}, {'positive' 'scalar'}));
inp.addParameter('MaxIter', Inf, @(x) ...
    checkattributes(x, {'numeric'}, {'positive' 'integer' 'scalar'}) ...
    || (isinf(x) && isscalar(x) && x > 0));
inp.parse(varargin{:});
tol = inp.Results.Tolerance;
maxiter = inp.Results.MaxIter;
maxiter = 1e4;
% first iteration - no test
iter = 1;
c = 1./sum(A);
r = 1./(A * c.');

% subsequent iterations include test
while iter < maxiter
    iter = iter + 1;
%     iter
    cinv = r.' * A;
    % test whether the tolerance was achieved on the last iteration
    if  max(abs(cinv .* c - 1)) <= tol
        break
    end
    c = 1./cinv;
    r = 1./(A * c.');
end

A = A .* (r * c);

end

%%
% A=abs(randn(size(A)));A=A+A';
% 
% N1=diag(randn(1,size(A,1)));
% N2=diag(randn(1,size(A,1)));
% 
% K=N1*A*N1;
% 
% K=diag(r)*A*diag(c);
% 
% 
% [V1,D1]=eig(A);[c,i]=sort(diag(D1),'descend');V1=V1(:,i);D1=D1(i,i);
% [V2,D2]=eig(K);D2=real(D2);[c,i]=sort(real(diag(D2)),'descend');V2=V2(:,i);D2=D2(i,i);
% 
% figure(); plot(diag(D1));
% hold on;plot(diag(D2))
% 
% figure(); myplot(V2(:,1:5)); hold on;myplot(V1(:,1:5));
% 
% 
% tmp=A*V2(:,2);
% figure(); plot(tmp/norm(tmp)); hold on; plot(V2(:,2),'*')
% 
% DegMat=diag(sum(A))
% L=DegMat-A;
% [VL,DL]=eig(L);[c,i]=sort(diag(DL),'ascend');VL=VL(:,i);DL=DL(i,i);
% 
% tmp=L*V2(:,10);
% figure(); plot(tmp/norm(tmp)); hold on; plot(V2(:,10),'*')
% 
% DegMat=diag(sum(A));
% NL=eye(size(A))-DegMat^(-1/2)*A*DegMat^(-1/2);
% [VNL,DNL]=eig(NL);[c,i]=sort(diag(DNL),'ascend');VNL=VNL(:,i);DNL=DNL(i,i);
% 
% tmp=K*V2(:,5);
% figure(); plot(tmp/norm(tmp)); hold on; plot(V2(:,5),'*')
% % 
% 
% 
% 
% tmp=L*VNL(:,5);
% figure(); plot(tmp/norm(tmp)); hold on; plot(VNL(:,5),'*')




% figure(); myplot(VNL(:,1:5))
