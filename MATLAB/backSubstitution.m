function x=backSubstitution(U,b,n)
% Solving an upper triangular system by back-substitution
% Input matrix U is an n by n upper triangular matrix
% Input vector b is n by 1
% Input scalar n specifies the dimensions of the arrays
% Output vector x is the solution to the linear system
% U x = b
% K. Ming Leung, 01/26/03

x=zeros(n,1);
for j=n:-1:1
    if (U(j,j)==0) error('Matrix is singular!'); end;
    x(j)=b(j)/U(j,j);
    b(1:j-1)=b(1:j-1)-U(1:j-1,j)*x(j);
end