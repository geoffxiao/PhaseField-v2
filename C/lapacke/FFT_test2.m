A = zeros(3,3,3);
c = 0;
for i = 1 : 3
for j = 1 : 3
for k = 1 : 3
    A(i,j,k) = c;
    c = c + 1;
end
end
end

Ak = fftn(A);

for i = 1 : 27
   fprintf('%d + %di\n', real(Ak(i)), imag(Ak(i)));
end