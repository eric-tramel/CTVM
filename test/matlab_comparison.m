L = 11;
N = L*L;
M = N;

A = randn(M,N)./sqrt(N);

X = double(imread('data/box.png'));
X = imresize(X,[L L]);

y = A*X(:);

opts.mu = 2^8;			% Should be small for high noise (i.e. here)
opts.beta = 2^5;
opts.tol = 1E-3;
opts.maxit = 32;
opts.TVnorm = 2;
opts.nonneg = false;
opts.disp = true;

t = cputime;
[U, out] = TVAL3(A,y,L,L,opts);
t = cputime - t;


figure(1);
    imagesc(U);
    axis image;
    colormap(gray);



