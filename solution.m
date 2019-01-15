% note: script takes roughly 40 seconds to publish

clear all;
close all;

%% Load Image 1 (Phantom)
n1 = 256; % image 1 (phantom) size is 256x256
X = phantom('Modified Shepp-Logan', n1); % load modified shepp-logan phantom 
rng('default')
rng(5);
% X = imnoise(X, 'gaussian', 0, 0.01); % add noise to image for recovery
% test
x = X(:);

%% Load Image 2 (MATLAB MRI dataset)
load mri;
X2 = D(:,:,:,:); % load mri dataset
X2 = squeeze(X2); % squeeze out redundant dimension
X2 = X2(:,:,16); % extract mid-horizontal slice
X2 = mat2gray(X2); % normalise the 2D image
% test
x2 = X2(:);

%% Parameter specification
n1 = 256; % image 1 (phantom) size is 256x256
n2 = 128; % image 2 (matlab mri) size is 128x128
% linenumber = 11;
% linenumber = 22;
linenumber = 55; % radial line number in transform domain for radial sampling

%% Gradient Plots showing sparsity (Image 1 and Image 2)

% % Image 1
% N = length(x);    
% Dh1 = spdiags([reshape([-ones(round(sqrt(N)),round(sqrt(N))-1)...
%                 zeros(round(sqrt(N)),1)],N,1) ...
%                 reshape([zeros(round(sqrt(N)),1) ones(round(sqrt(N)),...
%                 round(sqrt(N))-1)],N,1)], [0 round(sqrt(N))], N, N);
% Dv1 = spdiags([reshape([-ones(round(sqrt(N))-1,round(sqrt(N))); 
%                 zeros(1,round(sqrt(N)))],N,1), ...
%                 reshape([zeros(1,round(sqrt(N))); ones(round(sqrt(N))-1, ...
%                 round(sqrt(N)))],N,1)], [0 1], N, N);
% Dhy1 = Dh1*x;
% Dvy1 = Dv1*x;
% Dy1 = sqrt(Dvy1.^2 + Dhy1.^2);
% DHY1 = reshape(Dhy1,n1, n1); % horizontal gradient plot
% DVY1 = reshape(Dvy1,n1, n1); % vertical gradient plot
% DY1 = reshape(Dy1,n1, n1); % l2 norm of gradient plot
% imshow(DHY1); 
% figure; imshow(DVY1); 
% figure; imshow(DY1);
% % density of sparse matrices computed for different domains. 
% sparsity_pixel1 = nnz(X)/prod(size(X));
% sparsity_DHY1 = nnz(DHY1)/prod(size(DHY1));
% sparsity_DVY1 = nnz(DVY1)/prod(size(DVY1));
% sparsity_DY1 = nnz(DY1)/prod(size(DY1));
% % histograms of values in different domains to show sparsity
% figure; histogram(x); figure; histogram(Dhy1); figure; histogram(Dvy1);
% figure; histogram(Dy1);
% 
% % Image 2
% N = length(x2);    
% Dh2 = spdiags([reshape([-ones(round(sqrt(N)),round(sqrt(N))-1)...
%                 zeros(round(sqrt(N)),1)],N,1) ...
%                 reshape([zeros(round(sqrt(N)),1) ones(round(sqrt(N)),...
%                 round(sqrt(N))-1)],N,1)], [0 round(sqrt(N))], N, N);
% Dv2 = spdiags([reshape([-ones(round(sqrt(N))-1,round(sqrt(N))); 
%                 zeros(1,round(sqrt(N)))],N,1), ...
%                 reshape([zeros(1,round(sqrt(N))); ones(round(sqrt(N))-1, ...
%                 round(sqrt(N)))],N,1)], [0 1], N, N);
% Dhy2 = Dh2*x2;
% Dvy2 = Dv2*x2;
% Dy2 = sqrt(Dvy2.^2 + Dhy2.^2);
% DHY2 = reshape(Dhy2,n2, n2); % horizontal gradient plot
% DVY2 = reshape(Dvy2,n2, n2); % vertical gradient plot
% DY2 = reshape(Dy2,n2, n2); % l2 norm of gradient
% figure; imshow(DHY2); 
% figure; imshow(DVY2); 
% figure; imshow(DY2);
% % density of sparse matrices computed for different domains. 
% sparsity_pixel2 = nnz(X2)/prod(size(X2));
% sparsity_DHY2 = nnz(DHY2)/prod(size(DHY2));
% sparsity_DVY2 = nnz(DVY2)/prod(size(DVY2));
% sparsity_DY2 = nnz(DY2)/prod(size(DY2));
% % histograms of values in different domains to show sparsity
% figure; histogram(x2); figure; histogram(Dhy2); figure; histogram(Dvy2);
% figure; histogram(Dy2);

%% RECONSTRUCTION EXAMPLE (IMAGE 2)
% Image 2 (MRI MATLAB) will be reconstructed using the log-barrier
% algorithm with TV minimisation under quadratic constraints and 55 radial
% lines. 

%% ORIGINAL IMAGE
imshow(X2)


%% SUBSAMPLING (RANDOM UNIFORM OR RADIAL LINE SUBSAMPING)
% om contains locations of frequencies used in the sampling pattern
% for random uniform sampling, uncomment this
% random_uniform_mu = 0.14; % random uniform subsampling ratios (14%, 28%, 42%)
% random_uniform_mu = 0.28;
random_uniform_mu = 0.42;
% idx = randperm(n1*n1);
% top_index = round(random_uniform_mu*(n1*n1));
% om = idx(1:top_index)';

% radial sampling
om = LineMask(linenumber,n2); % function by Justin Romberg


%% FUNCTION HANDLES
% original image is real so first function returns real and im
% components of 2D fast fourier transform on upper half plane of 
% transform domain. second function does the same for bottom half.
phi = @(zz) phi_fourier_half_plane(zz, om); 
phi_t = @(zz) phi_t_fourier_half_plane(zz, om, n2);


%% BACK PROJECTION (MINIMUM ENERGY/L2 MINIMISATION)

z = phi(x2);
x_backprojection = phi_t(z);
X_backprojection = reshape(x_backprojection, n2, n2);
figure
imshow(X_backprojection)

%% RECONSTRUCTION

% parameters
nu = 2;
eps = 5e-3;
lb_tolerance = 1e-1; newton_tolerance = 1e-1; 
newton_iterlimit = 50; conjgrad_tol = 1e-8;
conjgrad_iterlimit = 200; init_newton_step = 0.9; alpha = 0.01; beta=0.5; 


% recovery of image
tic
%[x_reconstructed, container] = logbar_L1(x_backprojection, phi, phi_t, z, ...
%                              nu, eps, lb_tolerance, newton_tolerance, ...
%                              newton_iterlimit, conjgrad_tol, conjgrad_iterlimit,..._  
%                              init_newton_step, alpha, beta);

[x_reconstructed, container] = logbar_TV(x_backprojection, phi, phi_t, z, ...
                              nu, eps, lb_tolerance, newton_tolerance, ...
                              newton_iterlimit, conjgrad_tol, conjgrad_iterlimit,...  
                              init_newton_step, alpha, beta);
toc
X_RECON = reshape(x_reconstructed, n2, n2);
figure; 
imshow(X_RECON)

%% Duality Gap 
% observe tradeoff between nu and the duality gap

disp('Duality gap plots')

nu1 = 2;
nu2 = 5;
nu3 = 10;

[x_r1, container1] = logbar_TV(x_backprojection, phi, phi_t, z, ...
                              nu1, eps, lb_tolerance, newton_tolerance, ...
                              newton_iterlimit, conjgrad_tol, conjgrad_iterlimit,...  
                              init_newton_step, alpha, beta);

[x_r2, container2] = logbar_TV(x_backprojection, phi, phi_t, z, ...
                              nu2, eps, lb_tolerance, newton_tolerance, ...
                              newton_iterlimit, conjgrad_tol, conjgrad_iterlimit,...  
                              init_newton_step, alpha, beta);
                          
[x_r3, container3] = logbar_TV(x_backprojection, phi, phi_t, z, ...
                               nu3, eps, lb_tolerance, newton_tolerance, ...
                               newton_iterlimit, conjgrad_tol, conjgrad_iterlimit,...  
                               init_newton_step, alpha, beta);

figure;
plot(container1(:,1), container1(:,2))
hold on
plot(container2(:,1), container2(:,2))
plot(container3(:,1), container3(:,2))
xlabel('Cumulative number of Newton iterations')
ylabel('Duality gap')
legend('\nu=2','\nu=5', '\nu=10')
hold off                           

%% NESTA
 
% Psi = @(x2) x2; Psi_trans = Psi;
 
%          
% U = Psi_trans;
% Ut = Psi;
% mu = 0.0005; opts.maxintiter = 5; opts.TOlVar = 1e-5; opts.U = U; % etc etc
%
% 
% tic
% [x_NES,nesta_iters,residuals,errors] = NESTA(phi,phi_t,z,mu,delta,opts);
% toc
% XNES = reshape(X_NES,n2,n2);
% 
% figure; imshow(Xnesta)
%  



