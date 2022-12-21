%% Load data
clear
load hw3_problem1.mat % Pick right mat file to load for problem 1-2-3
data = b; % vector containing the measured data
% Note that m contains a binary mask describing the Fourier-space data
% points included in this problem (ie: the number of True entries in m
% is the length of b)
lambda = 0; % lambda=0 is problem 1, lambda>0 is problems 2 and 3
INIT_ZEROFILLED = 1; % Whether to initialize with the zero-filled ifft solution 
niter = 100; % Number of overall iterations
niterGS = 100; % Number of golden section search iterations for SD
s1 = sqrt(length(m));s2 = s1; % s1 and s2 are the dimensions of the image 
if ~exist('w') % If we don't have weights w defined, let's assume all 1's
    w = ones(s1*s2,1);
end
m2D = reshape(m,[s1,s2]);% If you wish to view the mask of sampled Fourier space as a 2D array

% Create the 2D finite difference taking matrix D2
D = 2*eye(s1) - circshift(eye(s1),[0, -1]) - circshift(eye(s1),[0, 1]);
D = sparse(D); 
I = speye(s1);
D2 = kron(I,D) + kron(D,I);

%% Zero-filled solution (also LS solution)
fx1 = zeros([s1,s2]);
fx1(m) = data;
x1 = fftshift(ifft2(ifftshift(fx1)));
figure(1);imagesc(abs(x1));axis equal tight off;title('Zero-filled solution','FontSize',18);

%% For each problem (1-2-3), need to implement SD and CG with several regularization parameters 
%% as described in the homework set
