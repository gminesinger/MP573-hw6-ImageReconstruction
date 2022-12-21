% let's start by setting up our problem
close all force

problem = "B";
if problem == "A"
    load("hw6_recon_problem1.mat")
    Nm2 = length(m);
    Nm = sqrt(Nm2); % m is Nm^2x1
    lambda = 0;
    D2 = 0;
    w = 0;
elseif problem == "B"
    load("hw6_recon_problem2.mat")
    Nm2 = length(m);
    Nm = sqrt(Nm2); % m is Nm^2x1
    lambda = 10^8;
    D = 2*eye(Nm) - circshift(eye(Nm),[0, -1]) - circshift(eye(Nm),[0, 1]);
    D = sparse(D); 
    I = speye(Nm);
    D2 = kron(I,D) + kron(D,I);
    w = ones([Nm2,1]);
else
    load("hw6_recon_problem3.mat")
    Nm2 = length(m);
    Nm = sqrt(Nm2); % m is Nm^2x1
    lambda = 10^8;
    D = 2*eye(Nm) - circshift(eye(Nm),[0, -1]) - circshift(eye(Nm),[0, 1]);
    D = sparse(D); 
    I = speye(Nm);
    D2 = kron(I,D) + kron(D,I);
    % w is loaded in with file
end

% initial guess of zero, called x0
x0 = zeros(Nm2,1);

% initial guess of zero-filled solution, called x1
fx1 = zeros([Nm,Nm]);
fx1(m) = b;
x1 = fftshift(ifft2(ifftshift(fx1)));
%figure(1);imagesc(abs(x1));axis equal tight off;title('Zero-filled solution','FontSize',18);
x1 = reshape(x1,[Nm2,1]);

calcf(x1,b,m,lambda,D2,w)
g=calcg(x1,b,m,lambda,D2,w);
disp(g(1))

% define your max iterations
MAX_ITER = 100;
MAX_ITER_GS = 100;

% % run steepest descent algorithm for either initial x choice
[xnew0_SD,cost0_SD] = steepestdescent(x0,b,m,lambda,D2,w,MAX_ITER,MAX_ITER_GS);
[xnew1_SD,cost1_SD] = steepestdescent(x1,b,m,lambda,D2,w,MAX_ITER,MAX_ITER_GS);

% run conjugate gradient for either initial x choice
[xnew0_CG,cost0_CG] = conjugategradient(x0,b,m,lambda,D2,w,MAX_ITER);
[xnew1_CG,cost1_CG] = conjugategradient(x1,b,m,lambda,D2,w,MAX_ITER);

figure;
plot(log10(1:MAX_ITER+1),log10(cost0_SD))
xlabel('log(Iteration)'); 
ylabel('log(Cost function)'); 

figure;
subplot(1,2,1)
imagesc(abs(reshape(x0,[Nm,Nm])));axis equal tight off;title('Zeroed guess','FontSize',18);
subplot(1,2,2)
imagesc(abs(reshape(xnew0_SD,[Nm,Nm])));axis equal tight off;title('Steepest Descent Solution','FontSize',18);

figure;
scatter(1:MAX_ITER+1,log(cost1_SD))
xlabel('Iteration'); 
ylabel('Cost function'); 

figure;
subplot(1,2,1)
imagesc(abs(reshape(x1,[Nm,Nm])));axis equal tight off;title('Zero-filled guess','FontSize',18);
subplot(1,2,2)
imagesc(abs(reshape(xnew1_SD,[Nm,Nm])));axis equal tight off;title('Steepest Descent Solution','FontSize',18);


figure;
plot(log10(1:MAX_ITER+1),log10(cost0_CG))
xlabel('log(Iteration)'); 
ylabel('log(Cost function)'); 

figure;
subplot(1,2,1)
imagesc(abs(reshape(x0,[Nm,Nm])));axis equal tight off;title('Zeroed guess','FontSize',18);
subplot(1,2,2)
imagesc(abs(reshape(xnew0_CG,[Nm,Nm])));axis equal tight off;title('Conjugate Gradient Solution','FontSize',18);


figure;
scatter(1:MAX_ITER+1,log(cost1_CG))
xlabel('Iteration'); 
ylabel('Cost function'); 

figure;
subplot(1,2,1)
imagesc(abs(reshape(x1,[Nm,Nm])));axis equal tight off;title('Zero-filled guess','FontSize',18);
subplot(1,2,2)
imagesc(abs(reshape(xnew1_CG,[Nm,Nm])));axis equal tight off;title('Conjugate Gradient Solution','FontSize',18);

