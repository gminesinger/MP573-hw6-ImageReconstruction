function f = calcf(x,b,m,lambda,D2,w)

% input x as vector (Nm^2x1), the image guess
% input b as data vector (Nb^2x1)

% f = cost function = ln2 norm squared of MFx-b where F is FFT and M
% samples the FFT of x, b is our data vectore

Nm = sqrt(length(m)); % m is Nm^2x1

Fx = fftshift(fft2(ifftshift(reshape(x,[Nm,Nm])))); % take 2D fourier transform of square x
Fx = reshape(Fx,[(Nm*Nm),1]); % reshape back to Nm^2x1 vector
MFx = Fx(logical(m)); % downsample so only using sample according to m mask

f1 = 0.5*(norm((MFx-b),2))^2; % this output is a scalar

f2 = 0.5*lambda*(norm((w.*(D2*x)),2))^2; % relevant for parts B and C

f = f1+f2;

end

