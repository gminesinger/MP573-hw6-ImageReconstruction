function g = calcg(x,b,m,lambda,D2,w)

% input x as vector (Nm^2x1), the image guess
% input b as data vector (Nb^2x1)

% f = cost function = ln2 norm squared of MFx-b where F is FFT and M
% samples the FFT of x, b is our data vector

% g = gradient of cost function = (MF)^H * MFx - (MF)^H * b
% g = (F^H * M^H * MFx) - (F^H * M^H * b)
% recall F^H = N^2*ifft and M^H is zero-filling based on mask m

Nm2 = length(m);
Nm = sqrt(Nm2); % m is Nm^2x1

% for now, do nothing with lambda, D2, or w

Fx = fftshift(fft2(ifftshift(reshape(x,[Nm,Nm])))); % take 2D fourier transform of square x
Fx = reshape(Fx,[(Nm*Nm),1]); % reshape back to Nm^2x1 vector
MFx = Fx(logical(m)); % downsample so only using sample according to mask m

fill_MFx(logical(m)) = MFx; % zero-fill MFx based on mask m
fill_b(logical(m)) = b; % zero-fill b based on mask m
fill_MFx = transpose(fill_MFx);
fill_b = transpose(fill_b);

fill_MFx = reshape(fill_MFx,[Nm,Nm]); % reshape into square version, NmxNm
fill_b = reshape(fill_b,[Nm,Nm]); % same for filled b
MFH_MFx = fftshift(ifft2(ifftshift(fill_MFx))); % take fft
MFH_b = fftshift(ifft2(ifftshift(fill_b))); % take fft

MFH_MFx = MFH_MFx*Nm2;
MFH_b = MFH_b*Nm2;

g1 = MFH_MFx - MFH_b;
g1 = reshape(g1,[Nm2,1]);

g2 = lambda*(transpose(D2)*(w.*(w.*(D2*x))));
g = g1+g2;

end

