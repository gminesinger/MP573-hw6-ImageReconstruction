function [x,cost] = conjugategradient(x,b,m,lambda,D2,w,MAX_ITER)

  Nm2 = length(m);
  Nm = sqrt(Nm2); % m is Nm^2x1

  cost = zeros(MAX_ITER+1,1);
  cost(1) = calcf(x,b,m,lambda,D2,w);

  xnew = x;
  g = calcg(x,b,m,lambda,D2,w);
  gnew = g;
  dnew = -g;

  if norm(g,2)==0
    disp("Your function is minimized, no need to iterate!")
  end
  
  if norm(g,2)~=0
    disp("Your gradient isn't initially zero, must iterate...")
    for k = 1:MAX_ITER
      g = gnew;
      d = dnew;
      x = xnew;
      gH = transpose(g); % matlab automatically calculates the hermitian transpose
      dH = transpose(d);

      % tricky part, calculating "Q" of d
      Qd = fftshift(fft2(ifftshift(reshape(d,[Nm,Nm])))); % reshape to NmxNm and take fft
      Qd = reshape(Qd,[Nm2,1]); % reshape back to vector
      Qd = Qd(logical(m)); % downsamples by m
      Qdfill = zeros(Nm2,1); % create empty array to fill
      Qdfill(logical(m)) = Qd; % fill back in with zeros
      Qdfill = reshape(Qdfill,[Nm,Nm]); % shape into square
      QdfillFFT = Nm*Nm*fftshift(ifft2(ifftshift(Qdfill))); % take IFFT
      Qd1 = reshape(QdfillFFT,[Nm2,1]); % reshape back into 1D array

      Qd2 = lambda*(transpose(D2)*(w.*(w.*(D2*d))));

      Qdfinal = Qd1+Qd2;

      a = -gH*d/(dH*Qdfinal);
      xnew = x + a*d;
      gnew = calcg(xnew,b,m,lambda,D2,w);
      if norm(gnew,2)==0
        break
      else
        gnewH = transpose(gnew);
        B = (gnewH*Qdfinal)/(dH*Qdfinal);
        dnew = -gnew + (B*d);
      end
      cost(k+1)=calcf(xnew,b,m,lambda,D2,w);
    end
  end
end

