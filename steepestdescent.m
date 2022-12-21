function [x,cost] = steepestdescent(x,b,m,lambda,D2,w,MAX_ITER,MAX_ITER_GS)

disp("Algorithm: Steepest Descent")

cost = zeros(MAX_ITER+1,1);
cost(1) = calcf(x,b,m,lambda,D2,w);

p = 0.382; % golden section

for k = 1:MAX_ITER
    
    g = calcg(x,b,m,lambda,D2,w); %,0,0,0)

    if norm(g,2)>0.000001 % if gradient not tiny

      d = -g/max(abs(g)); % direction of travel

      % Initializing golden search
      a0 = 0;
      b0 = 1;
      % performing golden search for alpha
      for k_GS = 1:MAX_ITER_GS
        a1 = x + d*a0 + d*(p*(b0-a0));
        b1 = x + d*a0 + d*((1-p)*(b0-a0));
        fa1 = calcf(a1,b,m,lambda,D2,w); %,0,0,0);
        fb1 = calcf(b1,b,m,lambda,D2,w); %,0,0,0);
        if fb1 < fa1
          a0 = a0 + p*(b0-a0); %# made range smaller, keep lower side
        else
          b0 = a0 + (1-p)*(b0-a0); %# made range smaller, keep upper side
        end
        %print("k_GS = ",k_GS,"alpha = ",(a0+b0)/2)
        %disp("alpha iteration GS: ")
        %disp(k_GS)
      end
      alpha = (a0+b0)/2;
 
      %# If we actually made the cost smaller, let's apply the step
      xtentative = x + alpha*d; %# Let us check that we actually descend
      ftentative = calcf(xtentative,b,m,lambda,D2,w); %,0,0,0)
      fold = calcf(x,b,m,lambda,D2,w); %,0,0,0)
      if ftentative < fold % if there's any improvement, set x to tentative x
          x = xtentative;
      else % if there's no improvement, keep old x
          %disp("No descent, may need more GS iterations/smaller step")
      end
    
    end % end of "if gradient not tiny"
  
  fnew = calcf(x,b,m,lambda,D2,w); % ,0,0,0); % calculate updated 
  cost(k+1) = fnew;
  disp("SD iteration #: ")
  disp(k)

end

