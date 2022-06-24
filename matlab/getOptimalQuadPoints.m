function [w, x, rec, it] = getOptimalQuadPoints(knot, p, w0, x0, knot0, rec)
% intended use: [w, x] = getOptimalQuadPoints(knot,p)
	warn_state  = warning();
	warning off;
  n  = numel(knot)-p-1; % dimension of our spline space

  if mod(n,2)==1 % need to have a space of even dimension
    i = find(diff(knot) == max(diff(knot))); 
    i = i(ceil(end/2));                    % insert new knot in middle of
    knot = sort([knot,mean(knot(i:i+1))]); % the largest, centermost knot span
    n = n+1;
  end

  % compute all greville points and integrals (used for initial guess)
  greville       = zeros(n,1);
  exact_integral = zeros(n,1);
  for i=1:n
    greville(i)       = sum(knot(i+1:i+p)) / p;
    exact_integral(i) = (knot(i+p+1)-knot(i))/(p+1);
  end

  if exist('x0') % if initial guess is provided, use these
    x = x0;
    w = w0;
  else           % else compute them based on greville points and integrals
    w   = (exact_integral(1:2:end) + exact_integral(2:2:end))  ;
    x   = (      greville(1:2:end) +       greville(2:2:end))/2;
    rec = 1;     % counter variable to count the number of recursive calls
  end

  newton_tol       = 1e-10;    % convergence tolerance
  newton_max_it    = 15;       % max iterations before divergence
  while true                   % recursive loop from algorithm 2
    for it = 1:newton_max_it   % newton iteration loop
      [N dN] = BSpline(knot, p, x);
      F      = N*w - exact_integral;
      dF     = [N, dN*diag(sparse(w))];

      dx = dF \ -F;
      w = w + dx(1:end/2);
      x = x + dx(end/2+1:end);

      % test for diverging (coarse heuristic, see section 3.3)
      if( min(x)<knot(1)  )     break;   end;
      if( max(x)>knot(end))     break;   end;

      % test for converging 
      if(norm(dx)<newton_tol)  warning(warn_state); return;  end;
    end 

    % at this point, newton iteration has diverged. solve recursively on easier knot
    if exist('knot0')
      [w, x, rec] = getOptimalQuadPoints((knot0 + knot)/2, p, w0, x0, knot0, rec);
      knot0       = (knot0 + knot)/2;
    else
      uniformKnot = linspace(knot(1),knot(end), n+p+1);
      [w, x, rec] = getOptimalQuadPoints(uniformKnot, p);
      knot0       = uniformKnot;
    end
    rec = rec + 1;
    x0 = x;
    w0 = w;
  end % loop up and start newton iteration with better initial guess
end
