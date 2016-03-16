function [N dN] = BSpline(knot, p, t)
% function [N dN] = BSpline(knot, p, t)
%     parameters:
%         knot - the knot vector 
%         p    - the polynomial order of the basis
%         t    - m component vector of points which is to be evaluated
%     returns:
%         N    - n by m matrix of the solution of all basis functions i 
%                evaluated at all points xi(j), given in N(i,j)
%         dN   - n by m matrix of the solution of all derivative of the 
%                basis functions i evaluated at all points xi(j), given 
%                in dN(i,j)

% pad knot vector, so we always compute with C^{-1} at the start/end
knot = [knot(1)*ones(1,p), knot, knot(end)*ones(1,p)];
n   = numel(knot)-p-1;          % number of basis functions +2p
Ni  = ones( numel(t)*(p+1),1);
Nj  = ones( numel(t)*(p+1),1);
Nv  = zeros(numel(t)*(p+1),1);
dNv = zeros(numel(t)*(p+1),1);
for i=1:numel(t),
  if t(i)==knot(end-p)     % evaluate right end-point from the left
    mu = find(knot>=t(i), 1);
  else                     % else evlauate in the limit from the right
    mu = find(knot>t(i), 1);
  end
  if numel(mu)==0 || mu==1 % evaluation outside domain
    continue;
  end
  mu = mu-1;               % index of last non-zero basis function
  
  N  = 1;
  for q=1:p,
    k = mu-q+1:mu;
    R = zeros(q+1,q);
    R(1:q+2:end) = (knot(k+q) - t(i)   ) ./ (knot(k+q)-knot(k));
    R(2:q+2:end) = (t(i)      - knot(k)) ./ (knot(k+q)-knot(k));
    if p==q
      dR = zeros(q+1,q);
      dR(1:q+2:end) = -p ./ (knot(k+q)-knot(k));
      dR(2:q+2:end) =  p ./ (knot(k+q)-knot(k));
      dN = dR*N;
    end
    N = R*N;
  end
  Ni( (i-1)*(p+1)+1:i*(p+1)) = mu-p:mu;
  Nj( (i-1)*(p+1)+1:i*(p+1)) = i*ones(1, p+1);
  Nv( (i-1)*(p+1)+1:i*(p+1)) = N;
  dNv((i-1)*(p+1)+1:i*(p+1)) = dN;
end
N  = sparse(Ni,Nj, Nv, n, numel(t));
dN = sparse(Ni,Nj,dNv, n, numel(t));

N  =  N(p+1:end-p,:); % remove extra functions from padding the knot vector
dN = dN(p+1:end-p,:);

