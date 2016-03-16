#include <mex.h>
#include <matrix.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs < 3)
    mexErrMsgTxt("Too few arguments");
  double *knot_in = mxGetPr(    prhs[0]);
  int    p        = mxGetScalar(prhs[1]);
  double *t       = mxGetPr(    prhs[2]);
  int    n        = mxGetN(prhs[0])*mxGetM(prhs[0]) - p - 1; // number of basis functions 
  int    m        = mxGetN(prhs[2])*mxGetM(prhs[2]);         // number of evaluation points
  // ensure open knot vector by appending with 2p knots (remove these at the end)
  vector<double> knot(n+p+1 + 2*p);
  fill(knot.begin(),   knot.begin()+p, knot_in[0]);
  copy(knot_in, knot_in + (n+p+1), knot.begin()+p);
  fill(knot.end()-p,  knot.end(), knot_in[n+p]);

  // allocate memory for result sparse matrix
  vector<int>     Ni(m*(p+1), 0);
  vector<int>     Nj(m+1    , 0);
  vector<double>  Nv(m*(p+1), 0.0);
  vector<double> dNv(m*(p+1), 0.0);
  vector<int>::iterator     iNi  = Ni.begin();
  vector<int>::iterator     iNj  = Nj.begin();
  vector<double>::iterator  iNv  = Nv.begin();
  vector<double>::iterator  idNv = dNv.begin();

  for(int j=0; j<m; ++j) {

    vector<double>::iterator mu_it = lower_bound(knot.begin(), knot.end(), t[j]);
    if(t[j] == knot.back())   // evaluation in the limit from below
      mu_it  = lower_bound(knot.begin(), knot.end(), t[j]);
    else                      // evaluation in the limit from above
      mu_it  = upper_bound(knot.begin(), knot.end(), t[j]);
    if(mu_it == knot.end()) { // evaluation outside domain
      iNj[1] = iNj[0];
      iNj++;
      continue;
    }
    int mu = ((int)(mu_it-knot.begin()))-1;

    vector<double>  N(p+1, 0.0);
    vector<double> dN(p+1, 0.0);
    N[p] = 1;
    for(int q=1; q<=p; ++q) {
      for(int k=p-q; k<p+1; ++k) {   // k is local index (0,p+1)
        int i = mu-p+k;              // i is global index (0,n)

        if( k != p-q ) {
          if( p==q )
            dN[k] = N[k]   *                p   / (knot[i+q]-knot[i]  );
          N[k]  =   N[k]   * (t[j] - knot[i]  ) / (knot[i+q]-knot[i]  );
        }

        if( k != p ) {
          if( p==q )
            dN[k]+= N[k+1] *              (-p)  / (knot[i+q+1]-knot[i+1]);
          N[k] +=   N[k+1] * (knot[i+q+1]-t[j]) / (knot[i+q+1]-knot[i+1]);
        }

      }


    }
    vector<double>::iterator iN  =  N.begin();
    vector<double>::iterator idN = dN.begin();
    int nnz = 0;
    for(int i=mu-p; i<=mu; ++i, ++iN, ++idN) {
      if(i>=p && i<n+p && (*iN != 0 || *idN != 0)) {
        *iNi++  = i-p;
        *iNv++  = *iN;
        *idNv++ = *idN;
        ++nnz;
      }
    }
    iNj[1] = iNj[0] + nnz;
    ++iNj;
  }

  if(nlhs > 0) {
    plhs[0] = mxCreateSparse(n,m, m*(p+1), mxREAL);
    size_t *ir = mxGetIr(plhs[0]);
    size_t *jc = mxGetJc(plhs[0]);
    double *pr = mxGetPr(plhs[0]);
    copy(Ni.begin(), Ni.end(), ir);
    copy(Nj.begin(), Nj.end(), jc);
    copy(Nv.begin(), Nv.end(), pr);
    if(nlhs > 1) {
      plhs[1] = mxCreateSparse(n,m, m*(p+1), mxREAL);
      ir = mxGetIr(plhs[1]);
      jc = mxGetJc(plhs[1]);
      pr = mxGetPr(plhs[1]);
      copy(Ni.begin(),  Ni.end(),  ir);
      copy(Nj.begin(),  Nj.end(),  jc);
      copy(dNv.begin(), dNv.end(), pr);
    }
  }

}
