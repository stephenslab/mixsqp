#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Scale each column A[,i] by b[i].
void scalecols (mat& A, const vec& b) {
  unsigned int n = A.n_cols;
  for (unsigned int i = 0; i < n; i++)
    A.col(i) *= b[i];
}

// Normalize each row of A so that the entries in each row sum to 1.
void normalizerows (mat& A) {
  vec b = sum(A,1);
  A.each_col() /= b;
}

// Scale each row of A so that the largest entry in each row is 1.
void normalizerowsbymax (mat& A) {
  vec b = max(A,1);
  A.each_col() /= b;
}
