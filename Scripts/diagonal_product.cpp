#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// Function to compute only the diagonal elements of a product of two square sparse matrtices in CSC format
// [[Rcpp::export]]
Eigen::VectorXd diagonal_producto_csc(const Eigen::SparseMatrix<double>& A, const Eigen::SparseMatrix<double>& B) {
  // Assure that both matrices have the same dimensions
  if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()) {
    throw std::invalid_argument("Matrices must be square and have the same dimensions.");
  }
  
  int n = A.rows();
  Eigen::VectorXd diagonal(n);
  diagonal.setZero(); // Initial values for the diagonal
  
  // Iteration over A columns
  for (int k = 0; k < n; ++k) {
    // Iterations over the non-zero elements of the k-th column of A
    for (Eigen::SparseMatrix<double>::InnerIterator it_A(A, k); it_A; ++it_A) {
      int i = it_A.row(); // Row of A in the k-th column
      
      // Iterations over the non-zero elements of the k-th row of B
      for (Eigen::SparseMatrix<double>::InnerIterator it_B(B, k); it_B; ++it_B) {
        int j = it_B.row(); // Row of B in the k-th column
        
        if (i == j) {
          diagonal(i) += it_A.value() * it_B.value();
        }
      }
    }
  }
  
  return diagonal;
}