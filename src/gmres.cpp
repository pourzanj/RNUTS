#include <Rcpp.h>
using namespace Rcpp;

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

// The code for this was taken from https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html
class MatrixReplacement;
using Eigen::SparseMatrix;
namespace Eigen {
namespace internal {
// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
template<>
struct traits<MatrixReplacement> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
{};
}
}
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };
  Index rows() const { return size_; }
  Index cols() const { return size_; }
  template<typename Rhs>
  Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
  // Custom API:
  MatrixReplacement(int size, Function Av) : size_(size), Av_(Av) {}

  Function get_Av() const { return Av_; }
private:
  int size_;
  Function Av_;
};
// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
template<typename Rhs>
struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixReplacement,Rhs,generic_product_impl<MatrixReplacement,Rhs> >
{
  typedef typename Product<MatrixReplacement,Rhs>::Scalar Scalar;
  template<typename Dest>
  static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
  {
    // This method should implement "dst += alpha * lhs * rhs" inplace,
    // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
    Function Av = lhs.get_Av();
    NumericVector rhs_vector(rhs.size());
    for(int i = 0; i < rhs.size(); i++) {
      rhs_vector[i] = rhs(i);
    }

    NumericVector out = Av(rhs_vector);

    for(int i = 0; i < dst.size(); i++)
      dst(i) += alpha * out(i);
  }
};
}
}


//' Solves Ax=b for x given a function to compute the product Av.
//'
//' Error returned is relative and is defined as |Ax-b| / |b| where
//' the two-norm is used.
//'
//' @param Av
//' @param b
//' @param x0
//' @param TOL
//' @export
// [[Rcpp::export]]
List gmres(Function Av, NumericVector b, NumericVector x0, double TOL) {

  int D = b.size();

  // declare and set Eigen versions of input
  MatrixReplacement A(D, Av);
  Eigen::VectorXd b_eigen(D), x0_eigen(D);

  for(int i = 0; i < D; i++) {
    b_eigen(i) = b[i];
    x0_eigen(i) = x0[i];
  }

  // solve using Eigen's GMRES implementation
  Eigen::GMRES<MatrixReplacement, Eigen::IdentityPreconditioner> solver;
  solver.setTolerance(TOL);
  solver.compute(A);
  Eigen::VectorXd x_eigen = solver.solveWithGuess(b_eigen, x0_eigen);

  // declare outputs and convert Eigen outputs back to Rcpp types and return
  NumericVector x(D);
  Rcpp::List out;

  for (int i = 0; i < D; ++i) {
    x(i) = x_eigen(i);
  }

  out["x"] = x;
  out["iterations"] = solver.iterations();
  out["error"] = solver.error();

  return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# A <- matrix(c(-0.594235967,  0.306972829, -0.006492769,
#               1.413205331,  0.479708437, -0.709804312,
#               0.308042007,  1.460280233, -1.115569795), nrow = 3, byrow = TRUE)
# b <- 1:3
# Av <- function(v) A %*% v
# x0 <- rep(0.0, 3.0)
#
# gmres(Av, b, x0, TOL = 1e-15)
*/
