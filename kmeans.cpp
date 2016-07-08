#include "RcppMLPACK.h"

using namespace mlpack::kmeans;
using namespace mlpack::metric;
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends(RcppMLPACK)]]

// [[Rcpp::export]]
List mlKmeans(const arma::mat& data, const int& clusters) {
  
  // size_t is the max number of iterations
  arma::Col<size_t> assignments;
  
  // The centroids will be stored in this matrix.
  arma::mat centroids;
  
  //mlpack::kmeans::RefinedStart::RefinedStart(const size_t samplings = 100 const double 	percentage = 0.02)	
  
  // The initialized Mahalanobis distance.
  
  // extern mlpack::metric::MahalanobisDistance distance;
  
  // Initialize with the default arguments.
  // More specification available
  // Like 
  //KMeans<mlpack::metric::MahalanobisDistance> k(1000, 1.0, distance);
  
  // this one works
  //KMeans<mlpack::metric::EuclideanDistance, RandomPartition> k(1000);
  
  
  
  KMeans<mlpack::metric::EuclideanDistance, RefinedStart> k(1000);
  
  
  k.Cluster(data, clusters, assignments, centroids); 
  
  return List::create(_["clusters"] = clusters,
                      _["result"]   = assignments,
                      _["centroids"] = centroids);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#mlKmeans(t(trees), 3)
*/

