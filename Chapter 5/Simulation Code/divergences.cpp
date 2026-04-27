#include <Rcpp.h>
#include <cmath>

/// this file was created with the help of Gemini
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Cjeffreys(NumericMatrix data) {
  int n = data.nrow();
  int p = data.ncol();

  // 1. Pre-calculate logs once to avoid repeating log() calls in the loops
  NumericMatrix log_data(n, p);
  for(int i = 0; i < n; i++) {
    for(int k = 0; k < p; k++) {
      log_data(i, k) = std::log(data(i, k));
    }
  }

  // 2. Prepare output matrix
  NumericMatrix D(n, n);

  // 3. Compute distances (Lower triangle + diagonal)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      double total = 0;
      for (int k = 0; k < p; k++) {
        total += (data(i, k) - data(j, k)) * (log_data(i, k) - log_data(j, k));
      }
      D(i, j) = total;
      D(j, i) = total; // Mirror to upper triangle
    }
  }

  return D;
}


// [[Rcpp::export]]
NumericMatrix Cjensensshannon(NumericMatrix data) {
  int n = data.nrow();
  int p = data.ncol();

  // 1. Pre-calculate p * log(p)
  NumericMatrix data_log_data(n, p);
  for(int i = 0; i < n; i++) {
    for(int k = 0; k < p; k++) {
      data_log_data(i, k) = data(i,k) * std::log(data(i, k));
    }
  }

  // 2. Prepare output matrix
  NumericMatrix D(n, n);


  // 3. Compute distances (Lower triangle + diagonal)
  double logm = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      double total = 0;
      for (int k = 0; k < p; k++) {
        logm = std::log(0.5 * (data(i, k) + data(j, k)));
        total += data_log_data(i, k) + data_log_data(j, k) - logm * (data(i, k) + data(j, k));
      }
      D(i, j) = 0.5 * total;
      D(j, i) = 0.5 * total; // Mirror to upper triangle
    }
  }

  return D;
}

// [[Rcpp::export]]
NumericMatrix Creny(NumericMatrix data, double alpha) {
  int n = data.nrow();
  int p = data.ncol();

  // 2. Prepare output matrix
  NumericMatrix D(n, n);


  // 3. Compute distances (Lower triangle + diagonal)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double total = 0;
      for (int k = 0; k < p; k++) {
        total += std::pow(data(i, k), alpha) / std::pow(data(j, k), (alpha - 1));
      }
      D(i, j) =  std::log(total) / (alpha - 1);
    }
  }

  return D;
}

// [[Rcpp::export]]
NumericMatrix Ckaniadakis(NumericMatrix data) {
  int n = data.nrow();
  int p = data.ncol();

  // 1. Pre-calculate logk(p) and A
  NumericMatrix logk_data(n, p);
  NumericMatrix A(n, p);
  for(int i = 0; i < n; i++) {
    double totalA = 0;
    for(int k = 0; k < p; k++) {
      logk_data(i, k) = 0.5 * (data(i, k) - 1) / data(i, k);
      A(i, k) = (2 * std::pow(data(i, k ), 2)) / (1 + std::pow(data(i, k ), 2));
      totalA += A(i, k);
    }
    // normalize A
    for(int k = 0; k < p; k++) {
      A(i, k) = A(i, k) / totalA;
    }
  }

  // 2. Prepare output matrix
  NumericMatrix D(n, n);


  // 3. Compute distances
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double total = 0;
      for (int k = 0; k < p; k++) {
        total += (logk_data(i, k) - logk_data(j, k)) * A(i, k);
      }
      D(i, j) = total;
    }
  }

  return D;
}
