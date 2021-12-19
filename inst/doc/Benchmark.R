## ----eval=FALSE---------------------------------------------------------------
#  gibbsR <- function(N, a, b, n, thin) {
#    mat <- matrix(nrow = N, ncol = 2)
#    x <- floor(n/2)
#    y <- 0.5
#    for (i in 1:N) {
#      for (j in 1:thin) {
#        x <- rbinom(1, size = n, prob = y)
#        y <- rbeta(1, x + a, n - x + b)
#      }
#      mat[i, ] <- c(x, y)
#    }
#    mat
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix gibbsC(int N, double a, double b, int n, int thin) {
#    NumericMatrix mat(N, 2);
#    double x = floor(n/2), y = 0.5;
#    for(int i = 0; i < N; i++) {
#      for(int j = 0; j < thin; j++) {
#        x = rbinom(1, n, y)[0];
#        y = rbeta(1, x+a, n-x+b)[0];
#      }
#      mat(i, 0) = x;
#      mat(i, 1) = y;
#    }
#    return(mat);
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp21017)
library(microbenchmark)
tm1 <- microbenchmark::microbenchmark(
  rnR = gibbsR(2000,1,4,50,10),
  rnC = gibbsC(2000,1,4,50,10)
)

knitr::kable(summary(tm1)[,c(1,3,5,6)])

