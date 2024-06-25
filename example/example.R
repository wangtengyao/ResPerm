vector.normalise <- function(v){v / sqrt(sum(v^2))}
n <- 600; p <- 100; k <- 5; coeff_sig <- 1; sigma <- 1; AR <- 0; 
coeff_shape='gaussian'
X <- matrix(rnorm(n*p), n)
alpha <- vector.normalise(c(rnorm(k), rep(0, p-k)))
beta <- vector.normalise(c(rnorm(k), rep(0, p-k)))
e <- rt(n, df=2)
eps <- rt(n, df=2)

b <- 0.5 # coefficient of interest
Z <- X %*% alpha + e
Y <- X %*% beta + b * Z + eps

X_aug <- cbind(Z, X)
NaivePermutationTest(X_aug, Y, K=99)
ResidualPermutationTest(X_aug, Y, K=99)