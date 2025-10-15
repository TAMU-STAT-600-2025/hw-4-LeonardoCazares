# Load the riboflavin data

# Uncomment below to install hdi package if you don't have it already; 
# install.packages("hdi") 
library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene expression
dim(riboflavin$x) # n = 71 samples by p = 4088 predictors
?riboflavin # this gives you more information on the dataset

# This is to make sure riboflavin$x can be converted and treated as matrix for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]


# Get matrix X and response vector Y
X = as.matrix(riboflavin$x)
Y = riboflavin$y

# Source your lasso functions
source("LassoFunctions.R")

# [ToDo] Use your fitLASSO function on the riboflavin data with 60 tuning parameters
fit60 <- fitLASSO(X, Y, n_lambda = 60, eps = 0.001) # Apply the fitLASSO on riboflavin data
lambda_seq <- fit60$lambda_seq # Get the lambdas
beta_mat   <- fit60$beta_mat # Get the betas
beta0_vec  <- fit60$beta0_vec # # Get the interceptions

# [ToDo] Based on the above output, plot the number of non-zero elements in each beta versus the value of tuning parameter
nnz <- colSums(beta_mat != 0)
plot(lambda_seq, nnz, type = "b",
     xlab = "lambda", ylab = "Number of nonzeros in beta",
     main = "LASSO sparsity path on riboflavin")

# [ToDo] Use microbenchmark 10 times to check the timing of your fitLASSO function above with 60 tuning parameters
library(microbenchmark) # Import microbenchmark
set.seed(2025)
mb <- microbenchmark(
  fitLASSO(X, Y, n_lambda = 60, eps = 0.001),
  times = 10 # Ten repetitions
)
print(mb) # Time statistics in miliseconds

mean_sec <- mean(mb$time) / 1e9 # Median in seconds
cat("Mean time (seconds):", round(mean_sec, 4), "\n")

# [ToDo] Report your median timing in the comments here: (~5.8 sec for Irina on her laptop)

median_sec <- median(mb$time) / 1e9 # Median in seconds
cat("Median time (seconds):", round(median_sec, 4), "\n")

# [ToDo] Use cvLASSO function on the riboflavin data with 30 tuning parameters (just 30 to make it faster)

set.seed(2025)  # CV folds randomized
cv30 <- cvLASSO(X, Y, n_lambda = 30, k = 5, eps = 0.001) # (5-folds and 30 lambdas)

# [ToDo] Based on the above output, plot the value of CV(lambda) versus tuning parameter. Note that this will change with each run since the folds are random, this is ok.

