
# Example boxes from Kokko's (2007) 'Modelling for Field Biologists...'
# adapted for R. Most instructions lifted verbatim from the text.


#####------------------------- CHAPTER 2 -------------------------#####

rm(list=ls())  # Clear all

## BOX 2.1 ##

# Create a vector with linearly equally spaced values
x <- seq(from = 0, to = 1, by = 1/101)

# Choose a value of m to be investigated
m <- 5

# Calculate y given x
y <- m * x/(m * x + (1 - x))

# Plot 'er
plot(x, y, type = 'l',
     xlab = 'Proportion of A males in the population',
     ylab = 'Proportion of A males among mating males')

# Repeat as above for three values of m
m1 <- 1.5
m2 <- 5
m3 <- 10

y1 <- m1 * x/(m1 * x + (1 - x))
y2 <- m2 * x/(m2 * x + (1 - x))
y3 <- m3 * x/(m3 * x + (1 - x))

# Plot in one 
plot(x, y1, type = 'l',
     xlab = 'Proportion of A males in the population',
     ylab = 'Proportion of A males among mating males')
lines(x, y2)
lines(x, y3)


## BOX 2.2 ##

# Make a vector of 101 values from 0 to 1
x <- seq(from = 0, to = 1, by = 1/101)

# Choose a value of m to be investigated
m <- 5

# Calculate how y depends on x
y <- (m * x)/(m * x + 1 - x)

# Calculate mating proportions
pAA <- x * y
pAa <- x * (1 - y)
paA <- (1 - x) * y
paa <- (1 - x) * (1 - y)

# Check these sum to 1
pAA + pAa + paA + paa

# Now to nA and na. We must also choose values for N, B, and b
N <- 1000
B <- 20
b <- 0.8

nA <- N * (pAA * b * B + 0.5 * pAa * B + 0.5 * paA * b * B)
na <- N * (0.5 * pAa * B + 0.5 * paA * b * B + paa * B)

# Then we get xnew, and plot it against x
xnew <- nA / (nA + na)

# Plot xnew against x, then x against x
plot(x, xnew, type = 'l')
lines(x, x)

# Converting all this jazz into a function

# Function 'sexconflict': Allele frequency change in one generation.
# Give m, b, N, B as per definitions in Ch. 2 and accuracy
# as an integer that specifies how many values between 0 and 1
# the variable x should take. 
# Function returns x & xnew & plots the results if plot = TRUE.

sexconflict <- function(m, b, N, B, accuracy, plot = FALSE){
  
  x <- seq(from = 0, to = 1, by = 1/accuracy)
  
  # Calculate  how y depends on x
  y <- (m * x)/(m * x + 1 - x)
  
  # Mating proportions
  pAA <- x * y
  pAa <- x * (1 - y)
  paA <- (1 - x) * y
  paa <- (1 - x) * (1 - y)
  
  nA <- N * (pAA * b * B + 0.5 * pAa * B + 0.5 * paA * b * B)
  na <- N * (0.5 * pAa * B + 0.5 * paA * b * B + (paa * B))
  
  xnew <- nA / (nA + na)
  
  if(plot == TRUE){
    plot(x, xnew, type = 'l')
    lines(x, x)
  }
  
  # Package up the output
  output <- as.data.frame(cbind(x, xnew))
  
  # Send it on back
  return(output)
  
}

# Check it works
woop <- sexconflict(m = 5, b = 0.9, N = 200, B = 20, accuracy = 101, plot = TRUE)


## BOX 2.3 ##

# A quick check to see that population size N does not seem to matter:
N_case1 <-  1
xnew_case1 <- sexconflict(5, 0.5, N_case1, 20, 101, plot = TRUE)
N_case2 <- 100000
xnew_case2 <- sexconflict(5, 0.5, N_case2, 20, 101, plot = TRUE)

# To check they are identical, we can subtract the solutions of the two cases from
# each other and look at the max absolute value of the difference between them.
# The difference should be ~ 0.

max(abs(xnew_case1 - xnew_case2))

# I'm not quite sure why this isn't exactly 0, maybe a quirk of how
# R handles numbers like this? I don't know enough about it...
# It's close enough to make the point anyway.


#####------------------------- CHAPTER 3 -------------------------#####

rm(list=ls())  # Clear all