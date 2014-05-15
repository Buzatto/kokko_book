
# Example boxes from Kokko's (2007) 'Modelling for Field Biologists...'
# adapted for R. Most instructions lifted verbatim from the text.


#####------------------------- CHAPTER 2 -------------------------#####

rm(list=ls())  # Clear all

## BOX 2.1 ##

# Create a vector with linearly equally spaced values
x <- seq(from = 0, to = 1, length = 101)

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
x <- seq(from = 0, to = 1, length = 101)

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
plot(x, xnew, type = 'l',
     xlab = 'Frequency of A in generation 0',
     ylab = 'Frequency of A males among mating males')
lines(x, x)

# Converting all this jazz into a function

# Function 'sexconflict': Allele frequency change in one generation.
# Give m, b, N, B as per definitions in Ch. 2 and accuracy
# as an integer that specifies how many values between 0 and 1
# the variable x should take. 
# Function returns x and xnew, and plots the results if plot = TRUE.

sexconflict <- function(m, b, N, B, accuracy, plot = FALSE){
  
  x <- seq(from = 0, to = 1, length = accuracy)
  
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
    plot(x, xnew, type = 'l',
         xlab = 'Frequency of A in generation 0',
         ylab = 'Frequency of A males among mating males')
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

# TOM:I'm not quite sure why this isn't exactly 0, maybe a quirk of how
# R handles numbers like this? I don't know enough about it. Suggestions welcome!


#####------------------------- CHAPTER 3 -------------------------#####

rm(list=ls())  # Clear all

## BOX 3.1 ##

# Function 'distr_shift'. A graphical explanation for the 
# shifting distributions, given the additive genetic variance,
# the initial mean of a distribution, and a and b of the function
# r(z) = a(z - b).
# Outputs are the magnitude of the shift, and the old and new
# distributions f0 & f1, which can be plotted against z.

distr_shift <- function(variance, initial_mean, a, b, plot = TRUE) {

  # We create 201 different values of z that lie within approximately
  # 12 standard deviations from the mean (6 for each side, to be similar to the 
  # range in the example). The 6.324555 (instead of just 6) is to make the grid z
  # be exactly like the one in the example, I arrived at that number through 
  # 2 / stand.d, because I wanted exactly two units to each side of the mean of 15
  # (to get from 13 to 17)
  stand.d <- sqrt(variance)
  z <- seq(from = (initial_mean - 6.324555*stand.d), to = (initial_mean + 6.324555*stand.d),
           length = 201)
  
  # The following is the density function of the normal distribution (simplified in R)
  f <- dnorm(z, mean=initial_mean, sd=sqrt(variance))
  
  # Then normalise it to 1. This might not work initially as we only have a selection
  # of discrete values of z, not 'all possible values'
  f0 <- f/sum(f)
  r <- a * (z - b)
  f1 <- f0 * exp(r) ## New distribution
  f1 <- f1 / sum(f1) ## Normalise
  
  new_mean <- sum(f1 * z)
  shift <- new_mean - initial_mean
  
  # Bundle up the output before plotting as it makes plotting vertical lines a bit easier
  output <- as.data.frame(cbind(shift, z, f0, f1))
  
  # Plot it (this got a bit out of hand...)
    if(plot == TRUE){
       plot(z, f0, type = 'l',
         xlab = 'Body size (z)',
         ylab = 'Frequency f(z)',
         ylim = c(0, max(f0) * 1.2))  # add a little extra room for lines atop the curves
       lines(z, f1)
       
       # Add lines & text atop curves showing the magnitude of the shift
       xf0 <- output$z[output$f0 == max(output$f0)]  # x coordinate for first line
       xf1 <- output$z[output$f1 == max(output$f1)]  # x coordinate for second line
       segments(x0 = xf0, x1 = xf0, y0 = max(f0), y1 = max(f0) * 1.1)  # line on first curve
       segments(x0 = xf1, x1 = xf1, y0 = max(f1), y1 = max(f1) * 1.1)  # line on second curve
       text(x = (xf0 + xf1)/2, y = max(f0) * 1.18, labels = round(shift, 1))  # Add text
       
    }
    
  return(shift)
}

# Test it
huzzah <- distr_shift(0.1, 15, 1, 15) 


# Now an alternative without normalizing, plotting numbers of individuals,
# and therefore showing population growth. Initial population size is necessary
# here.

# Function 'pop_shift'

pop_shift <- function(pop_size, variance, initial_mean, a, b, plot = TRUE) {
  
  # We create 201 different values of z that lie within approximately
  # 12 standard deviations from the mean (6 for each side, to be similar to the 
  # range in the example). The 6.324555 (instead of just 6) is to make the grid z
  # be exactly like the one in the example, I arrived at that number through 
  # 2 / stand.d, because I wanted exactly two units to each side of the mean of 15
  # (to get from 13 to 17)
  stand.d <- sqrt(variance)
  z <- seq(from = (initial_mean - 6.324555*stand.d), to = (initial_mean + 6.324555*stand.d),
           length = 201)
  
  # The following is the density function of the normal distribution (simplified in R)
  f <- dnorm(z, mean=initial_mean, sd=sqrt(variance))
  f0 <- f/sum(f)
  N0 <- f0*pop_size
  r <- a * (z - b)
  f1 <- f0 * exp(r) ## New distribution
  f1 <- f1 / sum(f1) ## Normalise
  N1 <- N0*exp(r)
  new_pop_size <- sum(N1) # get new population size
  
  new_mean <- sum(N1 * z) / new_pop_size
  shift <- new_mean - initial_mean
  
  # Bundle up the output before plotting as it makes plotting vertical lines a bit easier
  output <- as.data.frame(cbind(shift, z, N0, N1))
  
  # Plot it (this got a bit out of hand...)
  if(plot == TRUE){
    plot(z, N0, type = 'h',
    xlab = 'Body size (z)',
    ylab = 'Numbers of individuals',
    ylim = c(0, max(N1) * 1.2)) # add a little extra room for lines atop the curves
    points(z, N1, type = 'h', col="red")
    
    # Add lines & text atop curves showing the magnitude of the shift
    xN0 <- output$z[output$N0 == max(output$N0)] # x coordinate for first line
    xN1 <- output$z[output$N1 == max(output$N1)] # x coordinate for second line
    segments(x0 = xN0, x1 = xN0, y0 = max(N0), y1 = max(N1)*1.1) # line on first curve
    segments(x0 = xN1, x1 = xN1, y0 = max(N1), y1 = max(N1)* 1.1) # line on second curve
    text(x = (xN0 + xN1)/2, y = max(N1)*1.15, labels = round(shift, 1)) # Add text
    arrows(x0 = xN0, x1 = xN1, y0 = max(N1)*1.1, y1 = max(N1)*1.1, length = 0.1, lwd = 2)
  }
  
  return(shift)
}

# Test it
pop_shift(10000, 0.1, 15, 1, 15)


## BOX 3.2 ##

# Function 'barnacle'. 'Params' should contain b0, a0, a1, b1
# arranged in this weird order to make it easy to remember that
# the sequence of values should be decreasing. 
# addative_var is the amount of additive genetic variation.
# TOM: I've made 'generations' an argument, for fun.
# Output includes the mean value of z over time, and the 
# predicted equilibrium value (a0 - b0) / (b1 - b0 + a0 - a1).

barnacle <- function(params, additive_var, generations, plot = TRUE)
{
  
  z <- 0.1
  b0 <- params[1]
  a0 <- params[2]
  a1 <- params[3]
  b1 <- params[4]
  
  for (t in 1:(generations - 1)) {
    
    # Equation 3.4
    W <- z[t] * (z[t] / 2 * b1 + (1 - z[t] / 2) * b0) + (1 - z[t]) * ((1 + z[t]) / 2 * a1 + (1 - z[t]) / 2 * a0)
    
    # Equation 3.6
    deltaz <- additive_var / W * (z[t] * (b1 - a1) + (1 - z[t]) * (b0 - a0))
    
    # The new z is the old one plus the change that occurred
    z[t+1] <- z[t] + deltaz
    
  }
  
  if(plot == TRUE){
     plot(x = 1:generations, y = z, type = 'l',
          xlab = 'Generations', 
          ylab = 'Switchpoint (z)',
          ylim = c(0, 1))
  }
  
  equilibrium <- (a0 - b0) / (b1 - b0 + a0 - a1)
  
  output <- as.data.frame(cbind(z, equilibrium))
    
  return(output)
  
}

# Test it
squish <- barnacle(params = c(1.1, 1.05, 1.02, 1), additive_var = 1, generations = 1000)


#####------------------------- CHAPTER 4 -------------------------#####

rm(list=ls())  # Clear all









