
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

rm(list=ls()) # Clear all


## BOX 4.1 ##

Trait  <- seq(0,2, length=101)
b <- Trait
cL <-  2*Trait^2; cH <- Trait^2;

## Plotting figure 4.2
plot(Trait, b, frame.plot=F, 
     xlab="Male trait, T", ylab="Benefits and costs",
     type="l", ylim=c(0,8), lwd=2) # plots the benefits
points(Trait, cL, type="l") # adds the costs for low quality males
points(Trait, cH, type="l") # adds the costs for high quality males

# adds letters to identify the lines/curves
text(x=2, y=2.2, expression(bold(italic("b"))))
text(x=2, y=4.2, expression(bold(italic("c"))[H]))
text(x=2, y=8.2, expression(bold(italic("c"))[L]))


## Plotting figure 4.3
WL <-  b-cL; WH <- b-cH;
plot(Trait, WL, frame.plot=F, 
     xlab="Male trait, T", ylab="Net benefit",
     type="l", ylim=c(0,0.5), xlim=c(0,1),lwd=2) # plots WL
points(Trait, WH, type="l", lwd=2) # adds WH
grid(10,10) #adds a grid, Hanna mentions that in the book, but the grid isn't in Fig 4.3.

## bind Trait, WL and WH together, to find maxima and draw arrows
find.max <- data.frame(Trait, WL, WH)
arrows(x0 = mean(find.max$Trait[find.max$WL==max(WL)]),
       y0 = 0.08, y1 = max(WL), length = 0.1, lwd = 2)
arrows(x0 = mean(find.max$Trait[find.max$WH==max(WH)]),
       y0 = 0.15, y1 = max(WH), length = 0.1, lwd = 2)

text(x=mean(find.max$Trait[find.max$WL==max(WL)]),
     y=0.05, expression(bold(italic("W"))[L]))
text(x=mean(find.max$Trait[find.max$WH==max(WH)]),
     y=0.12, expression(bold(italic("W"))[H]))

# Hanna also included Matlab code to zoom in the plot, but you can do that
# in R by just increasing the window with the mouse.

## Plotting figure 4.4 (not in Box 4.1, but still useful!)

optimal.T <- function(alpha) {1/(2*alpha)}
alpha <- seq(0,2, length=101)

plot(alpha, optimal.T(alpha), type="l", frame.plot=F, ylim=c(0,20),
     xlab=expression(paste("Cost coefficient, ")~symbol(a)), 
     ylab=expression(paste("Optimal trait value,")~italic(T)~paste("*")))


## BOX 4.2 ##

maledisplay <- function(resource.range, resolution, plot=T) {
  
  Trait <- numeric(0)
  Seasons <- numeric(0)
  Fit <- numeric(0)
  R <- seq(resource.range[1], resource.range[2], length = 100)
  
  for(i in 1:length(R)) {
    trait <- seq(0, R[i], length=resolution)
    net.condition <- R[i] - trait
    survival <- 1 / (1 + exp(-10*(net.condition - 0.5)))
    future.mating.seasons <- survival / (1 - survival)
    m <- trait
    fitness <- m + m*future.mating.seasons
    finding.max <- data.frame(trait, net.condition, survival, future.mating.seasons, m, fitness)
    best.fitness <- finding.max[finding.max$fitness==max(finding.max$fitness),]
    
    #storing results
    Trait[i] <- as.numeric(best.fitness[1,1])
    Seasons[i] <- as.numeric(1 + best.fitness[1,4])
    Fit[i] <- as.numeric(best.fitness[1,6])
  }
  
  output <- data.frame(cbind(R, Trait, Seasons, Fit))

  if(plot == TRUE){
    par(mfrow=c(3,1))
    plot(output$R, output$Trait, type = 'l', xlab = '', ylab = 'Male trait', frame.plot=F)
    plot(output$R, output$Seasons, type = 'l', xlab = '', ylim=c(0, max(output$Seasons)+1),
         ylab = 'Expected number of mating seasons', frame.plot=F)
    plot(output$R, output$Fit, type = 'l',
         xlab = expression(paste("Male resources,")~italic(R)),
         ylab = 'Fitness', frame.plot=F)
  }
  
  return(output)
  
}

maledisplay(c(0, 0.2), resolution=1000) # generates Figure 4.7

maledisplay(c(0.8, 1), resolution=1000) # generates Figure 4.8
maledisplay(c(0.8, 1), resolution=10000) # smoother version of Figure 4.8

maledisplay(c(0, 1), resolution=1000) # generates Figure 4.9
