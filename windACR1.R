#' #Methods of Time Series Analysis Applied to Weather Measurements on a Wind Farm#

#' ## Introduction ##
#' The chosen source of data is a wind farm located in the North Sea
#'  (with coords: 54.01433°NL, 6.587667°EL) denoted as LES (1y 2015-01-01 2015-12-31 cfsr) 
#'  which comes from website http://www.vortexfdc.com/, specifically designed to provide its
#'  users with weather data for potential wind farm sites. This specific dataset has been 
#'  obtained (for free) from year 2015, and contains a variety of measurements, ranging from
#' wind velocity magnitude, through the wind direction angle, temperature up to dimensionless 
#' characteristics like  the Richardson number "Ri" which corresponds to the ratio of buoyancy
#' over shear flow.
#' 
#' In the first 7 parts of this project, I focus on the analysis of a single variable, namely
#' temperature T [°C]. Since the dataset consist of over 47 000 observations of 11 different 
#' variables, and thus might be difficult to visualize, I've chosen to include only 
#' the temperature measurements for each day at 12am.
#' 
#'---------------------------------------------------------------------------------------
#' 
#' ## 1. Visualization and Basic Stats ## 
#' We begin by extracting the data from a downloaded file

temperature <- read.table("vortex_TS.txt", header = T, skip = 3, 
                          colClasses = "character")

#' and tagging the 8th column with a mark "T" for temperature
names(temperature)[8] <- c("T")
#' then separating the time data into separate columns using the Tidyr package

temperature <- tidyr::separate(data = temperature, col = YYYYMMDD, 
                               into = c("year","month","day"), sep=c(4,6))
temperature <- tidyr::separate(data = temperature, col = HHMM, 
                               into = c("hour","minute"), sep=2)

#' and convert the time columns to numbers:

temperature <- as.data.frame(apply(temperature, 2, as.numeric))

#' Create a separate array with columns: "year","month","day","hour","minute","T"
#' filtering only the temperature at 12am each day:

temp <- temperature[temperature$hour==12 & temperature$minute==0,
                    c("year","month","day","hour","minute","T")] 

#' and uniting the time columns saving them as "YYYY-MM-DD-HH-MM" dates

temp <- tidyr::unite(temp, col=time, year, month, day, hour,
                     minute, sep="-", remove=T)
temp$time <- strptime(temp$time, format="%Y-%m-%d-%H-%M")
lastDate <- temp$time[length(temp$time)]
temp$time <- as.numeric(difftime(temp$time, temp$time[1], units = "days"))

#'Now we can visualize the data as followed

```{r basicDataPlot1, fig.width=9, fig.height=6}
par(mfrow=c(1,1))
plot(T ~ temp$time, temp, type="p",xlab="day", ylab="T[°C]",main="Daily temperature at 12:00",axes=F)
lines(temp$time, temp$T,col="black")
axis(side=1, at=seq(1, 365, 7))
axis(side=2, at=seq(min(temp$T), max(temp$T), by=2))
box()

#' Now we can clearly see the periodic character of the time series which corresponds to the 
#' seasonal changes in this latitude. For more details about the particular measurements we might
#' calculate the basic characteristics of the set of values:

stat <- summary(temp$T)
stddev <- sd(temp$T, na.rm=TRUE)
v <- as.numeric(c(stat[1],stat[6],stat[4],stat[3],stddev))
names(v) <- c("Min.","Max.","Mean","Median","Std.Dev.")
v

#' We see that the minimum temperature does not drop far below zero even during the winter, 
#' possibly due to warm currents in the northern Atlantic Ocean.
#' Aside from local fluctuations, the same pattern can be expected to appear in the next 
#' year, and in the year after that. Perhaps after several decades we might see a change in the 
#' annual average temperature due to global warming, for instance.
#' 
#'-----------------------------------------------------------------------------------------
#'    
#' ## 2 Time Series Decomposition: Trend & Seasonal Components ##
#'
#' Now we proceed to make the first step in the analysis of the extracted temperature data. 
#' It is more or less clear that we will not be able to observe any linear or exponential trends
#' on such small dataset. 
#' We can verify the absence of a trend by testing the series with the Mann-Kendall rank test:

randtests::rank.test(temp$T)

#' Which suggests that there is no significant trend in the sample as a whole. We will return to the test in 
#' section 3. 
#'
#' The only clearly visible systematic components are the periodic seasonal
#' changes due to Earth's tilt with respect to the Ecliptic plane. 
#' Clearly, the most visible cycle will repeat with a period of 365 days. Other fractions of this cycle 
#' might be present as well. The remaining cycles can be observed by examining the ACF (Autocorrelation Function)
#' which can be plotted using just the `acf` command:

```{r basicACFPlot2, fig.width=9, fig.height=4}
par(mfrow=c(1,2))
acf(temp$T, lag.max=365, main="Daily temp ACF")
acf(temp$T, lag.max=365, type="partial",main="Daily temp PACF")

#' The second plot shows a partial autocorrelation function which also accounts for lags 
#' shorter than the given lag k, that is: `k-1` , `k-2`, ... , `2`, `1`.
#' With common sense, one deduces that the data already has a 365-day cycle, but from the
#' correlogram we see that other smaller cycles are present in the time series as well. 
#' Hence we model the cycles with the following periods

n <- length(temp$T)
seasons <- c(365, 365/2, 365/4, 365/12)
names(seasons) <- c("S1.","S2.","S3","S4")
seasons

#' We may assume that there will be a month-long cycles which correspond to the rotation of the Moon, rather
#' than the calendar we use:

months = c(
  "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"
)
```{r monthlyPlot, fig.width=10, fig.height=3}
par(mfrow=c(1,1))
for (month in 1:12) {
  from = ((month - 1)*floor(seasons[4]) + 1)
  if (month == 12) {
    to = n;
  } else {
    to = (month * floor(seasons[4]));
  }
  seg = from:to
  plot(x=1:length(seg), y=temp$T[seg], type="l",
       main=paste("12:00 Temp. in ",months[month]), sub=paste("from day", from, " to ", to), xlab="[day]", ylab="T[°C]", lwd=2, lty=1.5)
}

#' The monthly data do not seem to have a repeating pattern. Hence we will only consider the annual 365-day period.
#' 
#' From now on we will separate our time series into a test part and an evaluation part. 
#' All regression analysis will be performed on the test part. Predictions and their quality will
#' be examined on the evaluation part. Generally the test part is taken as the first 80% of the original
#' time series:

nt <- floor(n * 8 / 10)
temp_test <- list()
temp_test$T <- temp$T[1:nt]
temp_test$time <- temp$time[1:nt]

temp_eval <- list()
temp_eval$T <- temp$T[(nt + 1):n]
temp_eval$time <- temp$time[(nt + 1):n]

par(mfrow=c(1,1))
plot(x=temp$time, y=temp$T, type="p", main="Test and Evaluaion Parts", xlab="time", ylab="T[°C]")
lines(x=temp_test$time, y=temp_test$T, lwd=1.5, col="blue")
lines(x=temp_eval$time, y=temp_eval$T, lwd=1.5, col="green")
legend("topleft", legend=c("test","eval"),
       col=c("blue","green"), lty=1,lwd=2 , cex=0.8)

#' Now we will take the chosen season and regress them against the test part:

model.sCos <- lm(T ~ cos(2*pi*temp$time/seasons[1]) + sin(2*pi*temp$time/seasons[1]) #+
                     #cos(2*pi*temp_test$time/seasons[2]) + sin(2*pi*temp_test$time/seasons[2]) +
                     #cos(2*pi*temp_test$time/seasons[3]*2) + sin(2*pi*temp_test$time/seasons[3]*2) +
                     #cos(2*pi*temp_test$time/seasons[4]) + sin(2*pi*temp_test$time/seasons[4])
                   , temp)
summary(model.sCos)

#' The summary of the regression model implies that the annual period is significant.
#' And now we can plot the residues after extracting the seasonal components

temp_test$sCos <- model.sCos$fitted.values[1:nt]
temp_eval$sCos <- model.sCos$fitted.values[(nt + 1):n] # fitted values for later use

temp_test$Tres <- model.sCos$residuals[1:nt]
temp_eval$Tres <- model.sCos$residuals[(nt + 1):n] # residuals for later use

```{r residPlot1, fig.width=10, fig.height=4}
par(mfrow=c(1,2))
plot(T ~ temp_test$time, temp_test,
     main="Fitted Annual Temperature Model",xlab="day",ylab="T[°C]",axes=F)
lines(temp_test$time, temp_test$sCos, col="blue", lwd=2)
axis(side=1, at=seq(1, n, 7))
axis(side=2, at=seq(min(temp_test$T), max(temp_test$T), by=2))
lines(temp_test$time, temp_test$sCos, col="blue", lwd=2)
box()

plot(x=temp_test$time, y=temp_test$Tres,
     main="Residuals",xlab="day",ylab="res(T)[°C]",axes=F,type="l")
axis(side=1, at=seq(1, n, 7))
axis(side=2, at=seq(round(min(temp_test$Tres), digits=1), round(max(temp_test$Tres), digits=1), by=2))
box()

#' We can see that the residuals after extracting the seasonal components still haven't 
#' lost some of their periodic behavior. Taking a look at the residual ACF and PACF we see 
#' that non-zero correlation extends from lag of approximately 50 days.

```{r ACFPlot2, fig.width=9, fig.height=4}
par(mfrow=c(1,2))
acf(temp_test$Tres, lag.max=365, main="Residual ACF")
acf(temp_test$Tres, lag.max=365, main="Residual PACF", type="partial")

#' When examining the residual ACF we notice remaining significant correlation for multiple lag
#' values. This is possible due to remaining cyclical components which will be found as significant
#' frequencies after Fourier analysis in the next section.
#'
#'---------------------------------------------------------------------------------------
#'             
#' ## 3 Time Series Decomposition: Cyclical Components & Randomness Tests on Residuals ##
#' ### 3.1 Cyclical Components and Fourier Analysis on Time Series ###
#'
#' Now we continue with another step in the decomposition of the original time series.
#' We separate the remaining oscillations via discrete Fourier analysis. A continuous spectrum
#' is generated from the ACF using the following function:

SpecDens <- Vectorize(  
  FUN = function(omega, acov=NULL, data=NULL) {
    if(is.null(acov)) {
      if(is.null(data)) {
        stop("Please provide either vector of autocovariance function values or time series data.")
      }
      else acov <- acf(data, type="covariance", plot=FALSE)
    }
    k <- seq(to=length(acov)-1)
    ( acov[1] + 2*sum(acov[-1]*cos(k*omega)) ) / (2*pi)
  },
  vectorize.args = "omega"
)

```{r SpecPlot2, fig.width=9, fig.height=4}
par(mfrow=c(1,2))
temp_test.ACF <- as.numeric( acf(temp_test$Tres, type="covariance", plot=FALSE)$acf )
plot(function(x) SpecDens(x, acov=temp_test.ACF), from=2*pi/(nt),
     to = 2 * pi / 4, xlab="frequency", ylab="spectral density")
plot(function(x) SpecDens(2*pi / x, acov=temp_test.ACF), from=4, 
     to = nt, xlab="period", ylab="spectral density")

#' the figure on the right provides a better image of the (continuous) distribution of individual
#' components with periods of oscillation. The hidden frequencies are obtained using discrete
#' Fourier transform (more specifically Fast Fourier Transform):

```{r SpecDens,fig.width=11, fig.height=3.5}
temp_test.FFT <- abs(fft(temp_test$Tres)) # Fourier transform
omega <- 2 * pi * seq(0, nt / 2 - 1) / nt
period <- 2 * pi / omega
par(mfrow=c(1,1))
plot(period, temp_test.FFT[seq_along(period)], main="Residual Spectral Density", 
     ylab="density", xlab="period (days)", type="h", log="x")

#' which we can then sort by density, and display the ones with the most weight:

spectrum <- data.frame(f=temp_test.FFT[seq_along(omega)], omega=omega, period=period)
spectrum <- spectrum[order(spectrum$f, decreasing = TRUE), ]
head(spectrum, 10)

#' One reliable way to find the most significant frequencies is by using Fischer's Periodicity 
#' Test:

signif <- cbind(spectrum[0,], 
                  data.frame(test.stat = numeric(0), crit.val = numeric(0)))
spec <- spectrum
repeat {
  test <- data.frame(test.stat = spec$f[1]/sum(spec$f),
                     crit.val = 1 - (0.05/nrow(spec))^(1/(nrow(spec)-1)))
  if(test$test.stat > test$crit.val) {
    signif <- rbind(signif, cbind(spec[1,],test));
    spec <- spec[-1,]
  } else break
}
rm(spec, test)
signif

#' As it appears, Fischer's Test yields no significant frequency, but since we
#' observe a sequence of significant frequencies that correspond with the ACF shown earlier.
#' Nonetheless, we assume that, say, the first 10 might contribute to the oscillatory behavior
#' of the time series. Thus we take:

( signif <- head(spectrum, 10) )

#' which we can then determine more precisely (by finding their local maxima) using the continuous spectral density:

newsignif <- sapply(
  signif$omega,
  function(x) optimize(SpecDens, 
                       interval = x + c(1,-1) * pi / nt, 
                       acov = temp_test.ACF, 
                       maximum = T  )$maximum)
( newSignifs <- cbind(signif, data.frame(omega_sm = newsignif, 
                                       period_sm = 2 * pi / newsignif)) )


#' We can try to complete the regression using the first 3 periods with their fractions, but as
#' it turns out, only the first one with its half is significant.

( periods <- head(newSignifs, 3)$period_sm )

model.cCos <- lm(Tres ~ cos(2*pi*temp_test$time/periods[1]) + sin(2*pi*temp_test$time/periods[1]) +
                     cos(2*pi*temp_test$time/periods[1]*2) + sin(2*pi*temp_test$time/periods[1]*2) #+
                   #cos(2*pi*temp_test$time/periods[2]) + sin(2*pi*temp_test$time/periods[2]) +
                   #cos(2*pi*temp_test$time/periods[3]) + sin(2*pi*temp_test$time/periods[3])
                   ,
                 temp_test)
summary(model.cCos)

#' The cyclical components defined by first two periods seem to be significant, but we cannot be sure if they truly
#' belong into the systematic model. For that reason we will save the residues of this model so that we can compare its
#' predictive properties with only `sCos` in section 5.

temp_test$cCos <- model.cCos$fitted.values

```{r cyclicalRegress,fig.width=10, fig.height=4}
par(mfrow=c(1,1))
plot(x=temp_test$time, y=temp_test$Tres,
     main="Temperature (Residues)",xlab="day",ylab="T[°C]",axes=F)
axis(side=1, at=seq(1, nt, 7))
axis(side=2, at=seq(round(min(temp_test$Tres), digits=1), round(max(temp_test$Tres), digits=1), by=2))
lines(x=temp_test$time, y=temp_test$cCos, col="blue", lwd=2)
box()

#' We can combine the seasonal and cyclical components into a single model:


model.SCcos <- lm(T ~ cos(2*pi*temp$time/seasons[1]) + sin(2*pi*temp$time/seasons[1]) +
                        cos(2*pi*temp$time/periods[1]) + sin(2*pi*temp$time/periods[1]) +
                        cos(2*pi*temp$time/periods[1]*2) + sin(2*pi*temp$time/periods[1]*2),
                      temp)
summary(model.SCcos)

temp_test$SCcos <- model.SCcos$fitted.values[1:nt]
temp_eval$SCcos <- model.SCcos$fitted.values[(nt + 1):n]

```{r finalRegress,fig.width=10, fig.height=5}
par(mfrow=c(1,1))
plot(T ~ temp_test$time, temp_test,
     main="Modeling Seasonal and Cyclical Components",xlab="day",ylab="T[°C]",axes=F)
lines(temp_test$time, temp_test$sCos, col="red", lwd=2)
lines(temp_test$time, temp_test$SCcos, col="blue", lwd=2)
axis(side=1, at=seq(1, nt, 7))
axis(side=2, at=seq(round(min(temp_test$T), digits=1), round(max(temp_test$T), digits=1), by=2))
box()
legend("topleft", legend=c("seasonal","seasonal + Fourier"),
       col=c("red","blue"), lty=1,lwd=2 , cex=0.8)

#length(temp_test$Tres)
#length(model.sCos$residuals)

temp_test$Tres <- model.sCos$residuals[1:nt]
temp_test$Tres2 <- model.SCcos$residuals[1:nt] # for future use

temp_eval$Tres2 <- model.SCcos$residuals[(nt + 1):n]

```{r finalResidues,fig.width=10, fig.height=3.5}
plot(x=temp_test$time, y=temp_test$Tres, type="l",
     main="Residuals After Seasonal",xlab="day",ylab="T[°C]",axes=F)
axis(side=1, at=seq(1, nt, 7))
axis(side=2, at=seq(round(min(temp_test$Tres), digits=1), round(max(temp_test$Tres), digits=1), by=2))
box()
plot(x=temp_test$time, y=temp_test$Tres2, type="l",
     main="Residuals After Seasonal & Cyclical",xlab="day",ylab="T[°C]",axes=F)
axis(side=1, at=seq(1, nt, 7))
axis(side=2, at=seq(round(min(temp_test$Tres), digits=1), round(max(temp_test$Tres), digits=1), by=2))
box()

#'
#'---------------------------------------------------------------------------------------
#'
#' ### 3.2 Radomness Tests on Residuals After Extracting Systematic Components ###
#'
#' Even though the original temperature time series has no trend, the residuals
#' of the resultant time series with extracted both seasonal (and cyclical) components, that is: the systematic
#' components might still contain some cyclical components. To verify that the residuals are a trajectory of a
#' (stationary) stochastic process we use the, so called, "Randomness tests":
#' 
#' 1.) Durbin-Watson autocorrelation test
#' 
#' 2.) Zero ACF test
#' 
#' 3.) Signed Rank Test 
#' 
#' 4.) Spearman Rho Test
#' 
#' 5.) Turning Point Test 
#' 
#' 6.) Median Test 
#' 
#' 
#' We begin by examining the residues after extracting all the systematic components in `model.sCos`, as well as their
#' ACF:

```{r resPlot1, fig.width=10, fig.height=5}
par(mfrow=c(1,2))
plot(x=temp_test$time, y=temp_test$Tres,
     main="Temperature",xlab="day",ylab="T[°C]",axes=F, type="l")
axis(side=1, at=seq(1, 356, 7))
axis(side=2, at=seq(round(min(temp_test$Tres), digits=1), round(max(temp_test$Tres), digits=1), by=2))
box()

acf(temp_test$Tres, lag.max=365, main="Residual ACF")

#' As we can clearly see, the residues still have a hidden oscillatory component with a period of around 30 days.
#' Hence we move ahead to test their randomness via the given tests.
#' 
#' #### 3.2.1 Durbin-Watson Test ####
#' The null hypothesis of the D-W test states that the residues of a time series after a least squares regression are 
#' uncorrelated, which is put against the alternative hypothesis: that the residuals follow a 1st order autoregressive (AR)
#' process (see section 5.). 

#' We can either use the inbuilt function in the "car" package:

#install.packages("car")


car::durbinWatsonTest(temp_test$Tres)

#' or write our own function:

# only unique values in the time series will be considered
values1 <- as.numeric(temp_test$Tres)
values2 <- unique(values1)

nv <- length(values2)
den <- sum(values2 ^ 2)
( DW <- (sum((values2[2:nv] - values2[1:(nv - 1)])^2))/den )

#' The resulting DW statistic is then compared with the D-W critical values for a sample of given size and a number
#' `k` of the terms of linear regression (intercept included), in our case `k = 2` and `size = 265`, hence in a 5% confidence interval:

#' dL         dU
# 1.78900  1.80444 (n = 260)
# 1.79306  1.80792 (n = 270)

dL <- 0.5 * (1.78900 + 1.79306)
dU <- 0.5 * (1.80444 + 1.80792)

if(DW <= dL) {
  message("Alternative Hypothesis: Positive Autocorrelation!")
} else if (4 - dL < DW && DW < 4) {
  message("Alternative Hypothesis: Negative Autocorrelation!")
} else if (dU < DW && DW < 4 - dU) {
  message("Null Hypothesis: No Autocorrelation!")
}

#' We can obtain p-values using the model matrix of the `Scos` model, simulating the stochastic process
#' `Y ~ X - 1`, computing a DW-statistic for each simulation, and counting the number of `DW < simDW`:

X <- model.matrix(model.sCos)[1:nt,]

reps = 1000
sig <- var(values2)
mu <- temp_test$sCos
Y <- matrix(rnorm(nv*reps, 0, sig), nv, reps) + matrix(mu, nv, reps)
E <- residuals(lm(Y ~ X - 1))
simDW <- apply(E, 2, function(e) (sum((e[2:length(e)] - e[1:(length(e) - 1)])^2) / sum(e ^ 2)) )

( p_val <- (sum(DW > simDW)) / reps ) # for negative correlation
( p_val <- (sum(DW < simDW)) / reps ) # for positive correlation
( p_val <- 2*min(p_val, 1 - p_val) ) # for two-sided correlation

#' The "rho" parameter appearing in the test is the, so called, autocorrelation coefficient for which -1 < rho < 1
#' and the null hypothesis translates to verifying that rho = 0 for some confidence interval. 
#' 
#' We can model correlation for other lags as well. Aside from the residual vector as an argument, the `durbinWatsonTest`
#' can process the entire regression model. It also has a `max.lag` argument which tells the function which lag to test for
#' when using a model `Y ~ X - lag`. Parameter `alternative` also tells the function which type of correlation to test for:

car::durbinWatsonTest(model.sCos, max.lag=10, alternative="two-sided")
car::durbinWatsonTest(model.sCos, max.lag=10, alternative="negative")
car::durbinWatsonTest(model.sCos, max.lag=10, alternative="positive")

#' We see that DW tests have zero p-value up to lag 3. This means that we could be dealing with an AR process of order 3 or
#' higher.
#'
#' #### 3.2.2 Zero ACF Test ####
#' The test reduces to just finding all the ACF lags k with ACF(k) > 2/sqrt(n) where n is the length of the time series.

```{r acfPlot1, fig.width=11, fig.height=5}
par(mfrow=c(1,1))
ACF <- acf(temp_test$Tres, lag.max=nt, plot=F)
lags <- which(abs(acf(temp_test$Tres, lag.max=nt, main="Residual ACF")$acf[-1]) > 2/sqrt(nt))
lagValues <- numeric()
for(i in 1:length(lags))
{
  lagValues[i] <- ACF$acf[lags[i] + 1]
  segments(lags[i], 0, lags[i], lagValues[i], col= "red", lwd=2)
}
lags

#' The result shows that the residual ACF exceeds the zero-value for the lags shown above.
#'
#' #### 3.2.3 Signed Rank Test ####
#' This is a non-parametric test for the presence of a trend in the residual time series, based on the analysis 
#' of differences of sets of three consecutive terms. The null hypothesis states that the resulting statistic is
#' asymptotically normally distributed.
#' using a custom approach:

# consecutive values with zero difference have to be omitted
dist_bool <- c(T,as.logical(diff(values2)))
Tt <- c()
for (i in 1:length(dist_bool)) 
{
  if (dist_bool[i]) Tt <- c(Tt, values2[i])
}

#now the ranking of differences
ranks <- c()
for (i in 2:length(Tt)) 
{
  if (Tt[i-1] < Tt[i]) ranks[i - 1] <- 1
  else 
  {
    ranks[i - 1] <- 0
  }
}

rank_sum <- sum(ranks)
rank_mean <- (length(Tt) - 1) / 2
rank_var <- (length(Tt) + 1) / 12
test_stat <- (rank_sum - rank_mean) / sqrt(rank_var)

alpha = 0.05

if (test_stat < qnorm(alpha) || test_stat > qnorm(1 - alpha)) {
  message(paste("test_stat = ",test_stat," alpha = ", alpha,
                " ::: Alternative Hypothesis! test_stat is not normally distributed for n-->infty ! Trend detected."))
} else {
  message(paste("test_stat = ",test_stat," alpha = ", alpha,
                " ::: Null Hypothesis! test_stat is normally distributed for n-->infty ! No trend detected."))
}

#' And we can compare the results with an inbuilt function from the "randtests" package:

require(randtests)
randtests::difference.sign.test(temp_test$Tres)

#' So as it appears, for a 5% confidence interval the residues still contain a trend. Substituting 0.01 for alpha, however,
#' gives an alternative hypothesis as a result. Thus the signed-rank test verifies its null hypothesis for 0.1 > alpha > 0.05.
#' This is probably caused by cutting off a portion of the data at the end of an annual sample.
#'
#'
#' #### 3.2.4 Spearman Rho Test ####
#' Another non-parametric test for the presence of a trend using concordant/discordant pairs of values can be implemented as:

#order the unique values in an increasing order

temp_ordered <- values2[order(values2, decreasing=F)]
q <- numeric()

for (i in 1:length(values2)) 
{
  q[i] <- match(c(values2[i]), temp_ordered) 
}

# and for the spearman rho coefficient
nv <- length(values2)
sum_term <- numeric()
for (i in 1:length(values2)) 
{
  sum_term[i] <- (i - q[i])^2 
}

rho <- 1 - (6 / (nv * (nv * nv - 1))) * sum(sum_term)

test_stat <- abs(rho) * sqrt(nv - 1)

if (test_stat <= qnorm(alpha) || test_stat >= qnorm(1 - alpha)) 
{
  message(paste("test_stat = ",test_stat," alpha = ", alpha,
                " ::: Alternative Hypothesis! test_stat is not normally distributed for n-->infty ! Trend detected."))
} else {
  message(paste("test_stat = ",test_stat," alpha = ", alpha,
                " ::: Null Hypothesis! test_stat is normally distributed for n-->infty ! No trend detected."))
}

#' And using the inbuilt `cor.test`:

# our rho:
rho

cor.test(temp_test$time, temp_test$Tres, method="spearman")

#' Clearly both the Signed-Rank, and Spearman Rho tests show that the residual time series is a realisation of independent 
#' identically distributed random variables.
#'
#' #### 3.2.5 Turning Point Test ####
#' This and the following test are both tests for the presence of (previously unfiltered) periodic components in the 
#' residual time series. Similarily to the signed rank test, it examines sets of three consecutive terms, except this time
#' looking for so called 'turning points', that is: triplets of values in which the middle value is locally extremal 
#' (in that triplet). The presence and periodic regularity of turning points implies oscillatory behavior or the residues.

dist_bool <- c(T,as.logical(diff(values2)))

#removing duplicate consecutive values
Tt <- c()
for (i in 1:length(dist_bool)) 
{
  if (dist_bool[i] == TRUE)  Tt <- c(Tt, values2[i])
}

# marking consecutive triplets with value 0 for upper and lower turning points
# and 1 otherwise

Mt <- c()
for (i in 2:(length(Tt) - 1)) 
{
  if ( ((Tt[i - 1] < Tt[i]) && (Tt[i] > Tt[i + 1])) || 
       ((Tt[i - 1] > Tt[i]) && (Tt[i] < Tt[i + 1])) ) 
  {
    Mt[i - 1] <- 1
  } else {
    Mt[i - 1] <- 0
  }
}

M_sum <- sum(Mt)
M_mean <- 2 * (length(Tt) - 2) / 3
M_var <- (16 * length(Tt) - 29) / 90

test_stat <- (M_sum - M_mean) / sqrt(M_var)

if (test_stat <= qnorm(alpha) || test_stat >= qnorm(1 - alpha)) {
  message(paste("test_stat = ",test_stat," alpha = ", alpha,
                " ::: Alternative Hypothesis! Periodic components detected."))
} else {
  message(paste("test_stat = ",test_stat," alpha = ", alpha,
                " ::: Null Hypothesis! No periodic components detected."))
}

( p_value <- pnorm(test_stat) ) # positive serial correlation
( p_value <- 2*min(p_value, 1 - p_value)) # (two-sided) non-randomness
( p_value <- 1 - p_value ) # negative serial correlation

#' And via the inbuilt function from randtests:

randtests::turning.point.test(temp_test$Tres, alternative="left.sided")
randtests::turning.point.test(temp_test$Tres, alternative="two.sided")
randtests::turning.point.test(temp_test$Tres, alternative="right.sided")

#' we, apparently, obtain the same result. 
#' 
#' #### 3.2.6 Median Test ####
#' This is another non-parametric test for the presene of periodic components which essentially divides the sample values
#' into distinct groups each of which corresponds to a group of values that are unilaterally above or below the sample median,
#' meaning that if the `i`-th value and its predecessors lie above the median, and the `(i+1)`-th value lies below, for instance, 
#' the `(i+1)`-th value is a member of a new group of values which lie below the sample median.

U <- temp_test$Tres
temp_med <- median(U)
U_ <- U - temp_med
U_ <- U_[U_ != 0]
M <- U_ > 0
head(as.numeric(M), 10)

c(below=sum(diff(M)<0), above=sum(diff(M)>0))

m <- sum(M)
P <- sum(diff(M)>0) + sum(diff(M)<0) + 1

```{r medianPlot1, fig.width=11, fig.height=5}
par(mfrow=c(1,1))
plot(temp_test$Tres ~ temp_test$time, main="Median", ylab="res")
abline(h=temp_med, col="green")
for(i in 1:length(temp_test$Tres)) {
  if (U[i] != temp_med ) {
    if (U[i] > temp_med){
      segments(temp_test$time[i], temp_med, temp_test$time[i], U[i], col= "red", lwd=1)
    } else {
      segments(temp_test$time[i], temp_med, temp_test$time[i], U[i], col= "blue", lwd=1)
    }
  }
}

( ZStatistic <- (P - (m+1)) / sqrt(m*(m-1)/(2*m-1)) )


if (ZStatistic > qnorm(alpha / 2) && ZStatistic < qnorm(1 - alpha/2)) {
  message(paste("qnorm(",alpha/2,") < Zstat < qnorm(",1 - alpha/2,")"))
  message("", round(qnorm(alpha / 2), digits=2)," < ", round(ZStatistic, digits=2)," < ",round(qnorm(1 - alpha/2), digits=2),"")
  message("Null Hypothesis! Randomness!")
  if (m <= 100) cat("Note: 'above median' group count: m = ", m," <= 100")
} else {
  message(paste("qnorm(",alpha/2,") <= Zstat || Zstat >= qnorm(",1 - alpha/2,")"))
  message("", round(qnorm(alpha / 2), digits=2)," <= ", round(ZStatistic, digits=2)," >= ",round(qnorm(1 - alpha/2), digits=2),"")
  message("Alternative Hypothesis! Non-Randomness ---> Periodicity!")
  if (m <= 100) cat("Note: 'above median' group count : m = ", m," <= 100")
}

( p_value <- pnorm(ZStatistic) ) # positive serial correlation
( p_value <- 2*min(p_value, 1 - p_value)) # (two-sided) non-randomness
( p_value <- 1 - p_value ) # negative serial correlation

randtests::runs.test(U)

# library(signmedian.test) # NOTE: signmedian.test does not do the same thing as the median test above
# signmedian.test(temp_test$Tres, alternative="two.sided")
# signmedian.test(temp_test$Tres, alternative="greater")
# signmedian.test(temp_test$Tres, alternative="less")

#' The median test implies randomness of the residual time series, yet its results are irrelevant for group count
#' `m <= 100`. The result of the Durbin-Watson Test implies that the residual time series may be a trajectory of an AR(3) 
#' (auto-regressive) process, and the residual ACF, of course, shows non-zero correlation for multiple lags. The signed-rank
#' test "almost" implies that the residues do not have a trend (more specifically, for `alpha = 0.1`), and the following
#' Spearman rho test does so as well. The presence of a trend may be caused by cutting off the tail of the annual series.
#' The Turning Point test, with the Median test which  followed show that the time series also contains unfiltered periodic
#' components, even though the Median test may be inconclusive since the group count m of groups of 
#' values above median is less than 100.
#' 
#'  Therefore, from all the tests carried out in this section, we can conclude the following:
#' 
#' - The residual time series may be a trajectory of an AR(3) (auto-regressive) process
#' - There is no trend present in the residual time series
#' - And finally the unfiltered periodic components have to be accounted for via other methods
#'
#'
#' ---------------------------------------------------------------------------------------------------------------
#' 
#' ## 4. Introducing ARMA Models and the Hannan-Rissanen Procedure ## 
#'
#' In general, a stationary (only lag-dependent) stochastic process is determined by p auto-regressive and q moving-average
#' terms. The auto-regressive terms (as by the name) determine the value of the stochastic process at each time from its past
#' values from up to p steps back, and the moving-average terms on the other hand imply dependence on a random process, such as
#' the white noise, for instance, up to q time steps back. 
#' 

```{r newRegressionAcf, fig.width=11, fig.height=5.5}
par(mfrow=c(1,3))
plot(x=temp_test$time, y=temp_test$Tres, main="Residues",xlab="day",ylab="T[°C]",axes=F, type="l")
axis(side=1, at=seq(1, 365, 7))
axis(side=2, at=seq(round(min(temp_test$Tres), digits=1), round(max(temp_test$Tres), digits=1), by=2))
box()

ACF <- acf(temp_test$Tres, lag.max=365, plot=F)
n<-length(temp_test$Tres)
lags <- which(abs(acf(temp_test$Tres, lag.max=365, main="Residual ACF")$acf[-1]) > 2/sqrt(n))
lagValues <- numeric()
for(i in 1:length(lags)) {
  lagValues[i] <- ACF$acf[lags[i] + 1]
  segments(lags[i], 0, lags[i], lagValues[i], col= "red", lwd=2)
}
acf(temp_test$Tres, lag.max=365, main="Residual PACF", type="partial")
 
#' Based on the result of the Durbin-Watson test, the additional oscillations
#'  (which can still be seen in red non-zero values of the last residual ACF) should be accounted for by modelling the 
#'  residues by a suitable AR(p) model (or ARMA(p,q) in general). The process for determining the proper orders p and q 
#'  of suitable ARMA model candidates is called the Hannan-Rissanen Procedure and it consists of multiple steps.
#'
#' First, we need to figure out whether the residues need additional transformation, so that they could be properly modelled.
#' For that we test whether:
#' 
#' (a) the time series has zero mean value
#' 
#' (b) the time series is stationary
#' 

mean(temp_test$Tres)

#' So the mean value approximation is close to zero, and the stationarity can be tested using the following test:

suppressMessages(require(tseries))
adf.test(temp_test$Tres, alternative="explosive") # Augmented Dickey-Fuller

#' According to the Augmented Dickey-Fuller tests the residuals the stationarity hypothesis is rejected.
#' However, since The tests have low statistical power in that they often cannot distinguish between true unit-root processes and 
#' near unit-root processes. This is called the "near observation equivalence" problem.
#'
#' ### 4.1. The Hannan-Rissanen Procedure ###
#' 
#' First, we need to set the maximum values for lag parameters based on the significant lags in residual ACF and PACF.
#' It will be easier to do so in a combined plot:

```{r plot03, echo=T, fig.width=8, fig.height=4}
par(mfrow=c(1,1))
acf(temp_test$Tres, ylab="(P)ACF", lag.max=100, main="Comparison of Residual ACF & PACF")
tmp <- pacf(temp_test$Tres, plot=F, lag.max=100)
points(tmp$lag+0.3, tmp$acf, col="red", type = "h")
legend("topright", legend=c("ACF","PACF"), col = c("black","red"), lty=c(1,1))

#' The partial correlogram (red) shows that the most significant correlations are for lags `L = 1`, and `L = 40`.
#' The second lag, however, does not exceed the zero region nearly as much as `L = 1`. It is then reasonable to assume
#' that the most viable model will be of AR order `p = 1`. Hence:

kmax = 50; pmax = 30; qmax = 30;

#' The Yule-Walker method based on solving regression equations against an increasing basis of AR terms can be done
#' through the inbuilt `ar()` function which automatically finds the model with the lowest AIC (Akaike's Information Criterion). 
#' And by plotting the `$aic` parameter we obtain differences `AIC_min - AIC_k` for all models.

model <- list()
model$ar.yw <- ar(temp_test$Tres, order.max=kmax)

#' Plotting the differences dAIC we notice that the highest drop in `dAIC` comes after `L=1` and then a much smaller 
#' drop follows after `L=15` and `L=40` after that. We can determine the maximum order `k` from the lags where this drop
#' in `dAIC` occurs.
#' An alternative to this approach is determining the AR order `p` from the residual variances:

```{r plot05, echo=T, fig.width=9.5, fig.height=4}
par(mfrow=c(1,2))
plot(0:(length(model$ar.yw$aic)-1), xlab="p", model$ar.yw$aic, ylab="dAIC", main="Differences in AIC")
rbind(coef=model$ar.yw$ar, se=sqrt(diag(model$ar.yw$asy.var.coef))) 
tmp <- sapply(1:kmax, function(x) ar(temp_test$Tres, aic=F, order.max=x)$var.pred)
plot(tmp, xlab="p", ylab="sigma^2", main="Residual variances")

#' The plot suggests that the highest order of the AR process should be around
#' `p = 20`. The maximum AR order can also be determined via a recursive Lewinson-Durbin algorithm.
#' For the sake of saving computation time, however, we choose the following maximum order parameters:
#'    
kmax = 20; pmax = 10; qmax = 10;

#LongAR method implementation:
LongAR = function(tser, k, p, q) {
  if (!is.ts(tser)) {
    tser = ts(tser)
  }
  
  tt <- data.frame(x=tser)
  tt$z <- ar.ols(tt$x, aic=F, order.max=k)$resid
  tt$z[is.na(tt$z)] <- 0
  
  # suppressMessages(library(dynlm))
  
  if (p > 0 & q > 0)  outmodel <- dynlm(x ~ L(x, 1:p) + L(z, 1:q), tt)
  else if (p > 0 & q == 0) outmodel <- dynlm(x ~ L(x, 1:p), tt)
  else if (q > 0 & p == 0) outmodel <- dynlm(x ~ L(z, 1:q), tt)
  else warning("LongAR::error! invalid order p,q")
  
  outmodel
}

# a method for returning model's AIC (Akaike's Information Criterion) and BIC (Bayesian Information Criterion)
# with a changeable parameter
InfCrit = function(ml_model, type="AIC") {
  if (type=="AIC") {
    AIC(ml_model)
  } else if (type=="BIC") {
    BIC(ml_model)
  } else {
    warning("Invalid type (AIC or BIC)", call=F)
  }
}


#determine the model's AIC using an inbuilt function ar()

AICs <- ar(temp_test$Tres, order.max=kmax)$aic

# find minimal AIC for pmax <= p <= kmax
minAIC <- min(AICs[pmax:kmax])
for (k in pmax:kmax) {
  if (AICs[k] == minAIC) {
    korder <- k
    break
  }
}

( k <- korder )

models.arma <- list() 

suppressMessages(require(dynlm))
models.arma[[paste(1,0,sep=",")]] <- LongAR(temp_test$Tres, k, 1, 0)
models.arma[[paste(0,1,sep=",")]] <- LongAR(temp_test$Tres, k, 0, 1)

for(p in 1:pmax) {  
  for(q in 1:qmax) {
    models.arma[[paste(p,q,sep=",")]] <- LongAR(temp_test$Tres, k, p, q)
  }
}

bic <- cbind(
  BIC = sapply(models.arma, function(x) InfCrit(x, type="BIC")),
  AIC = sapply(models.arma, function(x) InfCrit(x, type="AIC"))
)
bic <- bic[order(bic[,"BIC"]),]

head(bic, n=5)

#' Now we have 5 ARMA models with the lowest BIC with orders "p,q".

bestArma <- head(bic, n=5)
orders <- rownames(bestArma)
topModels <- list() # here I put the top models marked by their orders "p,q" as a key

for (j in 1:length(orders)) {
  key <- orders[[j]]
  topModels[[key]] <- models.arma[[key]]
}

( bestArma <- cbind(bestArma, topModels) ) #checking if the list contains models
#and then accessing them in the following way:
bestArma[[1,3]] #first in the list and the third column for the model itself


#' ### 4.2. Adjusting the Regression Coefficients Using the Maximum Likelihood Estimate ###
#' 
#' The Maximum Likelihood Estimate (MLE) optimizes the, so called, likelihood-function with respect to
#' the regression coefficients. One can find estimates using the Residual Square Sum (RSS). 
#' We will optimize the model parameters via an inbuilt function arima, specifying parameter `method = "ML"`
#' and also using the former `dynlm`-type model and directly calculating the RSS.

MaxLikelihoodOptim = function (tmpmodel, use_arima=TRUE) {
  params <- names(tmpmodel$coefficients)
  ( p <- sum(grepl("x",params)) )
  ( q <- sum(grepl("z",params)) )
  
  pars0 <- tmpmodel$coefficients
  
  #inbuilt function
  if (use_arima) {
    coefs <- arima(temp_test$Tres, order = c(p, 0, q), 
                   init=c(pars0[-1], pars0[1]), method = "ML", transform.pars = F)$coef
    
    # making the intercept coefficient first in the list
    intersect <- coefs[length(coefs)]
    coefs <- coefs[1:(length(coefs)-1)]
    coefs <- c(intersect, coefs)
    
    return(coefs)
  } else {
    #other version
    qmp <- max(0,q - p)
    x <- as.numeric(temp_test$Tres)
    ntest <- length(x)
    
    #residual square sum
    RSS = function(pars) {
      z <- rep(0,length(x) + qmp)
      for(i in (p + 1):ntest) {
        xz <- c(1, x[i:(i - p)], z[(i:(i - q)) + qmp])
        z[i + qmp] <- x[i] - c(append(head(pars, n = p + 1), 0, after = 1), 0, tail(pars, n = q)) %*% xz
      }
      z <- z[-(1:(p + qmp))]
      sum(z^2)
    }
    
    optim(pars0, RSS)$par
  }
}

#model comparison function
CompareModels = function(models) {
  orders <- rownames(models)
  
  output <- list()
  
  for (i in 1:length(orders)) {
    HannanRissanen_Result <- models[[i,3]]$coefficients
    MLE_using_ARIMA <- MaxLikelihoodOptim(models[[i,3]])
    MLE_not_ARIMA <- MaxLikelihoodOptim(models[[i,3]], use_arima=F)
    output[[ orders[[i]] ]] <- cbind(HannanRissanen_Result, MLE_using_ARIMA, MLE_not_ARIMA)
  }
  
  output
}

( comparison <- CompareModels(bestArma) )

#changing model coefficients
OptimizeModels = function (models, use_arima=F) {
  orders <- rownames(models)
  
  for (i in 1:length(orders)) {
    if (use_arima) {
      result_coeffs <- MaxLikelihoodOptim(models[[i,3]])
    } else {
      result_coeffs <- MaxLikelihoodOptim(models[[i,3]], use_arima=F)
    }
    
    models[[i,3]]$coefficients <- result_coeffs
    
  }
  models
}

newBestArma <- OptimizeModels(bestArma, use_arima=F)

# checking if the coefficients were changed
bestArma[[5,3]]$coefficients
newBestArma[[5,3]]$coefficients

#' ### 4.3. Plotting the Resulting Models With their Residuals  ###

```{r plot11csaa, echo=T, include=T, fig.width=9, fig.height=6}
par(mfrow=c(2,2))

for (i in 1:5) {
  plot(temp_test$Tres, type="p", main=paste("ARMA(",orders[[i]],")"), ylab="T")
  lines(newBestArma[[i,3]]$fitted.values, col="blue")
  plot(newBestArma[[i,3]]$residuals, main=paste("ARMA(",orders[[i]],") residuals"), ylab="res")
}

#' ### 4.4. Plotting the Resulting Models In the Original Time Series ###

```{r plot1212, echo=T, fig.width=8, fig.height=4}
par(mfrow=c(1,1))
n <- length(temp_test$Tres)
for(i in 1:5) {
  o <- unlist(strsplit(orders[[i]],",")); p <- as.numeric(o[[1]]); q <- as.numeric(o[[2]])
  m <- length(newBestArma[[i,3]]$fitted.values)
  plot(temp$T[1:m], type="p", main=paste("Systematic lm + ARMA(",orders[[i]],")"), xlab="day",ylab="T")
  
  syst_fit <- model.sCos$fitted.values[1:m]
  lines(syst_fit, col="red", lwd=2)
  lines(newBestArma[[i,3]]$fitted.values + syst_fit, col="blue", lwd=2)
  legend("topleft", legend=c("(1): Seasonal", paste("(2): (1) + ARMA(",orders[[i]],")")),
         col=c("red","blue"), lty=1,lwd=2 , cex=0.75)
}

#'
#' The Hannan-Rissanen procedure with maximum orders: `kmax = 20`, `pmax = 10`, and `qmax = 10`, determined the optimal ARMA models
#' to be: 
row.names(newBestArma)
#' Combining the systematic components, i.e.: the seasonal and cyclical 
#' components, with the linear regression of a given ARMA model we have found 5 most accurate linear models for the test part
#' of the original time series. The evaluation part will be used in the following section where we will test the given models
#' and carry out predictions.
#' 
#' ## 5. ARMA Model Diagnostics and Predictions ##
#' 
#' Prior to carrying out (single-step and multiple step) predictions, we need to test the accuracy of the given models. For that we use the 
#' following series of diagnostic tests:
#' 
#' - Zero Autocorrelation
#' 
#' - Normality of Residues
#' 
#' - Conditional Heteroskedasticity
#' 
#' ### 5.1. Zero Autocorrelation ###
#' 
#' The Autocorrelation function (ACF) shows a reasonably clear image of the dependency of individual random variables with
#' respect to time lag. This means that if the ACF has non-zero values for some lags greater than zero, then it is likely
#' that the residues after removing the ARMA regression's fitted values still contain some systematic or ARMA component.
#' To find out we test for a null hypothesis: `ACF(k) = 0` for all `k > 0`, against an alternative `ACF(k) != 0` for some `k > 0`.

```{r plot1zaa, echo=T, fig.width=9, fig.height=4}
ZeroACF = function (x, zeroValue) {
  x = ifelse (abs(x) < zeroValue, 0, x)
}

NonZeroValues = function (series, zero) {
  indices <- which(series != 0)
  nonzeros <- list()
  
  nonzeros[["zero val"]] <- zero
  
  for (i in 2:length(indices)) {
    nonzeros[[paste("k",indices[i], sep = '_')]] <- series[indices[i]]
  }
  
  data.frame(head(nonzeros))
}

ACFs <- list()
Nonzeros <- list()

par(mfrow=c(1,2))
for (i in 1:5) {
  acf(newBestArma[[i,3]]$residuals, ylab="ACF", lag.max=nt, main=paste("ARMA(",orders[[i]],") res"))
  
  ACFs[[i]] <- acf(newBestArma[[i,3]]$residuals, ylab="ACF", lag.max=nt, plot=F)$acf
  zero <- 2 / sqrt(nt)
  ACFs[[i]] <- sapply(ACFs[[i]], function(x) ZeroACF(x, zero))
  lines(ACFs[[i]], col="red", lwd=2)
  
  Nonzeros[[i]] <- NonZeroValues(ACFs[[i]], zero)
}



#' Using the zero-value approximation: `2/sqrt(nt)`, where `nt` is the length of the residual time series,
#' we notice that there seems to be a resilient 23-day lag which was not removed (even after including the significant
#' residual lags such as 14, and 23 in the model.sCos in section 3). The significance of the 23-day lag, however,
#' is not nearly as high as it would have otherwise been (with additional significant lags) if additional lags had not been 
#' introduced prior to carrying out the Hannan-Rissanen procedure. Thus we can assume that all models, having
#' no other significant lags than the 23-day lag, describe the behavior of the test part of the residual time series 
#' reasonably well.
#' In the following list:

Nonzeros

#' we see the exact lags for which the ACF exceeds the zero-value.
#' 
#' ### 5.2. Normality of Residues ###
#' 
#' If a given ARMA model describes the time series well enough, the residues should behave as if they were generated
#' by a stochastic process with values at individual times `t` being the realizations of independent identically distributed
#' random variables. For sufficiently large number of time steps, the residues are assumed to be the results of a normally 
#' distributed random variable. 
#' 
#' To test that, we use the, so called, Jarque-Bera test:

```{r JBTest, echo=T}
JBToutput <- list()

for (i in 1:5) {
  JBToutput[[i]] <- tseries::jarque.bera.test(newBestArma[[i,3]]$residuals)
  
  pValue <- as.numeric(JBToutput[[i]]$p.value)
  
  JBToutput1 <- c("1,0",JBToutput[[i]]$statistic, JBToutput[[i]]$parameter, pValue)
  names(JBToutput1) <- c("p,q","statistic", "DOF", "p.value")
  data.frame(head(JBToutput1))
  if (pValue > 0.05) {
    message(paste("ARMA(",orders[[i]],") ::: p-value = ", pValue,", Normality hypothesis rejected"))
  } else { message(paste("ARMA(",orders[[i]],") ::: p-value = ", pValue,", Normality hypothesis accepted")) }
}


#' As it appears, the normality of the residues of all top 5 ARMA models is accepted. 
#'
#' ### 5.3. Conditional Heteroskedasticity ###
#' 
#' The term "heteroskedasticity" of a time series describes the tendency of its variance to change with time.
#' A time series whose variance (up to a given time) remains the same are called "homoskedastic".
#' The following Breusch-Pagan test assesses the homoskedasticity of a given model:

```{r Skedast, echo=T}
alpha = 0.05
SkedtestResults <- list()
for (i in 1:5) {
  (SkedtestResults[[i]] <- lmtest::bptest(newBestArma[[i,3]]) )
  
  pValue <- as.numeric(SkedtestResults[[i]]$p.value)
  
  if (pValue > alpha) {
    message(paste("ARMA(",orders[[i]],") ::: p-value = ", pValue,", Null hypothesis rejected: model is heteroskedastic"))
  } else {
    message(paste("ARMA(",orders[[i]],") ::: p-value = ", pValue,", Null hypothesis accepted: model is homoskedastic"))
  }
} 

#' Thus on the 5% significance level we accept the null hypothesis for all models,
#' i.e.: all models are heteroskedastic.
#'
#'-----------------------------------------------------------------------------------------------------------------
#'
#' ### 5.4. Forecasting and Single-Step Predictions ###
#' 
#' Now we make use of the evaluation part of the time series which we have separated from the original series earlier
#' in section 4. Since the `dynlm` model we used earlier does not support `predict()` function from the `forecast` library
#' we use equivalent `ARIMA(p,0,q)` models instead of `ARMA(p,q)`, which have an inbuilt function `arima()`, the result of which
#' is an object that can be processed by the `predict()` call.

( neval <- length(temp_eval$Tres) )

```{r predictplot1, echo=TRUE, fig.width=9, fig.height=4}
par(mfrow=c(1,1))
library(forecast)
eval.predictions <- list()
merrors <- list()


for (i in 1:5) {
  #extract model orders from the top models from Hannan-Rissanen:
  o <- unlist(strsplit(orders[[i]],",")); p <- as.numeric(o[[1]]); q <- as.numeric(o[[2]])
  #estimate a suitable model on the test part of the time series
  model_estim <- arima(temp_test$Tres, order=c(p,0,q))
  #use the estimated model as an argument, fitting the same model to the whole time series without re-estimating any parameters
  eval.predictions[[i]] <- Arima(model.sCos$residuals, order=c(p,0,q), model=model_estim)
  #and take the evaluation part for examination
  eval.predictions[[i]] <- tail(fitted(eval.predictions[[i]]), neval)
  merrors[[i]] <- mean((temp_eval$Tres - eval.predictions[[i]])^2)

  plot(x=temp_eval$time, y=temp_eval$T, 
       main=paste("Predictions ARMA(",p,",",q,"), RMSE = ", round(sqrt(merrors[[i]]), 4)),
       xlab="day", ylab="[°C]", type="l", lwd=1.5)
  lines(x=temp_eval$time, y=(eval.predictions[[i]] + temp_eval$sCos), col="blue", lwd=2)
  legend("bottomleft", legend=c("prediction", "actual values"),
         col=c("blue","black"), lty=1,lwd=2 , cex=0.75)
}

#' We have also computed Root Mean Square Errors for each model, with the smallest one being:

minE <- which.min(merrors)

#' ### 5.5. Comparison of Predictive Abilities Using the Diebold-Mariano Test ###
#' 
#' Taking two different ARMA models, the null hypothesis of the test states that both models have equal predictive
#' abilities. We will be comparing 5 different models against each other, which will essentially amount to testing
#' the predictive abilities of 10 distinct pairs of models. Practically, the test results can be expressed in a 5x5
#' "DieboldMarianoMatrix" with values 0 (if model i has the same predictive abilities as model j), 1 (if model i has
#' better predictive abilities than model j), and -1 (the other way around):

DieboldMarianoMatrix <- matrix(0, ncol=5, nrow=5)
pvalueMatrix <- matrix(1, ncol=5, nrow=5)
for(i in 1:length(eval.predictions)) {
  for(j in 1:length(eval.predictions)) {
    if(i==j) next
    test <- forecast::dm.test(temp_eval$Tres - eval.predictions[[i]], temp_eval$Tres - eval.predictions[[j]])
    DieboldMarianoMatrix[i,j] <- (test$p.value < alpha) * sign(test$statistic) 
    pvalueMatrix[i,j] <- round(test$p.value, 3)
  }
}
DieboldMarianoMatrix

pvalueMatrix

#' As the results of the Diebold-Mariano test suggest, models 3 and 5 have better predictive abilities than all the remaining models
#' on a 5% significance level.
#' 
#' A quantitative characteristic of predictive ability is the Root Mean Square Error (RMSE).
#' Here we can see how the predictions match the observed values from the evaluation part of the time series:

```{r predictPlot2, echo=T, fig.width=9, fig.height=6}
par(mfrow=c(1,1));
o <- unlist(strsplit(orders[[minE]],",")); p <- as.numeric(o[[1]]); q <- as.numeric(o[[2]])
n <- length(model.sCos$fitted.values)
fullpredict <- as.numeric(eval.predictions[[minE]] + temp_eval$sCos)
plot(x=temp$time, y=temp$T, main=paste("Prediction + Systematic , ARMA(", orders[[minE]],"):"), ylab="T",xlab="day")
lines(x=temp_eval$time, y=fullpredict, col="blue", lwd=2)

k = max(p, q)
syst_fit <- model.sCos$fitted.values[(k + 1):nt]
test_model <- newBestArma[[minE,3]]$fitted.values + syst_fit

lines(x=temp_test$time[(k + 1):nt], y=test_model, col="red",lwd=2)
lines(x=temp$time[nt:(nt + 1)], y=c(test_model[nt - 1], fullpredict[1]), col="red", lwd=2) #line connecting first and second value set
legend("topleft", legend=c(paste("Model ARMA(",orders[[minE]],")"), paste("Model ARMA(",orders[[minE]],") prediction")),
       col=c("red","blue"), lty=1,lwd=2 , cex=0.75)

#'
#' The model with the smallest RMSE is ARMA(1,0). Except for the outlying values, the predictive capabilities of the model
#' can also be verified visually from the plot above.
#' Now recall that we have fully neglected the effect of cyclical components. The residues after removing cyclical components as well
#' have been stored in `temp_test$Tres2`. We can run the whole H-R procedure again for these residues to see if an ARMA type model can
#' be found which would have better predictive abilities:
#' 

AICs <- ar(temp_test$Tres2, order.max=kmax)$aic

minAIC <- min(AICs[pmax:kmax])
for (k in pmax:kmax) {
  if (AICs[k] == minAIC) {
    korder <- k
    break
  }
}

( k <- korder )

models.arma2 <- list() 

suppressMessages(require(dynlm))
models.arma2[[paste(1,0,sep=",")]] <- LongAR(temp_test$Tres2, k, 1, 0)
models.arma2[[paste(0,1,sep=",")]] <- LongAR(temp_test$Tres2, k, 0, 1)

for(p in 1:pmax) {  
  for(q in 1:qmax) {
    models.arma2[[paste(p,q,sep=",")]] <- LongAR(temp_test$Tres2, k, p, q)
  }
}

bic2 <- cbind(
  BIC = sapply(models.arma2, function(x) InfCrit(x, type="BIC")),
  AIC = sapply(models.arma2, function(x) InfCrit(x, type="AIC"))
)
bic2 <- bic2[order(bic2[,"BIC"]),]

head(bic2, n=5)
head(bic, n=5)

#' Using `temp_test$Tres2` does not change the ARMA model selection at all. Hence we will only need to compare the resulting 
#' error when cyclical components are included (in model `model.SCcos`).

bestArma2 <- head(bic2, n=5)
orders2 <- rownames(bestArma2)
topModels2 <- list()

for (j in 1:length(orders2)) {
  key <- orders2[[j]]
  topModels2[[key]] <- models.arma2[[key]]
}

( bestArma2 <- cbind(bestArma2, topModels2) ) 

```{r acf2plot, echo=T, fig.width=9, fig.height=3.5}
acf(bestArma2[[1,3]]$residuals, ylab="ACF", lag.max=nt, main=paste("ARMA(",orders2[[1]],") res"))

ACFs <- acf(bestArma2[[1,3]]$residuals, ylab="ACF", lag.max=nt, plot=F)$acf
zero <- 2 / sqrt(nt)
ACFs <- ZeroACF(ACFs, zero)
lines(ACFs, col="red", lwd=2)

( Nonzero <- NonZeroValues(ACFs, zero) )

#' the significant lags in the acf of the first ARMA model's residuals are still present.
#' 

eval.predictions2 <- list()
merrors2 <- list()

```{r predictPlot3, echo=T, fig.width=9, fig.height=4}
for (i in 1:5) {
  o <- unlist(strsplit(orders[[i]],",")); p <- as.numeric(o[[1]]); q <- as.numeric(o[[2]])
  model_estim2 <- arima(temp_test$Tres2, order=c(p,0,q))
  eval.predictions2[[i]] <- Arima(model.SCcos$residuals, order=c(p,0,q), model=model_estim2)
  eval.predictions2[[i]] <- tail(fitted(eval.predictions2[[i]]), neval)
  
  merrors2[[i]] <- mean((temp_eval$Tres2 - eval.predictions2[[i]])^2)
  
  plot(x=temp_eval$time, y=temp_eval$T, 
       main=paste("Predictions ARMA(",p,",",q,"), RMSE = (", round(sqrt(merrors[[i]]), 4), ",",round(sqrt(merrors2[[i]]), 4)," )"),
       xlab="day", ylab="[°C]", type="l", lwd=1.5)
  lines(x=temp_eval$time, y=(eval.predictions2[[i]] + temp_eval$SCcos), col="red", lwd=2)
  lines(x=temp_eval$time, y=(eval.predictions[[i]] + temp_eval$sCos), col="blue", lwd=2)
  legend("bottomleft", legend=c("prediction1", "prediction2", "actual values"),
         col=c("blue","red","black"), lty=1,lwd=2 , cex=0.75)
}

errorCmp <- data.frame(
  RMSE1=round(sqrt(as.numeric(merrors)), 4),
  RMSE2=round(sqrt(as.numeric(merrors2)), 4)
)

row.names(errorCmp) <- orders

errorCmp

#' The results show that even though adding cyclical components decreases the `RMSE` of the model predictions, it does not
#' change their orders. To see if the new models have significantly better predictive abilities, we can perform a Diebold-Mariano
#' test on a twice as large set of models:

DieboldMarianoMatrix <- matrix(0, ncol=10, nrow=10)
pvalueMatrix <- matrix(1, ncol=10, nrow=10)
for(i in 1:10) {
  for(j in 1:10) {
    if(i==j) next
    if (i <= 5 && j <= 5) {
      test <- forecast::dm.test(temp_eval$Tres - eval.predictions[[i]], temp_eval$Tres - eval.predictions[[j]])
    } else if (i <= 5 && j > 5) {
      test <- forecast::dm.test(temp_eval$Tres - eval.predictions[[i]], temp_eval$Tres2 - eval.predictions2[[j - 5]])
    } else if (i > 5 && j <= 5) {
      test <- forecast::dm.test(temp_eval$Tres2 - eval.predictions2[[i - 5]], temp_eval$Tres - eval.predictions[[j]])
    } else {
      test <- forecast::dm.test(temp_eval$Tres2 - eval.predictions2[[i - 5]], temp_eval$Tres2 - eval.predictions2[[j - 5]])
    }
    DieboldMarianoMatrix[i,j] <- (test$p.value < alpha) * sign(test$statistic) 
    pvalueMatrix[i, j] <- round(test$p.value, 3)
  }
}
DieboldMarianoMatrix

pvalueMatrix

#' Which only shows that there is no significant difference between the models since the `DieboldMarianoMatrix` values
#' for rows 6 to 10 and columns 1 to 5 suggest none of the models from the second set are significantly better than the ones in the first set.
#' Thus we will consider only ARMA + seasonal systematic component models.
#' 
#' ### 5.6 Multiple-Step Predictions ###
#' 

hmax = 5; # prediction horizon

```{r multipredictPlot, echo=T, fig.width=10, fig.height=5}
par(mfrow=c(1,1));

for (i in 1:5) {
  # i <- 1
  x_predictions <- list()
  t_predictions <- list()
  
  o <- unlist(strsplit(orders[[i]],",")); p <- as.numeric(o[[1]]); q <- as.numeric(o[[2]])
  x <- as.numeric(model.sCos$residuals)
  
  phi0 = newBestArma[[i, 3]]$coefficients[1]
  phi = newBestArma[[i, 3]]$coefficients[2:(p + 1)] # AR coeffs
  theta = newBestArma[[i, 3]]$coefficients[(p + 2):(1 + p + q)] # MA coeffs
  
  for (h in 1:hmax) {
    # h <- 1
    stepmax <- max(p, q, h)
    t0 <- nt + 1
    x_preds <- numeric() # h-step prediction estimates
    
    for (t in t0:n) {
      # t <- t0
      for (j in 1:h) {
        # j <- 1
        if (q > 0) {
          zz <- (x[(t + j - q):(t + j - 1)] - ifelse(j < 2, phi %*% x[(t + j - p - 1):(t + j - 2)] + phi0, x_preds[t - t0 + j - 1]))
          zt <- sapply(1:q, function(qq) ifelse(j - qq > 0, 0, zz[qq]))
        } else {
          zt <- matrix(0)
        }
        if (p > 0) {
          xt <- as.matrix(sapply(1:p, function(pp) ifelse(j - pp > 0, x_preds[t - t0 + j - pp], x[t + j - pp] - phi0))  )
        } else {
          xt <- matrix(0)
        }
        
        if (p == 0 && q > 0) {
          xp <- as.numeric(phi0 + theta %*% zt)
        } else if (p > 0 && q == 0) {
          xp <- as.numeric(phi0 + phi %*% xt)
        } else if (p > 0 && q > 0) {
          xp <- as.numeric(
              phi0 + #  ARMA intercept
              phi %*% xt + # AR part
              theta %*% zt # MA part
          )          
        }

        x_preds[t - t0 + j] <- xp
      }
    }
    
    m <- length(x_preds)
    ne <- length(temp_eval$time)
    
    x_preds <- as.matrix(na.omit(data.frame(x_preds[1:ne] + temp_eval$sCos ) ) )
    np <- length(x_preds)
    
    #if (ne > np) {
    #  message(paste( "ne > np : ", ne ," > ", np, ", h = ", h,", ARMA(", p, ",", q,")"));
    #}
    
    x_predictions[[h]] <- x_preds
    t_predictions[[h]] <- temp_eval$time[(ne - np + 1):ne]
  }
  
  plot(x=temp_eval$time, y=(temp_eval$Tres + temp_eval$sCos), type="b", lwd=1.5, main=paste("ARMA(",p,",",q,") Multiple-Step predictions"),
       xlab="day", ylab="res(T)[°C]")
  for (h in 1:hmax) {
    lines(x=t_predictions[[h]], y=x_predictions[[h]], col=rgb(0, (1 - h / hmax), h / hmax), lwd=2)
  }
  legend("bottomleft", legend=sapply(1:hmax, function(h) paste("h = ", h)),
         col=sapply(1:hmax, function(h) rgb(0, (1 - h / hmax), h / hmax)), lty=1,lwd=2 , cex=0.75)
  
}



#'
#' --------------------------------------------------------------------------------------------------------------
#' 
#' ## 6. Non-Stationary Models ##
#' 
#' ARMA type models can be used when we safely assume stationarity of the time series. Namely, when the series has a constant variance 
#' its mean values and covariances only depend on the lag, not on time. Although ARMA models seem to sufficiently describe the given 
#' time series, we may want to try out the so called Integrated processes, and conditionally heteroskedastic processes.
#' 
#' ### 6.1. Diagnostic Tests ###
#' 
#' First we will test whether the original time series after extracting the seasonal component has the properties of a non-stationary process.
#' The most common approach involves testing for unit-roots of the process' characteristic polynomial. We begin by testing the residuals 
#' after filtering seasonal components via the Dickey-Fuller (DF) test.
#' 
#' #### 6.1.1 Dickey-Fuller Test (DF) ####
#' 
#' The null hypothesis of the DF test states: `x_t = x_{t-1} + e_t`, that is, the time series is a realization of a random walk process.
#' It is carried out by setting `dx_t = x_t - x_{t-1} = rho * x_{t-1} + e_t` where `e_t` is an i.i.d process. It suffices to show that the `rho`
#' parameter of such regression equals zero. In general the regression may assume `dx_t = mu + delta * t + rho * x_{t-1} + e_t` with additional
#' parameters `mu` (mean value) and `delta` (deterministic drift).

x <- temp_test$Tres
x <- as.numeric(x)
dx <- diff(x) # first differences
m <- length(dx)

xt1 <- x[1:m]
tt <- 1:m
reg <- lm(dx ~ xt1 + 1 + tt)

( reg.sum <- summary(reg) )

#' And after the regression of the difference time series we proceed to compute the `DF`-statistic from the `rho` parameter

rho <- reg.sum$coefficients[2,1]
sigma_rho <- reg.sum$coefficients[2,2]
( DFstat <- rho / sigma_rho )

#' The value has to be compared with the critical values for the `DF`-statistic, whose distribution is estimated via simulations
#' for particular sample sizes:

critVals <- -c(3.60, 3.50, 3.45, 3.43, 3.42, 3.41) # for alpha = 0.05
sampleSizes <- c(25, 50, 100, 250, 500, 100000) # sample size

( critVal <- approx(sampleSizes, critVals, m, rule=2)$y )

if (DFstat < critVal) {
  message(paste(" DF = ", round(DFstat, digits=4), " < ", critVal ,", rho = ", round(rho, digits=4),
                " ::: Null Hypothesis accepted. Possible stochastic trend."))
} else {
  message(paste(" DF = ", round(DFstat, digits=4), " > ", critVal, ", rho = ", round(rho, digits=4),
                " ::: Null Hypothesis rejected. Stochastic trend insignificant."))
}

#' The presence of a stochastic trend for process `dx_t = mu + delta * t + rho * x_{t-1} + e_t` is accepted on a 5% significance
#' level. This means that the time series may be modelled as an integrated process.
#' We can compare the results with the inbuilt `adf.test` function from `tseries` package:

require(tseries)
adf.test(x, k = 0) # k = 0 corresponds to lag = 1

#' #### 6.1.2 Augmented Dickey-Fuller Test (ADF) ####
#' 
#' The same procedure can be generalized for lags `k > 0` as well. For that we can just use the method from `tseries` package:

DFi <- list()
for (i in 1:5) DFi[[i]] <- suppressWarnings(adf.test(temp_test$Tres, k = i)$statistic)
DFi <- data.frame(matrix(unlist(DFi), nrow=length(DFi), byrow=T))
names(DFi) <- c("DF")
DFi

#' The ADF method sets a default lag value to 

( k <- trunc((length(x)-1)^(1/3))  )

adf.test(x)

#' #### 6.1.3 Phillips-Perron Test ####
#' 
#' Similarily to the DF test, the PP-test is, again, a unit-root test which relies on regressing `dX ~ 1 + t + X_{t - 1}`. 
#' While the ADF test adresses the issue of the time series being potentially generated by a process of higher orders by using 
#' lags `k = 1, 2, ...` and expanding the regression equation, PP-test makes a non-parametric correction to the t-statistic.
#' 

pp.test(x)

#' This makes the PP-test robust with respect to unspecified autocorrelation and heteroskedasticity.
#' 
#' #### 6.1.4 KPSS Test ####
#' 
#' Published by Kwiatkovski, Philips, Schmidt, and Shin (1992), another test approaches the stationarity-vs-non-stationarity issue
#' from an opposite point of view. Using KPSS we can test for stationarity hypothesis (about a deterministic trend). Contrary to the 
#' previously used tests, the presence of a unit root is the alternative. It can be used in cases when unit roots hypothesis cannot be
#' rejected and the time series coult be non-stationary (or we have a sample of insufficient length).
#' 
#' The test relies on testing for zero variance of an i.i.d process present in the stochastic trend part in `x_t ~ ST_t + eps_t`, or
#' alternatively with deterministic trend: `x_t ~ ST_t + delta * t + eps_t`. 

kpss.test(x, null="Level")
kpss.test(x, null="Trend")
