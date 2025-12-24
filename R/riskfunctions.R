## necessary data:

## Saunders et al linear parameters
## risk per one unit increase in BMI was 14.8% (95%CI: 13.3-16.3)
t <- log(1 - 0.148) # risk function parameter
1 - exp(t) # risk increase with 1 unit decrease


## fits from bilinear model
C <- fread(here("rawdata/general_population_piecewise_parameters.csv"))
D <- fread(here("rawdata/general_population_vcov_matrix.csv"))
## 18.0% (95%CI: 16.4-19.6) for BMI<25.0kg/m2 and 6.9% (95%CI: 4.6-9.2) for BMI>=25.0kg/m2 in
exp(C$Value[4:5]) # corresponds to above
mut <- C$Value[4:5]
t1 <- mut[1]
t2 <- mut[2]
V <- as.matrix(D[, .(slope_below_breakpoint, slope_change_above_breakpoint)])
## slope_change_above_breakpoint
## slope above = slope below + change
## a = b + c
## cov(b,a) = cov(b,b+c) = cov(b,c) + var(b)
## cov(a,a) = cov(b+c,b+c) = var(b) + var(c) + 2*cov(b,c)
W <- V
W[2, 2] <- V[2, 2] + # var(c)
  V[1, 1] + # var(b)
  2 * V[1, 2] # cov(b,c)
W[1, 2] <- W[2, 1] <- V[1, 2] + V[1, 1]
## 18.0% (95%CI: 16.4-19.6) for BMI<25.0kg/m2 and 6.9% (95%CI: 4.6-9.2) for BMI>=25.0kg/m2 in
## check -- OK: matches above
S <- exp(mvrnorm(n = 1e4, mu = mut, Sigma = W))
1e2 * (1 - colMeans(S))
1e2 * quantile(1 - S[, 1], c(0.025, 1 - 0.025))
1e2 * quantile(1 - S[, 2], c(0.025, 1 - 0.025))


## risk ratio calculator
bmirefpop <- data.table(k = 25.0, theta = 1.0) #for testing
RRfun <- function(k, theta, t) {
  (1 - t * bmirefpop$theta)^bmirefpop$k / (1 - t * theta)^k
}

## check
K <- 1e5
bmi1 <- rgamma(K, shape = 24, scale = 0.7)
bmi0 <- rgamma(K, shape = bmirefpop$k, scale = bmirefpop$theta)
## E_1[exp(t*X)] / E_0[exp(t*X)] with a shift for numerics
mean(exp(t * (bmi1 - 30))) / mean(exp(t * (bmi0 - 30)))
RRfun(24, 0.7, t) # OK

## bilinear version
## see: https://search.r-project.org/CRAN/refmans/expint/html/gammainc.html
## Γ(a,x)=Γ(a)(1−P(a,x))
## https://search.r-project.org/R/refmans/stats/html/GammaDist.html
RRfunBL0 <- function(k, theta, t1, t2) {
  x0 <- 25
  a1 <- -t1
  a2 <- -t2
  ## stuff x Γ(k,x0*(a2+1/theta))/Γ(k) =
  ## stuff x (1-P(k,x0*(a2+1/theta))) = stuff x P(k,x0*(a2+1/theta),lower=FALSE)
  lans2 <- a2 * x0 +
    pgamma(k, x0 * (a2 + 1 / theta), log = TRUE, lower = TRUE) -
    k * log(1 + a2 * theta)
  lans1 <- a1 * x0 +
    pgamma(k, x0 * (a1 + 1 / theta), log = TRUE, lower = FALSE) -
    k * log(1 + a1 * theta)
  exp(lans1) + exp(lans2)
}
## NOTE works in tests

RRfunBL <- function(k, theta, t1, t2){
 RRfunBL0(k, theta, t1, t2) /
   RRfunBL0(bmirefpop$k, bmirefpop$theta, t1, t2)
}

## check
BL <- function(x, t1, t2) {
  ans <- (x - 25)
  less <- ans < 0
  ans[less] <- t1 * ans[less]
  ans[!less] <- t2 * ans[!less]
  ans
}

## xx <- seq(from = 10, 40, by = 0.1)
## plot(xx, BL(xx, t1, t2), type = "l")
## plot(xx, BL(xx, t1, t2 - 0.1), type = "l")

mean(exp(BL(bmi1, t1, t1))) / mean(exp(BL(bmi0, t1, t1)))
## E_1[exp(t*X)] / E_0[exp(t*X)] with a shift for numerics
mean(exp(t1 * (bmi1 - 30))) / mean(exp(t1 * (bmi0 - 30)))
RRfunBL0(24, 0.7, t1, t1) / RRfunBL0(bmirefpop$k, bmirefpop$theta, t1, t1)
RRfunBL(24, 0.7, t1, t1)

mean(exp(BL(bmi1, t1, t2))) / mean(exp(BL(bmi0, t1, t2)))
RRfunBL(24, 0.7, t1, t2) ## OK

## correct magnitude of change?
mean(bmi1)
mean(bmi0)
exp(t * (mean(bmi1) - mean(bmi0)))
RRfun(24, 0.7, t) # about as good as one might expect


## this file contains the various RR functions
RRlopoff0 <- function(k, theta, t1, t2, L) {
  x0 <- 25
  a1 <- -t1
  a2 <- -t2
  ans1 <- exp(a1 * x0) * (
    pgamma(x0, k, scale = theta / (1 + a1 * theta)) -
      pgamma(L, k, scale = theta / (1 + a1 * theta))
  ) / (1 + a1 * theta)^k
  ans2 <- exp(a2 * x0) *
    pgamma(x0, k, scale = theta / (1 + a2 * theta), lower.tail = FALSE) /
    (1 + a2 * theta)^k
  ans <- ans1 + ans2
  ans / pgamma(L, k, scale = theta, lower.tail = FALSE)
}

## ## test
## RRlopoff0(24, 0.7, t1, t1, 17)

RRlopoff <- function(k, theta, t1, t2, L) {
  RRlopoff0(k, theta, t1, t2, L) /
    RRlopoff0(bmirefpop$k, bmirefpop$theta, t1, t2, 0)
}

## ## test
## RRlopoff(24,0.7,t1,t1,0)
## RRlopoff(24,0.7,t1,t1,17)
## RRlopoff(rep(24,10),rep(0.7,10),rep(t1,10),rep(t1,10),rep(17,10))

## --- flat17
flat0 <- function(k, theta, t1, t2, L, H) { # NOTE assumes H <= 25
  ## NOTE no dependence on k, theta
  x0 <- 25
  a1 <- -t1
  exp(x0 * a1) * (exp(-L * a1) - exp(-H * a1)) / (a1 * (H - L))
}

## ## delta fn test
## exp(t1 * (17.01 - 25))
## flat0(bmirefpop$k, bmirefpop$theta, t1, t1, 17, 17.01) # OK


RRflat <- function(k, theta, t1, t2, L, H) { # NOTE assumes H <= 25
  w <- pgamma(L, k, scale = theta, lower.tail = TRUE)
  dnmntr <- RRlopoff0(k, theta, t1, t2, 0)
  nmrtr <- w * flat0(k, theta, t1, t2, L, H) +
    (1 - w) * RRlopoff0(k, theta, t1, t2, L)
  nmrtr / dnmntr
}


## --- shift17
RRshift <- function(k, theta, t1, t2, L, H) { # NOTE assumes H <= 25
  x0 <- 25
  a1 <- -t1
  ans1 <- exp(a1 * (x0 + L - H)) *
    pgamma(L, k, scale = theta / (1 + a1 * theta)) /
    (1 + a1 * theta)^k
  w <- pgamma(L, k, scale = theta, lower.tail = TRUE)
  nmrtr <- (1 - w) * RRlopoff0(k, theta, t1, t2, L) +
    ans1
  dnmntr <- RRlopoff0(k, theta, t1, t2, 0)
  nmrtr / dnmntr
}

## ========= testing RR calculations
## analytical vs sampling

## 'lopoff'
bmi17lopoff <- bmi0
under <- bmi17lopoff < 17 # those under threshold
bmi17lopoff[under] <- sample(
  bmi17lopoff[!under], # those over threshold
  size = sum(under),
  replace = TRUE
)

## compare
mean(exp(BL(bmi17lopoff, t1, t2))) / mean(exp(BL(bmi0, t1, t2)))
RRlopoff(bmirefpop$k, bmirefpop$theta, t1, t1, 17) # OK

## --- flat17
## sample
bmi17flat <- bmi0
bmi17flat[under] <- runif(sum(under), 17, 25)

## compare
mean(exp(BL(bmi17flat, t1, t2))) / mean(exp(BL(bmi0, t1, t2)))
RRflat(bmirefpop$k, bmirefpop$theta, t1, t1, 17, 25) # OK


## --- shift17
## sample
bmi17shift <- bmi0
bmi17shift[under] <- bmi17shift[under] + 25 - 17

## compare
mean(exp(BL(bmi17shift, t1, t2))) / mean(exp(BL(bmi0, t1, t2)))
RRshift(bmirefpop$k, bmirefpop$theta, t1, t1, 17, 25) # OK


## === wrappers for main CF analyses

## --- use 25 as upper bound
RRflat_hi <- function(k, theta, t1, t2, L) {
  H <- 25
  ifelse(L > 0,
    RRflat(k, theta, t1, t2, L, H),
    RRlopoff(k, theta, t1, t2, 0)
  )
}

RRshift_hi <- function(k, theta, t1, t2, L) {
  H <- 25
  ifelse(L > 0,
    RRshift(k, theta, t1, t2, L, H),
    RRlopoff(k, theta, t1, t2, 0)
  )
}

## --- use (L + 25)/2 as upper bound
RRflat_lo <- function(k, theta, t1, t2, L) {
  H <- (L + 25) / 2
  ifelse(L > 0,
    RRflat(k, theta, t1, t2, L, H),
    RRlopoff(k, theta, t1, t2, 0)
  )
}


RRshift_lo <- function(k, theta, t1, t2, L) {
  H <- (L + 25) / 2
  ifelse(L > 0,
    RRshift(k, theta, t1, t2, L, H),
    RRlopoff(k, theta, t1, t2, 0)
  )
}

