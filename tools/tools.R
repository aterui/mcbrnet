devtools::load_all(".")
Rcpp::sourceCpp("tools/src/deriv.cpp")

s <- 4
b <- 2
l <- 5
p <- 2
ts <- 100

alpha <- foodweb(s, b, l, theta = 1)
A <- kronecker(diag(p), alpha)

c0 <- matrix(c(0,0,0,0), 2, 2)
C <- kronecker(c0, diag(s))

R <- findr(alpha, k0 = 10)
r <- rep(R[, 1], 2)

# pseudo time steps
psi <- cumsum(rexp(n = ceiling(1 + 100), 0.1))
psi <- psi[psi <= ts]
time <- seq(0, ts, by = 0.5)
t_psi <- sapply(psi, FUN = function(x) which.min(abs(time - x)))

# create dummy time-series for signals
signal <- data.frame(time = time,
                     psi = rep(0, length(time)))

signal$psi[t_psi] <- 1

# linear interpolation function
input <- stats::approxfun(signal, rule = 2)

parms <- list(r = r,
              e = rep(runif(p, 0, 5), each = s),
              A = A,
              C = C)

# initial values for a state variable
n_init <- runif(s * p, 0, 5)

# time-series
times <- seq(0, ts, by = 0.01)


# ode ---------------------------------------------------------------------

tr <- 0.01
eventfun <- function(t, n, parms) {
  n <- ifelse(n <= tr, 0, n)
  return(n)
}

rootfun <- function(t, n, parms) {
  return(n - tr)
}

out <- deSolve::ode(y = n_init, times = times, func = deriv, parms = parms,
                    events = list(func = eventfun, root = TRUE),
                    rootfun = rootfun)
plot(out)
cbind(rep(R[,2], 2), out[nrow(out), -1])
