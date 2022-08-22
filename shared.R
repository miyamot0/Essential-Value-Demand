
createFormattingString = function(evFrame, row = 1) {

  evVal = round(evFrame[row, 'estimate'], 6) %>% as.numeric()
  sVal  = round(evFrame[row + 6, 'estimate'], 2) %>% as.numeric()

  if (row > 1) {
    evVal = round(evFrame[row, 'estimate'] + evFrame[1, 'estimate'], 6) %>% as.numeric()
    sVal  = round(evFrame[row + 6, 'estimate'] + evFrame[4, 'estimate'], 2) %>% as.numeric()
  }

  stringOut = bquote(alpha == .(evVal)*k^.(sVal))
  return(stringOut)
}

lambertW = function(z,b=0,maxiter=10,eps=.Machine$double.eps,min.imag=1e-9) {
  if (any(round(Re(b)) != b))
    stop("branch number for W must be an integer")
  if (!is.complex(z) && any(z<0)) z=as.complex(z)
  ## series expansion about -1/e
  ##
  ## p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
  ## w = (11/72)*p;
  ## w = (w - 1/3).*p;
  ## w = (w + 1).*p - 1
  ##
  ## first-order version suffices:
  ##
  w = (1 - 2*abs(b))*sqrt(2*exp(1)*z + 2) - 1
  ## asymptotic expansion at 0 and Inf
  ##
  v = log(z + as.numeric(z==0 & b==0)) + 2*pi*b*1i;
  v = v - log(v + as.numeric(v==0))
  ## choose strategy for initial guess
  ##
  c = abs(z + exp(-1));
  c = (c > 1.45 - 1.1*abs(b));
  c = c | (b*Im(z) > 0) | (!Im(z) & (b == 1))
  w = (1 - c)*w + c*v
  ## Halley iteration
  ##
  for (n in 1:maxiter) {
    p = exp(w)
    t = w*p - z
    f = (w != -1)
    t = f*t/(p*(w + f) - 0.5*(w + 2.0)*t/(w + f))
    w = w - t
    if (abs(Re(t)) < (2.48*eps)*(1.0 + abs(Re(w)))
        && abs(Im(t)) < (2.48*eps)*(1.0 + abs(Im(w))))
      break
  }
  if (n==maxiter) warning(paste("iteration limit (",maxiter,
                                ") reached, result of W may be inaccurate",sep=""))
  if (all(Im(w)<min.imag)) w = as.numeric(w)
  return(w)
}

GetAnalyticPmaxFallback <- function(K_, A_, Q0_) {
  result <- optimx::optimx(par = c((1/(Q0_ * A_ * K_^1.5)) * (0.083 * K_ + 0.65)),
                           fn = function(par, data) {
                             abs((log((10^data$K)) * (-data$A * data$Q0 * par[1] * exp(-data$A * data$Q0 * par[1]))) + 1)
                           },
                           data = data.frame(Q0 = Q0_,
                                             A = A_,
                                             K = K_),
                           method = c("BFGS"),
                           control=list(maxit=2500))

  return(result$p1)
}

GetAnalyticPmax <- function(Alpha, K, Q0) {
  if (K <= exp(1)/log(10)) {
    return (NaN);
    return (GetAnalyticPmaxFallback(K, Alpha, Q0));
  } else {
    return (-lambertW(z = -1/log((10^K))) / (Alpha * Q0));
  }
}
