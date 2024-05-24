#----------------------------------------------------------------------------
#' Generate univariate signals of different type
#'
#' Using code from the archived \code{wmtsa} package
#'
#' @param name name of the signal to be generated
#' @param n length of the series
#' @param snr desired signal-to-noise ratio
"make.signal" <- function(name, n=1024, snr=Inf)
{

  ".wave.demo.signals" <- c("dirac", "kronecker", "heavisine", "bumps", "blocks",
                            "doppler", "ramp", "cusp", "crease", "sing", "hisine",
                            "losine", "linchirp", "twochirp", "quadchirp",
                            "mishmash1", "mishmash2", "mishmash3", "levelshift",
                            "jumpsine", "gauss", "patches",
                            "linear", "quadratic", "cubic")

  x <- (0:(n-1.))/n
  z <- switch(name,
              dirac=n*(x == floor(.37*n)/n),
              kronecker=(x == floor(.37*n)/n),
              heavisine=4*sin(4*pi*x)-sign(x-.3)-sign(.72-x),
              bumps={
                pos <- c(.1, .13, .15, .23, .25, .4, .44, .65, .76, .78, .81)
                hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
                wth <- c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008,.005)
                y <- rep(0, n)
                for(j in 1:length(pos)) y <- y+hgt[j]/(1+abs((x-pos[j]))/wth[j])^4
                y
              },
              blocks={
                pos <- c(.1, .13, .15, .23, .25, .4, .44, .65, .76, .78, .81)
                hgt <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1,2.1, -4.2)
                y <- rep(0, n)
                for(j in 1:length(pos)) y <- y+(1+sign(x-pos[j]))*hgt[j]/2
                y
              },
              doppler=sqrt(x*(1-x))*sin((2*pi*1.05)/(x+.05)),
              ramp=x-(x >= .37),
              cusp=sqrt(abs(x-.37)),
              crease=exp(-4*abs(x-.5)),
              sing=1/abs(x-(floor(n*.37)+.5)/n),
              hisine=sin(pi*n*.6902*x),
              midsine=sin(pi*n*.3333*x),
              losine=sin(pi*n*.03*x),
              linchirp=sin(.125*pi*n*x^2),
              twochirp=sin(pi*n*x^2) + sin((pi/3)*n*x^2),
              quadchirp=sin((pi/3)*n*x^3),
              # QuadChirp + LinChirp + HiSine
              mishmash1=sin((pi/3)*n*x^3) + sin(pi*n*.6902*x) + sin(pi*n*.125*x^2),
              # QuadChirp + LinChirp + HiSine + Bumps
              mishmash2={		# wernersorrows
                y   <- sin(pi*(n/2)*x^3)+sin(pi*n*.6902*x)+sin(pi*n*x^2)
                pos <- c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
                hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
                wth <- c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008,.005)
                for(j in 1:length(pos)) y <- y + hgt[j]/(1+abs((x-pos[j])/wth[j]))^4
                y
              },
              # QuadChirp + MidSine + LoSine + Sing/200.
              mishmash3=sin((pi/3)*n*x^3) + sin(pi*n*.3333*x) + sin(pi*n*.03*x) +
                (1/abs(x-(floor(n*.37)+.5)/n))/(200.*n/512.),
              gauss=dnorm(x, .3, .025),
              jumpsine=10.*(sin(4*pi*x) + as.numeric(x >= 0.625 & x < 0.875)),
              levelshift=as.numeric(x >= 0.25 & x < 0.39),
              patches={
                if(n<16) stop("n must be >= 16 to generate patches\n")
                J <- logb(n, base=2)
                y <- rep(0., n)
                for(j in 0:(J-4.)) y[(1:2^j)+3.*2.^(j+2.)] <- 1.
                y
              },
              linear=2.*x-1.,
              quadratic=4. * (1. - x) * x,
              cubic=64. * x * (x - 1.) * (x - .5) / 3.,
              stop("Unknown signal name.  Allowable names are:\n",
                   paste(.wave.demo.signals, collapse = ", ")))

  if (snr > 0)
    z <- z + rnorm(n) * sqrt(var(z)) / snr

  z
}
