setwd('D:/Documents/R/Statistics and Econometrics/StatisticR')

library(devtools)
load_all()

# test for large sample
n <- 50000
df <- data.frame(x=rnorm(n),
                 y=rnorm(n))
ols.estimate(df$y, df['x'], T)

# ----- test for dummy generation
x <- c(1, 2, 3, 2, 4, 1)
print(gen.dummy(x, prefix='x'))
print(gen.dummy(x, prefix='x', first.rm=TRUE))
x <- c(1, 2, 1, 2, 1)
print(gen.dummy(x, prefix='x', first.rm=TRUE))


# ----- test for within estimate
n <- 500
df <- data.frame(x=rnorm(n),
                 indic=round(runif(n, 1, 4)))
y <- apply(df, 1, sum)

df$x[5] <- NA

print(ols.estimate(y, df['x']))

dummy <- gen.dummy(df$indic, prefix='indic')
X <- cbind(df['x'], dummy[, 1:(dim(dummy)[2]-1)])
print(ols.estimate(y, X))
print(ols.within.estimate(y, df['x'], df$indic))


# ----- test for numeric.gr
foo <- function(theta, c){
  return(theta[1]*c+theta[2]**2)
}

foo.gr <- numeric.gr(foo)
print(foo.gr(c(2,3), 2))
# 2 6

# ----- test for clean4regression
X <- data.frame(x=1:6,
                y=c(2, 2, 2, 1, 3, 2),
                z=c(1, 6, 7, NA, NA, 3))
W <- data.frame(w=1:6)
args <- list(X=as.matrix(X), W=as.matrix(W))
print(clean4regression(args))


# ----- test for linear regression, model.res.export
ols.model.res <- ols.estimate(regdata$y, regdata[c('x1', 'x2', 'x3', 'x4')], robust = T)
mle.model.res <- mle.linear.estimate(regdata$y, regdata[c('x1', 'x2', 'x3', 'x4')])

res.table <- model.res.table.export(list(ols.model.res))
print(res.table)
res.table <- model.res.table.export(list(ols.model.res, mle.model.res), digits=4)
print(res.table)


# ----- test for binary models
all.params <- c('x1', 'x2', 'x3', 'x4')
y <- ifelse(regdata$y > 18, 1, 0)

probit.res <- mle.probit.estimate(y, regdata[c('x1', 'x3', 'x4')])
logit.res <- mle.logit.estimate(y, regdata[c('x1', 'x3', 'x4')])

res.table <- model.res.table.export(list(probit.res, logit.res), all.params)
print(res.table)


# ----- test for duration models
all.params <- c('workprg', 'priors', 'tserved', 'felon', 'alcohol',
                'drugs', 'black', 'married', 'educ', 'age', '_const', 'lnp', 'gamma')

X <- recid[c('workprg', 'priors', 'tserved', 'felon', 'alcohol',
             'drugs', 'black', 'married', 'educ', 'age')]

exp.res <- mle.exp.estimate(recid$durat, X, 1-recid$cens)
weibull.res <- mle.weibull.estimate(recid$durat, X, 1-recid$cens)
gompertz.res <- mle.gompertz.estimate(recid$durat, X, 1-recid$cens)

# exponential model result from Stata
# args <- list(t=df$durat, X=as.matrix(X), d=1-df$cens)
# stata_params <- c(.09558006, .09133702, .01440089, -.31222409, .46767064, .29416476,
#                   .47567491, -.1519512, -.02421234, -.00391121, -4.169213)
# print(exp.lnlike(stata_params, args), digits=10)

# gompertz model result from stata
# stata_params <- c(.08532224, .08691561, .01288472, -.28549804, .42928766, .27253692,
#                   .43281454, -.15295268, -.02185928, -.00356317, -3.6110097, -.02170298)
# print(gompertz.lnlike(stata_params, args), digits=10)

res.table <- model.res.table.export(list(exp.res, weibull.res, gompertz.res), all.params)
print(res.table)


# ----- test for heckman
y <- womenwk$lw
X <- womenwk[c('education', 'age', 'children')]
z <- ifelse(is.na(y), 0, 1)
W <- womenwk[c('age', 'married', 'children', 'education')]

# in this model, optimization result would be better if setting the initial params=[1e-4, 1e-4, ...]
# lnlike: 1052.857 -> 1044.652
heckman.res <- mle.heckman.estimate(y,X,z,W)
print(heckman.res)


# ----- test for duration models with selection
all.params <- c('W.democ', 'W.autoc', 'W.tennewL', 'W.tennewAL', 'W.total', 'W.majpow',
                'W.lntropen', 'W._const', 'X.tenl', 'X.democ', 'X.tendem', 'X.rdpopl',
                'X.rwin', 'X._const', 'atanhalpha', 'lnp')

t <- wardata$wsurv
X <- wardata[c('tenl', 'democ', 'tendem', 'rdpopl', 'rwin')]
z <- wardata$enter
W <- wardata[c('democ', 'autoc', 'tennewL', 'tennewAL', 'total', 'majpow', 'lntropen')]
d <- wardata$rtcensor

exp.res <- mle.exp.selection.estimate(t, X, z, W, d)
weibull.res <- mle.weibull.selection.estimate(t, X, z, W, d)
res.table <- model.res.table.export(list(exp.res, weibull.res), all.params)
print(res.table)
