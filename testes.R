#############
library('Rcpp')
library('ggplot2')
sourceCpp(file = "clownet.cpp", verbose = T)
source('clownet.R')
library('glmnet')

set.seed(142)
n = 300000
p = 1000
X = matrix(rnorm(n*p), n, p)

y = rbinom(n, 1, 1/(1 + exp(-(2*X[,1] + 2*X[,3] + 2*X[,3]))) )
t1 = Sys.time()
clo = cv.clownet(X, y, trace = T, n_lambdas = 100)
t2 = Sys.time()

t3 = Sys.time()
clo = cv.glmnet(X, y, family = "binomial", lambda = clo$lambda)
t4 = Sys.time()

dados = wooldridge::alcohol
X = model.matrix(abuse ~ ., dados)
y = dados$abuse
# 32.23152 mins
clo = cv.clownet(X, y, n_lambdas = 64)

caso1 = microbenchmark::microbenchmark(
  CLOW = cv.clownet(X, y, n_lambdas = 64),
  GLM = cv.glmnet(X, y, family = "binomial", lambda = clo$lambdas),
  times = 20
)

dados = carData::Mroz
X = model.matrix(lfp ~., dados)
y = ifelse(dados$lfp == 'yes', 1, 0)

clo = cv.clownet(X, y, n_lambdas = 64)

caso2 = microbenchmark::microbenchmark(
  CLOW = cv.clownet(X, y, n_lambdas = 64),
  GLM = cv.glmnet(X, y, family = "binomial", lambda = clo$lambdas),
  times = 100
)

dados = stevedata::TV16
dados = dados[complete.cases(dados) ,!(colnames(dados) %in% "uid")]

X = model.matrix(votetrump ~., dados)
y = dados$votetrump

clo = cv.clownet(X, y, n_lambdas = 80)

caso3 = microbenchmark::microbenchmark(
  CLOW = cv.clownet(X, y, n_lambdas = 80),
  GLM = cv.glmnet(X, y, family = "binomial", lambda = clo$lambdas),
  times = 10
)

pd_speech_features = read.csv("~/pd_speech_features.csv")[,-1]

X = model.matrix(class ~ . - 1,pd_speech_features)
y = pd_speech_features$class

clo = cv.clownet(X, y, trace = F, n_lambdas = 100)
net = glmnet(X, y, family = "binomial")

p = profmem::profmem({
  cv.clownet(X,y,trace = F, n_lambdas = 100)
})

clo_mem = p$bytes |> sum(na.rm = T)

g =  profmem::profmem({
  cv.glmnet(X,y, family = "binomial")
})

glm_mem = g$bytes |> sum(na.rm = T)