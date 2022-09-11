#############
library('Rcpp')
library('ggplot2')
sourceCpp(file = "clownet.cpp", verbose = T)
source('clownet.R')
library('glmnet')

set.seed(142)

######### Cuidado ########

# n = 300000
# p = 1000
# X = matrix(rnorm(n*p), n, p)
# 
# y = rbinom(n, 1, 1/(1 + exp(-(2*X[,1] + 2*X[,3] + 2*X[,3]))) )
# t1 = Sys.time()
# clo = cv.clownet(X, y, trace = T, n_lambdas = 100)
# t2 = Sys.time()
# 
# t3 = Sys.time()
# clo = cv.glmnet(X, y, family = "binomial")
# t4 = Sys.time()
# 
# betas_esparsos = coef(clo, s = "lambda_min")
# 
# sum(abs(betas_esparsos))
# which(abs(betas_esparsos) >  0.05)

############### CONJUNTOS DE DADOS ###############
dados = wooldridge::alcohol
X = model.matrix(abuse ~ ., dados)
y = dados$abuse
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

#### Parkison disease - Utilizado no relatório ####
#FONTE: https://archive-beta.ics.uci.edu/ml/datasets/parkinson+s+disease+classification#Descriptive

pd_speech_features = read.csv("~/pd_speech_features.csv")[,-1]

amostra = sample(nrow(pd_speech_features), trunc(0.7*nrow(pd_speech_features)))
X = model.matrix(class ~ . - 1,pd_speech_features)
y = pd_speech_features$class

caso4 = microbenchmark::microbenchmark(
  CLOW = cv.clownet(X, y, n_lambdas = 100),
  GLM = cv.glmnet(X, y, family = "binomial", lambda = clo$lambdas),
  times = 10
)

### Modelos
clo = cv.clownet(X[amostra,], y[amostra], trace = F, n_lambdas = 100)
net = cv.glmnet(X[amostra,], y[amostra], family = "binomial")

### PERFORMANCE

preditos_cn = predict(clo, s = "lambda_min", newdata = cbind(1, X[-amostra,]))
classes_cn = ifelse(1/(1 + exp(-preditos_cn)) > 0.5, 1, 0)

clocc = caret::confusionMatrix(as.factor(pd_speech_features$class[-amostra]), 
                       as.factor(classes_cn))


clocc$table |> xtable::xtable()

preditos_glm = predict(net, s = "lambda.min", newx = cbind(X[-amostra,]))
classes_glm = ifelse(1/(1 + exp(-preditos_glm)) > 0.5, 1, 0)

clogl = caret::confusionMatrix(as.factor(pd_speech_features$class[-amostra]), 
                       as.factor(classes_glm))

clogl$table |> xtable::xtable()

### ZEROS

(coef(net, s = "lambda.min") == 0) |> sum()
(coef(clo, s = "lambda_min") == 0) |> sum()

coef(net, s = "lambda.min") |> abs() |> sum()
coef(clo, s = "lambda_min") |> abs() |> sum()


### MEMÓRIA

p = profmem::profmem({
  cv.clownet(X,y,trace = F, n_lambdas = 100)
})

clo_mem = p$bytes |> sum(na.rm = T)

g =  profmem::profmem({
  cv.glmnet(X,y, family = "binomial")
})

glm_mem = g$bytes |> sum(na.rm = T)
