library('Rcpp')
library('ggplot2')

sourceCpp(file = "clownet.cpp")

cv.clownet = function(X, y, n_lambdas = 50, n_blocos = 10, trace = FALSE){

  ####First time
  obj = clownet(X, y, n_lambdas = n_lambdas)
  lambdas = obj$lambdas
  
  n = nrow(X)
  
  compo_de_desvio = matrix(0, n_blocos, n_lambdas)
  ep = trunc(n/n_blocos)
  preditos = matrix(NA, nrow = n, ncol = n_lambdas)
  
  prob_min = 10^(-4)
  prob_max = 1 - prob_min
  
  #Matriz
  X = X[, !colnames(X) %in% "(Intercept)"]
  X = scale(X)
  attr_center = attr(X, "scaled:center")
  attr_scale = attr(X, "scaled:scale")
  X = as.matrix(cbind(Intercept = rep(1, n), X))
  
  shuffle = sample(nrow(X))
  X = X[shuffle,]
  y = y[shuffle]
  
  for(i in 1:n_blocos){
    if(i == n_blocos) bloco = ((i-1)*ep+1):n
    else bloco = ((i-1)*ep+1):(i*ep)
      
    X_treino = X[-bloco,]
    y_treino = y[-bloco]
    
    aux = CPP_clownet(X_treino, y_treino, 
                                attr_scale, attr_center, lambdas = lambdas,
                                cv = TRUE,
                                n_lambda = n_lambdas, trace = trace) 
    
    X_bloco = X[bloco,]

    aux = CPP_predict(X_bloco, aux$betas, n_lambdas, y[bloco])
    
    preditos[bloco,] = aux
    
    compo_de_desvio[i,] = CPP_calc_deviance_vector(aux, y[bloco])
  }
  
  preditos = preditos[-1,]
  preditos = pmin(pmax(preditos,prob_min),prob_max)
  
  lu = obj$lp/n + qchisq(0.975, 1)*apply(compo_de_desvio, 2, sd)/n
  ld = obj$lp/n - qchisq(0.975, 1)*apply(compo_de_desvio, 2, sd)/n
  
  out = list(lambdas = lambdas, dev.x = obj$lp/n, dev.up = lu,
             dev.lo = ld, centro_coef = attr_center, escala_coef = attr_scale)
  class(out) = c('cv.clownet', 'list')
  return(out)
}

X = model.matrix(lfp ~ ., data = carData::Mroz)
y = ifelse(carData::Mroz$lfp == "yes", 1, 0)

modelo = cv.clownet(X, y)

plot(modelo$lambdas |> log(), modelo$dev.x)

###### PLOTAR

df_plot = data.frame(dev_min = modelo$dev.lo,
dev_mean = modelo$dev.x,
dev_max = modelo$dev.up,
lambda = modelo$lambdas)

ggplot(df_plot, aes(x = log(lambda), y = dev_mean)) + 
  geom_point(size = 0.9) +
  #geom_ribbon(aes(ymin=dev_min,ymax=dev_max), fill = "tomato", col = "red", alpha=0.4) +
  geom_errorbar(aes(ymin=dev_min,ymax=dev_max), col = "red") +
  geom_line(aes(x = log(lambda), y = dev_mean), col = "black", linetype = "dashed") +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x = expression(log(lambda)), y = "Deviance")

library(glmnet)

x = cv.glmnet(X,y, family = binomial())

X = model.matrix(lfp ~ ., data = carData::Mroz)
y = ifelse(carData::Mroz$lfp == "yes", 1, 0)

profvis::profvis({clownet(X, y)})

library(microbenchmark)

bm = microbenchmark(
  GLMNET = glmnet(X, y, family = binomial()),
  CLOWNET = clownet(X, y),
  times = 5
)


clownet = function(X, y, n_lambdas = 50, trace = F){
  prob_min = 10^(-4)
  prob_max = 1 - prob_min
  
  X = X[, !colnames(X) %in% "(Intercept)"]
  X = scale(X)
  attr_center = attr(X, "scaled:center")
  attr_scale = attr(X, "scaled:scale")
  X = as.matrix(cbind(Intercept = rep(1, nrow(X)), X))
  
  obj = CPP_clownet(X, y, attr_scale, attr_center, n_lambda = n_lambdas, trace = F)

  preditos = CPP_predict(X, obj$betas, n_lambdas, y)
  preditos = pmin(pmax(preditos,prob_min),prob_max)
  #lp = CPP_calc_deviance_vector(preditos, y)
  lp = apply(preditos, 2, \(x) sum((-2*((1-y)*log(1-x) + y*log(x)))))
  
  out = list(lambdas = obj$lambdas, betas = obj$betas, preditos = preditos, lp = lp,
             centro_coef = attr_center, escala_coef = attr_scale)
  class(out) = c('list', 'clownet')
  return(out)
}

plot(log(clo$lambdas), clo$lp)

coef.clownet = function(obj, s){
  
}

X = model.matrix(lfp ~ ., data = carData::Mroz)
X = X[, !colnames(X) %in% "(Intercept)"]
X = scale(X)
attr_center = attr(X, "scaled:center")
attr_scale = attr(X, "scaled:scale")
X = as.matrix(cbind(Intercept = rep(1, nrow(X)), X))
y = ifelse(carData::Mroz$lfp == "yes", 1, 0)
X_net = model.matrix(lfp ~ ., data = carData::Mroz)

obj = CPP_clownet(X, y, attr_scale, attr_center, trace = F)


library('microbenchmark')

microbenchmark(
  CLOWNET = CPP_clownet(X, y, attr_scale, attr_center, trace = F),
  GLMNET = glmnet(X_net, y, family = binomial()),
  times = 100
) |> plot()
