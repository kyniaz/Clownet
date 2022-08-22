library('Rcpp')
library('ggplot2')

sourceCpp(file = "rcc_lambda.cpp")

X = model.matrix(lfp ~ ., data = carData::Mroz)
X = X[, !colnames(X) %in% "(Intercept)"]
non_scaled_X = X
X = scale(X)

attr_center = attr(X, "scaled:center")
attr_scale = attr(X, "scaled:scale")

X = as.matrix(cbind(Intercept = rep(1, nrow(X)), X))

y = ifelse(carData::Mroz$lfp == "yes", 1, 0)

obj = clownet(X, y, attr_scale, attr_center, trace = F)

obj = as.data.frame(obj)
colnames(obj) = c(colnames(X), "lambda")

l1 = rowSums(abs(obj[,-ncol(obj)]))
obj$l1 = l1

ci = length(colnames(obj))

plt = eval(parse(text = paste0('ggplot(obj) + ', paste0('geom_line(aes(x = l1, y = ',
                                                              colnames(obj)[-c(1, ci - 1, ci)],',
                                                                    color = "',colnames(obj)[-c(1, ci - 1, ci)],'"))',
                                                              collapse = '+'))))

plt = plt + theme_minimal() + scale_colour_brewer(palette = "Set1")

plt

####################
##
soft = function(z, y){
  if (z > y) return(z - y)
  else if (z < -y) return(z + y)
  else return(0)
}

lq = function(par, X, w, z){
  return((-1/(2*nrow(X)))*sum(w*(z - par[1] - X[,-1]%*%par[-1])^2))
}

log_veros = function(par,y, X){
  return((1/nrow(X))*sum(y*(par[1] + X[,-1]%*%par[-1])  - log(1 + exp(par[1] + X[,-1]%*%par[-1]))))
}

### CLOWNET
clownet = function(X, y, maxit = 100, tol = 10^(-8), alpha = 1){
  X = X[, !colnames(X) %in% "(Intercept)"]
  non_scaled_X = X
  X = scale(X)
  
  attr_center = attr(X, "scaled:center")
  attr_scale = attr(X, "scaled:scale")
  
  X = as.matrix(cbind(Intercept = rep(1, nrow(X)), X))
  
  betas = c(glm.fit(X[,1], y, family = binomial())$coefficients, rep(0, ncol(X) - 1))
  
  max_l = (1/nrow(X)) *  abs(sapply(2:ncol(X), function(col) sum(X[,col]%*%y)))
  
  lambdas = seq(0, max(max_l), length.out= 200)
  
  storage_lambdas = data.frame(matrix(NA, nrow = length(lambdas), ncol = ncol(X)))
  colnames(storage_lambdas) = colnames(X)
  risco = vector(mode = "numeric", length = length(lambdas))
  
  #funs
  
  veros_parcial = function(par, j, s) {
    if(j != 1) lasso = - s*abs(par)
    else lasso = 0
    betas_lasso[j] = par
    return((-1/(2*nrow(X)))*sum(w*(z - betas_lasso[1] - X[,-1]%*%betas_lasso[-1])^2) + lasso)
  }
  
  grad_veros_parcial = function(j, betas_lasso){
    if(j==1){
      return((-1/(2*nrow(X)))*sum(2*w*(z - X%*%betas_lasso)))
    }
    else{
      return((-1/(2*nrow(X)))*sum(2*w*(z - X%*%betas_lasso)*X[,j]))
    }
  }
  
  for(l in length(lambdas):1){
    cat('Lambda : ', lambdas[l], '\n')
    erro = 1
    tol = 10^(-8)
    it = 1
    maxit = 100
    alpha = 1
    s = lambdas[l]
    
    while(erro > tol & it < maxit){
      
      eta = X%*%betas
      
      mu = 1/(1 + exp(-eta))
      
      logit_mu = log(mu/(1 - mu))
      
      #Z e W
      z = logit_mu + (y - mu)* 1/(mu*(1-mu))
      w = (mu*(1 - mu)) #canonica
      
      
      #betas_novos = betas
      betas_lasso = betas
      erro_lasso = 1
      it_lasso = 1
      
      while(erro_lasso > tol & it_lasso < maxit){
        
        betas_ant = betas_lasso
        
        for(j in 1:length(betas_lasso)){
          #zed = optimise(veros_parcial, interval = c(-10, 10), maximum = T, j = j, s = s)$maximum
          zed = betas_lasso[j] - 1*grad_veros_parcial(j, betas_lasso)
            
          if( j != 1) betas_lasso[j] = soft(zed, s)
          else { 
            if(l >= 199){
              cat("mu: ", mean(mu), "\n")
              cat("eta: ", mean(eta), "\n")
              cat("z: ", mean(z), "\n")
              cat("w: ", mean(w), "\n")
            }
            betas_lasso[j] = zed 
          }
          
          eta = X%*%betas_lasso
          
          mu = 1/(1 + exp(-eta))
          
          logit_mu = log(mu/(1 - mu))
          
          #Z e W
          z = logit_mu + (y - mu)* 1/(mu*(1-mu))
          w = (mu*(1 - mu)) #canonica
        }
        
        erro_lasso = max(abs(betas_lasso - betas_ant))
        
        it_lasso = it_lasso + 1
      }
      
      betas_novos = betas_lasso
      
      erro = sum(((betas_novos - betas)/ifelse(betas == 0, betas + 0.001, betas))^2)
      
      if(erro < tol){
        break
      }
      
      betas = betas_novos
      it =  it + 1
    }
    
    storage_lambdas[l,] = c(betas[1], betas[-1]/attr_scale)
    
    pred = 1/(1 + exp(- (betas[-j]%*%t(X[,-j]))))
    risco[l] = mean((y - pred)^2)
  }
  out = structure(storage_lambdas, lambdas = lambdas, it_fs = it, erro = erro, risco = risco,
                  matriz_planejamento = cbind(non_scaled_X, y))
  class(out) = c(class(out), "clownet")
  return(out)
}

plot.clownet = function(modelo, type = ""){
  
}

library('car')
library('ggplot2')

data('Mroz')

X = model.matrix(lfp ~ ., data = carData::Mroz)
y = ifelse(carData::Mroz$lfp == "yes", 1, 0)

modelo_cn = clownet(X, y)

modelo_cn$l1 = rowSums(abs(modelo_cn))
modelo_cn$lambdas = attr(modelo_cn, 'lambdas')

ci = length(colnames(modelo_cn))

plt = eval(parse(text = paste0('ggplot(modelo_cn) + ', paste0('geom_line(aes(x = l1, y = ',
                                                              colnames(modelo_cn)[-c(1, ci - 1, ci)],',
                                                                    color = "',colnames(modelo_cn)[-c(1, ci - 1, ci)],'"))',
                                                              collapse = '+'))))

plt = plt + theme_minimal() + scale_colour_brewer(palette = "Set1")

plt

### GLMNET
library(glmnet)

X = model.matrix(lfp ~ ., data = Mroz)
X = X[, !colnames(X) %in% "(Intercept)"]

modelo = glmnet(X, y, family = binomial())

coef(modelo, s = 0)

plot(modelo, label = T)

############ TESTE

X = model.matrix(lfp ~ ., data = carData::Mroz)
X = X[, !colnames(X) %in% "(Intercept)"]
X = scale(X)
attr_center = attr(X, "scaled:center")
attr_scale = attr(X, "scaled:scale")
X = as.matrix(cbind(Intercept = rep(1, nrow(X)), X))
y = ifelse(carData::Mroz$lfp == "yes", 1, 0)

obj = CPP_clownet(X, y, attr_scale, attr_center, trace = F)


library('microbenchmark')

microbenchmark(
  CLOWNET = clownet(X, y, attr_scale, attr_center, trace = F),
  GLMNET = glmnet(X_net, y, family = binomial()),
  times = 10
)


cv.clownet = function(X, y, n_lambdas = 200, n_blocos = 10, trace = FALSE){
  
  #Matriz
  X = X[, !colnames(X) %in% "(Intercept)"]
  X_nao_escalada = X
  X = scale(X)
  
  #Mudar
  attr_center = attr(X, "scaled:center")
  attr_scale = attr(X, "scaled:scale")
  n = nrow(X)
  
  X = as.matrix(cbind(Intercept = rep(1, n), X))
  
  lambda = CPP_calc_lambda_double(X, y)
  
  lambdas = seq(0, max(lambda), length.out= n_lambdas)
  amostra = sample(n, n, replace = F)
  
  compo_de_desvio = matrix(0, n_blocos, n_lambdas)
  
  ep = trunc(n/n_blocos)
  
  preditos = vector(mode = "list", n_lambdas)
  
  for(i in 1:n_blocos){
    if(i == n_blocos) bloco = ((i-1)*ep+1):n
    else bloco = ((i-1)*ep+1):(i*ep)
      
    X_treino = X[-bloco,]
    y_treino = y[-bloco]
    
    aux = CPP_clownet(X_treino, y_treino, attr_center, 
                                attr_scale, lambdas = lambdas,
                                cv = TRUE,
                                n_lambda = 50, trace = trace) 
    
    preditos[[i]] = CPP_predict(X_nao_escalada[-bloco,], aux$betas, n_lambdas, y_treino)
    compo_de_desvio[i,] = CPP_calc_deviance_vector(preditos[[i]], y_treino)
  }
  
  compo_de_desvio = as.data.frame(t(compo_de_desvio))
  colnames(compo_de_desvio) = paste0("bloco: ", 1:n_blocos)
  compo_de_desvio$lambda = rev(lambdas)
    
  return(compo_de_desvio)
}

X = model.matrix(lfp ~ ., data = carData::Mroz)
y = ifelse(carData::Mroz$lfp == "yes", 1, 0)

modelo = cv.clownet(X, y)

###### PLOTAR

erro_padrao = apply(modelo[,1:10], 1, sd)
media = apply(modelo[,1:10], 1, mean)

df_plot = data.frame(dev_min = media - qnorm(0.95)*erro_padrao,
dev_mean = apply(modelo[,1:10], 1, mean),
dev_max =  media + qnorm(0.95)*erro_padrao,
lambda = modelo$lambda)

ggplot(df_plot, aes(x = log(lambda), y = dev_mean)) + 
  geom_point(size = 0.9) +
  geom_ribbon(aes(ymin=dev_min,ymax=dev_max), fill = "tomato", col = "red", alpha=0.4) +
  geom_line(aes(x = log(lambda), y = dev_mean), col = "black", linetype = "dashed") +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5)) +
  labs(x = expression(log(lambda)), y = "Deviance")

