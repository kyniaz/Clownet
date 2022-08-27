library('Rcpp')
library('ggplot2')

sourceCpp(file = "clownet.cpp")

clownet = function(X, y, n_lambdas = 50, trace = F){
  prob_min = 10^(-4)
  prob_max = 1 - prob_min
  chamado = colnames(X)
  
  X = X[, !colnames(X) %in% "(Intercept)"]
  X = scale(X)
  attr_center = attr(X, "scaled:center")
  attr_scale = attr(X, "scaled:scale")
  X = as.matrix(cbind(Intercept = rep(1, nrow(X)), X))
  
  obj = CPP_clownet(X, y, attr_scale, attr_center, n_lambda = n_lambdas, trace = F)
  
  preditos = CPP_predict(X, obj$betas, n_lambdas + 1, y)
  preditos = pmin(pmax(preditos,prob_min),prob_max)
  #lp = CPP_calc_deviance_vector(preditos, y)
  lp = apply(preditos, 2, \(x) sum((-2*((1-y)*log(1-x) + y*log(x)))))
  
  out = list(chamado = chamado,
             lambdas = obj$lambdas, betas = obj$betas, preditos = preditos, lp = lp,
             centro_coef = attr_center, escala_coef = attr_scale)
  class(out) = c('list', 'clownet')
  return(out)
}

cv.clownet = function(X, y, n_lambdas = 50, n_blocos = 10, trace = FALSE){

  ####First time
  obj = clownet(X, y, n_lambdas = n_lambdas)
  lambdas = obj$lambdas
  
  n = nrow(X)
  
  compo_de_desvio = matrix(0, n_blocos, n_lambdas + 1)
  ep = trunc(n/n_blocos)
  preditos = matrix(NA, nrow = n, ncol = n_lambdas + 1)
  
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

    aux = CPP_predict(X_bloco, aux$betas, n_lambdas + 1, y[bloco])
    
    preditos[bloco,] = aux
    
    compo_de_desvio[i,] = CPP_calc_deviance_vector(aux, y[bloco])
  }
  
  preditos = preditos[-1,]
  preditos = pmin(pmax(preditos,prob_min),prob_max)
  
  lu = obj$lp/n + qchisq(0.975, 1)*apply(compo_de_desvio, 2, sd)/n
  ld = obj$lp/n - qchisq(0.975, 1)*apply(compo_de_desvio, 2, sd)/n
  
  out = list(lambdas = lambdas, dev.x = obj$lp/n, dev.up = lu, betas = obj$betas,
             dev.lo = ld, centro_coef = attr_center, escala_coef = attr_scale)
  class(out) = c('cv.clownet', 'clownet', 'list')
  return(out)
}

X = model.matrix(lfp ~ ., data = carData::Mroz)
y = ifelse(carData::Mroz$lfp == "yes", 1, 0)

modelo = cv.clownet(X, y)


library(glmnet)

x = glmnet(X,y, family = binomial())

X = model.matrix(lfp ~ ., data = carData::Mroz)
y = ifelse(carData::Mroz$lfp == "yes", 1, 0)
clo = clownet(X, y)

plot(log(clo$lambdas), clo$lp)

extract.clownet = function(obj, s){
  
  if(s > max(obj$lambda) | s < min(obj$lambda)){
    stop("choose a valid value for s.")
  }
  
  ls = which.min(abs(s - obj$lambdas))
  
  betas = obj$betas[ls,]
  betas_out = numeric(length(betas))
  centro = obj$centro_coef
  escala = obj$escala_coef
  
  #desescalar coeficientes
  val = 0;
  for(l in 2:length(betas)){
    val = val + betas[l]*(centro[l-1]/escala[l-1]);
    betas_out[l] = betas[l]/escala[l-1];
  }
  
  betas_out[1] = betas[1] - val;
  
  return(betas_out)
}

#Extrair coeficientes;
coef.clownet = function(obj, s){
  do.call("extract.clownet", list(obj, s))
}

coef(clo, s = 0.001) |> abs() |> sum()

plot.cv.clownet = function(obj){
  df_plot = data.frame(dev_min = obj$dev.lo,
                       dev_mean = obj$dev.x,
                       dev_max = obj$dev.up,
                       lambda = obj$lambdas)
  
  ggplot2::ggplot(df_plot, aes(x = log(lambda), y = dev_mean)) + 
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_errorbar(aes(ymin=dev_min,ymax=dev_max), col = "red") +
    ggplot2::geom_line(aes(x = log(lambda), y = dev_mean), col = "black", linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title=element_text(hjust=0.5)) +
    ggplot2::labs(x = expression(log(lambda)), y = "Deviance")
}

plot(modelo)

plot.clownet = function(obj){
  aux = data.frame(obj$betas)
  aux = aux[,-1]
  chamado = obj$chamado[-1]
  colnames(aux) = chamado
  
  aux$l1 = rowSums(abs(aux))
  
  plt = eval(parse(text = paste0('ggplot2::ggplot(aux) + ', 
                                 paste0('ggplot2::geom_line(aes(x = l1, y = ',
                                    chamado,',
                                    color = "',chamado,'"))',
                                 collapse = '+'))))
  
  plt = plt + ggplot2::theme_minimal() + ggplot2::scale_colour_brewer(palette = "Set1") +
  ggplot2::labs(x = "Soma L1", y = "Valor de Coeficientes")
  
  plt
}

plot(clo)
#plot cv

#plot