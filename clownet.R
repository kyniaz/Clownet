clownet = function(X, y, n_lambdas = 50, trace = F, maxit = 10){
  prob_min = 10^(-5)
  prob_max = 1 - prob_min
  chamado = colnames(X)
  
  #Matriz
  if(any(colnames(X) %in% "(Intercept)")){
    X = X[, !colnames(X) %in% "(Intercept)"]
  }
  
  X = scale(X)
  attr_center = attr(X, "scaled:center")
  attr_scale = attr(X, "scaled:scale")
  X = as.matrix(cbind(Intercept = rep(1, nrow(X)), X))
  
  obj = CPP_clownet(X, y, attr_scale, attr_center, n_lambda = n_lambdas, trace = trace,
                    maxit)
  
  preditos = CPP_predict(X, obj$betas, n_lambdas + 1, y)
  preditos = pmin(pmax(preditos,prob_min),prob_max)

  lp = apply(preditos, 2, \(x) sum(-2*((1-y)*log(1-x) + y*log(x))))
  
  zeros = apply(obj$betas, 1, \(x) sum(x == 0))
  colnames(obj$betas) = colnames(X)
  
  out = list(chamado = chamado,
             lambdas = obj$lambdas, betas = obj$betas, preditos = preditos, dev.x = lp,
             centro_coef = attr_center, escala_coef = attr_scale, zeros = zeros)
  class(out) = c('list', 'clownet')
  return(out)
}

cv.clownet = function(X, y, n_lambdas = 50, n_blocos = 10, trace = FALSE,
                      maxit = 10){

  obj = clownet(X, y, n_lambdas = n_lambdas)
  lambdas = obj$lambdas
  
  n = nrow(X)
  
  compo_de_desvio = matrix(0, n_blocos, n_lambdas + 1)
  ep = trunc(n/n_blocos)
  preditos = matrix(NA, nrow = n, ncol = n_lambdas + 1)
  
  prob_min = 10^(-5)
  prob_max = 1 - prob_min
  
  #Matriz
  if(any(colnames(X) %in% "(Intercept)")){
    X = X[, !colnames(X) %in% "(Intercept)"]
  }
  
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
                                n_lambda = n_lambdas, trace = trace,
                                maxit = maxit) 
    
    X_bloco = X[bloco,]

    aux = CPP_predict(X_bloco, aux$betas, n_lambdas + 1, y[bloco])
    
    preditos[bloco,] = aux
    
    compo_de_desvio[i,] = CPP_calc_deviance_vector(aux, y[bloco])
  }
  
  preditos = preditos[-1,]
  preditos = pmin(pmax(preditos,prob_min),prob_max)
  
  lp = (apply(compo_de_desvio, 2, mean)/n)*n_blocos
  lu = lp + qchisq(0.975, 1)*apply(compo_de_desvio, 2, sd)/n
  ld = lp - qchisq(0.975, 1)*apply(compo_de_desvio, 2, sd)/n
  
  zeros = apply(obj$betas, 1, \(x) sum(x == 0))
  
  out = list(lambdas = lambdas, dev.x = lp, dev.up = lu, betas = obj$betas,
             dev.lo = ld, centro_coef = attr_center, escala_coef = attr_scale,
             zeros = zeros)
  class(out) = c('cv.clownet', 'clownet', 'list')
  return(out)
}

extract.clownet = function(obj, s){
  if(s == "lambda_min"){
    ls = which.min(obj$dev.x)
  }
  else if(s == 0){
    ls = which.min(obj$lambda)
  }
  else if(s > max(obj$lambda) | s < min(obj$lambda)){
    stop("choose a valid value for s.")
  }
  else if(is.numeric(x)){
    ls = which.min(abs(s - obj$lambdas))
  }
  else{
    stop("choose a valid value for s.")
  }
  
  betas = obj$betas[ls,]
  betas_out = numeric(length(betas))

  centro = obj$centro_coef
  escala = obj$escala_coef
  
  #browser()
  #desescalar coeficientes
  val = 0;
  for(l in 2:length(betas)){
    val = val + betas[l]*(centro[l-1]/escala[l-1]);
    betas_out[l] = betas[l]/escala[l-1];
  }
  
  betas_out[1] = betas[1] - val;
  names(betas_out) = names(betas)
  
  return(betas_out)
}

#Extrair coeficientes;
coef.clownet = function(obj, s){
  do.call("extract.clownet", list(obj, s))
}

plot.cv.clownet = function(obj){
  df_plot = data.frame(dev_min = obj$dev.lo,
                       dev_mean = obj$dev.x,
                       dev_max = obj$dev.up,
                       lambda = obj$lambdas)
  
  lambda_min = log(df_plot$lambda[which.min(df_plot$dev_mean)])
                
  ggplot(df_plot, aes(x = log(lambda), y = dev_mean)) + 
    geom_point(size = 1, col = "red") +
    geom_line(col = "red") + 
    geom_errorbar(aes(ymin=dev_min,ymax=dev_max), col = "black") +
    geom_segment(aes(x = lambda_min, xend = lambda_min,
                     y = min(dev_min), yend = max(dev_max)), col = "royalblue",
                 linetype = "dashed") +
    #geom_line(aes(x = log(lambda), y = dev_mean), col = "black", linetype = "dashed") +
    theme_minimal() +
    theme(plot.title=element_text(hjust=0.5)) +
    labs(x = expression(log(lambda)), y = "Deviance")
}


plot.clownet = function(obj){
  aux = data.frame(obj$betas)
  
  aux = aux[,-1]
  centro = obj$centro_coef
  escala = obj$escala_coef
  
  val = numeric(nrow(aux));
  for(l in 1:ncol(aux)){
    val = val + aux[,l]*(centro[l]/escala[l]);
    aux[,l] = aux[,l]/escala[l];
  }
  
  chamado = obj$chamado[-1]
  if(is.null(chamado)){
    chamado = paste0("V",1:ncol(chamado))
  }
  colnames(aux) = chamado
  
  aux$l1 = rowSums(abs(aux))
  
  plt = eval(parse(text = paste0('ggplot(aux) + ', 
                                 paste0('geom_line(aes(x = l1, y = ',
                                    chamado,',
                                    color = "',chamado,'"))',
                                 collapse = '+'))))
  
  plt = plt + theme_minimal() + scale_colour_brewer(palette = "Set1") +
  labs(x = "Soma L1", y = "Valor de Coeficientes")
  
  plt
}

deviance.clownet = function(obj) {
  return(obj$lp)
}

predict.clownet = function(obj, s, newdata){
  betas = coef.clownet(obj, s)
  return(newdata%*%betas)
}
