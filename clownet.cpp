#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void CPP_calc_lambda(double& lambda, NumericMatrix& m, NumericVector& y){
  double nrow = m.rows();
  double ncol = m.cols();
  
  //double lambda;
  NumericVector lambdas_calc(ncol);
  double inner_prod;
  double max = 0;
  double pond = 1/nrow;
  
  for(int j = 1; j < ncol; j++){ //Sem levar o intercepto em conta
    for(int i = 0; i < nrow; i++){
      inner_prod = inner_prod + m(i,j)*y[i];
    }
    lambdas_calc[j] = pond*fabs(inner_prod);
    if(lambdas_calc[j] > max){
      max = lambdas_calc[j];
    }
    inner_prod = 0;
  }
  
  lambda = max;
}

// [[Rcpp::export]]
double CPP_veros_parcial(double par, 
                     NumericMatrix& X, NumericVector& w,
                     NumericVector& z, NumericVector B,
                     int j, double s)
{
  double lasso;
  double nrow = X.rows();
  double ncol = X.cols();
  double out;
  NumericVector soma(nrow);
// 
//   if(j != 0){
//     lasso = -s*fabs(par);
//   }
//   else {
//     lasso = 0;
//   }
//   
//   B[j] = lasso;

  for(int i = 0; i < ncol; i++){
    if(i == 0){
      for(int k = 0; k < nrow; k++){
        soma[k] = soma[k] + z[k] - B[0];
      }
    }
    else{
      for(int k = 0; k < nrow; k++){
        soma[k] = soma[k] + B[i]*X(k,j);
      }
    }
  }
  
  for(int k = 0; k < nrow; k++){
    out += pow(w[k]*soma[k],2);
  }
  
  return((-1/(2*nrow)*out));
}

// [[Rcpp::export]]
double calc_grad(NumericMatrix& X,
               NumericVector& z,
               NumericVector& w,
               NumericVector& B,
               NumericVector& eta,
               int j
               )
{
  double nrow = X.rows();
  double ncol = X.cols();
  NumericVector out(nrow);
  
  if(j == 0){
    for(int k = 0; k < nrow; k++){
      out[k] = out[k] + 2*w[k]*(z[k] - eta[k]);
    }
  }
  else{
    for(int k = 0; k < nrow; k++){
      out[k] = out[k] + 2*w[k]*(z[k] - eta[k])*X(k,j);
    }
  }
  
  return(sum(out)*(-1/(2*nrow)));
}
  
// [[Rcpp::export]]
double soft_threshold(double z, double s){
  if (z > s){
    return(z - s);
  } 
  else if (z < -s){
    return(z + s);  
  }
  else{
    return(0);
  } 
}
  
// [[Rcpp::export]]
void calc_eta(NumericMatrix& X,
              NumericVector& y,
              NumericVector& B,
              NumericVector& eta,
              NumericVector& mu,
              NumericVector& z,
              NumericVector& w
              ){
  int i,j,k;
  double nrow = X.rows();
  double ncol = X.cols();
  
  // for(i = 0; i < nrow; i++){
  //   eta[i] = 0;
  // }
  
  std::fill(eta.begin(), eta.end(), 0);
  
  for(j = 0; j < ncol; j++){ //new
    if(B[j] != 0){
      for(i = 0; i < nrow; i++){
        eta[i] = eta[i] + B[j]*X(i,j);
      }
    }
  }
  
  for(k = 0; k < nrow; k++){
    mu[k] = 1/(1 + exp(-eta[k]));
    w[k] = mu[k]*(1-mu[k]);
    z[k] = log(mu[k]/(1 - mu[k])) + (y[k] - mu[k])*(1/w[k]);
  }
}

// [[Rcpp::export]]
List CPP_clownet(
    NumericMatrix X, 
    NumericVector y, 
    NumericVector escala,
    NumericVector centro,
    NumericVector lambdas = NA_REAL,
    bool cv = false,
    int n_lambda = 50, 
    bool trace = true)
{
  
  double aux;
  double lambda;
  double nrow = X.rows();
  double ncol = X.cols();
  double lambda_dist;
  double tol = pow(10, -7);
  double erro_lasso;
  double s;
  double step_size = 1;
  double val;
  double lambda_min;
  int maxit = 10;
  int it;
  int it_lasso;
  
  NumericVector eta(nrow);
  NumericVector mu(nrow);
  NumericVector z(nrow);
  NumericVector w(nrow);
  NumericVector seq_lambda(n_lambda);
  NumericVector betas_atuais(ncol);
  NumericVector betas_antigos(ncol);
  NumericVector erros_lasso(ncol);
  NumericVector grad(ncol);
  NumericVector out(nrow);
  NumericMatrix betas(n_lambda, ncol);
    
  if(cv == false){
    CPP_calc_lambda(lambda, X, y);
    
    lambda_min = log((1/n_lambda)*lambda + 0.001);
    lambda_dist = (log(lambda) - lambda_min)/n_lambda;
    
    seq_lambda[0] = exp(lambda_min) + 0.0001; //n deixar ir para zero
    
    for(int i = 1; i < seq_lambda.length(); i++){
      seq_lambda[i] = exp(log(seq_lambda[i-1]) + lambda_dist);
    }
  }
  else{
    seq_lambda = lambdas;
  }
  
  //Intercepto
  betas_atuais[0] = log(mean(y)/(1- mean(y)));
  
  for(int j = (n_lambda - 1); j >= 0; j--){
      Rcpp::checkUserInterrupt();
    
      if(trace == true){
        Rcout << "Lambda: " << seq_lambda[j] << ";\n";
      }
      
      erro_lasso = 1;
      it_lasso = 1;

      s = seq_lambda[j];
      
      while( (erro_lasso > tol) & (it_lasso < maxit)){
        
        
        calc_eta(X, y, betas_atuais, eta, mu, z, w);
        betas_antigos = clone(betas_atuais);
        
        for(int k = 0; k < ncol; k++){
          
          grad[k] = calc_grad(X, z, w, betas_atuais, eta, k);
          
          aux = betas_atuais[k] - step_size*grad[k];

          if(k == 0){
            betas_atuais[k] = aux; 
          }
          else {
            betas_atuais[k] = soft_threshold(aux, s);
          }
 
          if(fabs(betas_atuais[k] - betas_antigos[k]) > tol){
            calc_eta(X, y, betas_atuais, eta, mu, z, w);
            erros_lasso[k] = fabs(betas_atuais[k] - betas_antigos[k]);
          }
        }
        
        //Depois de atualizar compute o máximo os erros
        erro_lasso = max(erros_lasso);
        
        if(erro_lasso < tol){
          break;
        }
        
        it_lasso = it_lasso + 1;
      }
      
      betas(j, _) = betas_atuais;
      
      //desescalar coeficientes
      // val = 0;
      // for(int l = 1; l < ncol; l++){
      //   val = val + betas_atuais[l]*(centro[l-1]/escala[l-1]);
      //   betas(j, l) = betas_atuais[l]/escala[l-1];
      // }
      // 
      // betas(j, 0) = betas_atuais[0] - val;
      
  }
  // End it all;
  return(List::create(Named("betas") = betas,
                      Named("lambdas") = seq_lambda));
}

// [[Rcpp::export]]
void CPP_calc_deviance(NumericVector& y,
                       NumericVector& mu,
                       NumericVector& desvio)
{
  int i;
  for(i = 0; i < y.length(); i++){ 
    desvio[i] = y[i]*log(mu[i]) + (1 - y[i])*log(1 - mu[i]);
  }
}
  
  
// [[Rcpp::export]]
NumericVector CPP_calc_deviance_vector(NumericMatrix M,
                                       NumericVector y)
{
  int i, j, k;
  double nrow = M.rows();
  double ncol = M.cols();
  double val;  
  NumericVector out(ncol);

  for(j = 0; j < ncol; j++){
    val = 0;
    for(i = 0; i < y.length(); i++){ 
      if(fabs(1 - M(i,j)) < 0.001) {
        M(i,j) = 1 - 0.001;
      }
      else if (M(i,j) < 0.001){
        M(i,j) = 0.001;
      }
      val = val + y[i]*log(M(i,j)) + (1 - y[i])*log(1 - M(i,j));
    }
    out[j] = -2*val;
  }
 
  return(out);
}

// [[Rcpp::export]]
double CPP_calc_lambda_double(NumericMatrix X, NumericVector y){
  double nrow = X.rows();
  double ncol = X.cols();
  
  //double lambda;
  NumericVector lambdas_calc(ncol);
  double inner_prod;
  double max = 0;
  double pond = 1/nrow;
  
  for(int j = 1; j < ncol; j++){ 
    for(int i = 0; i < nrow; i++){
      inner_prod = inner_prod + X(i,j)*y[i];
    }
    lambdas_calc[j] = pond*fabs(inner_prod);
    if(lambdas_calc[j] > max){
      max = lambdas_calc[j];
    }
    inner_prod = 0;
  }
  
  return(max);
}

// [[Rcpp::export]]
NumericMatrix CPP_predict(NumericMatrix X,
              NumericMatrix B,
              int n_lambdas,
              NumericVector y)
{
  double nrow = X.rows();
  double ncol = X.cols();
  int i, j, k;
  NumericMatrix out(nrow, n_lambdas);
  
  for(k = 0; k < n_lambdas; k++){
    for(i = 0; i < ncol; i++){
      for(j = 0; j < nrow; j++){
        out(j,k) = out(j,k) + B(k,i)*X(j,i);
      }
    }
  }
  
  //apply de exp
  for(k = 0; k < n_lambdas; k++){
    for(j = 0; j < nrow; j++){
      out(j,k) = 1/(1+ exp(-out(j,k)));
    }
  }
  
  return(out);
}

// // [[Rcpp::export]]
// List CPP_cv_clownet(
//     NumericMatrix X,
//     NumericVector y,
//     NumericVector escala,
//     NumericVector centro,
//     int n_blocos,
//     int n_lambda = 50,
//     bool trace = true)
// {
//   int i, j;
//   double nrow = X.rows();
//   double ncol = X.cols();
//   double acumulado_bloco;
//   double lambda;
//   double lambda_dist;
//   NumericVector seq_lambda(n_lambda);
//   NumericVector desvio(nrow);
//   double tamanho_bloco;
// 
//   CPP_calc_lambda(lambda, X, y);
// 
//   lambda_dist = lambda/n_lambda;
// 
//   for(int i = 1; i < seq_lambda.length(); i++){
//     seq_lambda[i] = seq_lambda[i-1] + lambda_dist;
//   }
// 
//   tamanho_bloco = nrow/n_blocos;
// 
//   for(j = 0; j < n_blocos; j++){
//     if(j == n_blocos){
//       NumericVector bloco_retirar = Range((j-1)*tamanho_bloco, nrow - 1);
//       NumericMatrix X_treino = X(-bloco_retirar, _)
//     }
//     else{
//       NumericVector bloco_retirar = Range((j-1)*tamanho_bloco, j*tamanho_bloco - 1);
// 
//     }
//   }
// }