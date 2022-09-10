#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppClock)]]
#include <RcppClock.h>
#include <thread>

// [[Rcpp::export]]
void CPP_calc_lambda(double& lambda, const NumericMatrix& m, const IntegerVector& y)
//****************************************************************************
//
//  Purpose:
//
//    CPP_calc_lambda is a subroutine (void type) that calculates the maximum 
//    value of lambda in which all coefficients are 0.
//
//  Discussion:
//
//  Parameters:
//
//    Input/Output, double& lambda, value of lambda, to be used later, in CPP_clownet.
//
//    Input, NumericMatrix& (double **) m, design matrix of the regression.
//
//    Input, IntegerVector& (int *) y, numeric vector of response, expected to
//    be integer (0's and 1's).
//
{
  double nrow = m.rows();
  double ncol = m.cols();
  
  //double lambda;
  NumericVector lambdas_calc(ncol);
  double inner_prod;
  double max = 0;
  double pond = 1/nrow;
  
  for(int j = 1; j < ncol; j++){ //Sem levar o intercepto em conta
    for(int i = 0; i < nrow; i++){
      inner_prod += m(i,j)*y[i];
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
double CPP_partial_log_likehood( 
                     NumericMatrix& X, NumericVector& w,
                     NumericVector& z, NumericVector& B,
                     int j)
//****************************************************************************
//
//  Purpose:
//
//    CPP_partial_log_likehood calculates the partial log-likelihood
//    for one coefficient.
//
//  Discussion:
//
//    Could be used if we solved the LASSO problem by Gradient Descent.
//
//    Input/Output, double& lambda, value of lambda, to be used later, in CPP_clownet.
//
//    Input, NumericMatrix& (double **) m, design matrix of the regression.
//
//    Input, IntegerVector& (int *) y, numeric vector of response, expected to
//    be integer (0's and 1's).
//
{
  double nrow = X.rows();
  double ncol = X.cols();
  double out = 0;
  NumericVector soma(nrow);

  for(int i = 0; i < ncol; i++){
    if(i == 0){
      for(int k = 0; k < nrow; k++){
        soma[k] += z[k] - B[0];
      }
    }
    else{
      for(int k = 0; k < nrow; k++){
        soma[k] += B[i]*X(k,j);
      }
    }
  }
  
  for(int k = 0; k < nrow; k++){
    out += pow(w[k]*soma[k],2);
  }
  
  return((-1/(2*nrow)*out));
}

// [[Rcpp::export]]
double calc_grad(const NumericMatrix& X,
                 const NumericVector& z,
                 const NumericVector& w,
                 const NumericVector& B,
                 const NumericVector& eta,
                 int j
               )
//****************************************************************************
//
//  Purpose:
//
//    calc_grad is an internal C++ function to compute the gradient (or fist
//    deritivate) of the log likelihood in relation to one Bj.
//
//  Discussion:
//    
//    This derivative is symbolic and it's calculation are exact.
//
//  Parameters:
//
//    Input, NumericMatrix& (double **) X, design matrix of the regression.
//
//    Input, NumericVector& (double) w, weights (the same as Fisher's Scoring) for
//    the logistic regression with logit linkage.
//    
//    Input, NumericVector& (double) z, working response for
//    the logistic regression with logit linkage.
//
//    Input, NumericVector& (double) B, current estimates for B.
//
//    Input, NumericVector& (double) eta, linear predictor, or XB.
//
//    Input, int j, the coordinate index.
//
//    Output, double out, gradient for the coordinate Bj.
{
  double nrow = X.rows();
  double out = 0;
  
  if(j == 0){
    for(int k = 0; k < nrow; k++){
      out += 2*w[k]*(z[k] - eta[k]);
    }
  }
  else{
    for(int k = 0; k < nrow; k++){
      out += 2*w[k]*(z[k] - eta[k])*X(k,j);
    }
  }
  
  out = out*(-1/(2*nrow));
    
  return(out);
}
  
// [[Rcpp::export]]
double calc_hess(NumericMatrix& X,
                 NumericVector& z,
                 NumericVector& w,
                 NumericVector& B,
                 NumericVector& eta,
                 int j
)
//****************************************************************************
//
//  Purpose:
//
//    calc_hess is an internal C++ function to compute the hessian (or second
//    deritivate) of the log likelihood in relation to one Bj.
//
//  Discussion:
//    
//    This derivative is symbolic and it's calculation exact.
//
//  Parameters:
//
//    Input, NumericMatrix& (double **) X, design matrix of the regression.
//
//    Input, NumericVector& (double) w, weights (the same as Fisher's Scoring) for
//    the logistic regression with logit linkage.
//    
//    Input, NumericVector& (double) z, working response for
//    the logistic regression with logit linkage.
//
//    Input, NumericVector& (double) B, current estimates for B.
//
//    Input, NumericVector& (double) eta, linear predictor, or XB.
//
//    Input, int j, the coordinate index.
//
//    Output, double out, hessian (just second derivative) for the coordinate Bj.
{
  double nrow = X.rows();
  //double ncol = X.cols();
  double out = 0;
  
  if(j == 0){
    return(1);
  }
  else{
    for(int k = 0; k < nrow; k++){
      out += 2*w[k]*(pow(X(k,j),2));
    }
  }
  
  return(out*(1/(2*nrow)));
}

// [[Rcpp::export]]
double soft_threshold(double z, double s)
//****************************************************************************
//
//  Purpose:
//
//    soft_threshold is an internal C++ function to compute the soft thresholding
//    for the value z in respect to s.
//
//  Discussion:
//    
//
//  Parameters:
//
//    Input, double z, design matrix of the regression.
//
//    Input, double w, weights (the same as Fisher's Scoring) for
//    the logistic regression with logit linkage.
//
//    Output, double out, hessian (just second derivative) for the coordinate Bj.  
  {
  if (z > s){
    return(z - s);
  } 
  else if (z < (-s)){
    return(z + s);  
  }
  else{
    return(0.0);
  } 
}
 
 
// [[Rcpp::export]]
double exp_approx(double x) 
//****************************************************************************
//
//  Purpose:
//
//    exp_approx is C++ routine to calculate the exponential of X.
//
//  Discussion:
//    The std::exp and std::log functions can be somewhat slow compared to their
//    FORTRAN 77 counterparts (up to 2X slow in our systems).
//    Some programmers rewrite those functions to be less acurate but faster, this one
//    makes a Taylor approximation and it's around 2x faster.
//
//  Parameters:
//
//    Input, double x, the value for which compute exp(x).
//
//    Output, double x, exponential of x.   
{
 x = 1.0 + x / 1024;
 x *= x; x *= x; x *= x; x *= x;
 x *= x; x *= x; x *= x; x *= x;
 x *= x; x *= x;
 return x;
}

// [[Rcpp::export]]
void calc_eta(const NumericMatrix& X,
              const IntegerVector& y,
              const NumericVector& B,
              NumericVector& eta,
              NumericVector& mu,
              NumericVector& z,
              NumericVector& w
              )
//****************************************************************************
//
//  Purpose:
//
//    calc_eta is an internal C++ subroutine to update the linear predictor
//    and it's dependencies (mu,z and w).
//
//  Discussion:
//    
//    void functions are awesome.
//
//  Parameters:
//
//    Input, NumericMatrix& (double **) X, design matrix of the regression.
//
//    Input, IntegerVector& (int) y, response vector (0's and 1's).
//
//    Input, NumericVector& (double) B, current estimates for B.
//
//    Input, NumericVector& (double) eta, linear predictor, or XB.
//
//    Input, NumericVector& (double) mu, current estimated response.
//
//    Input, NumericVector& (double) w, weights (the same as Fisher's Scoring) for
//    the logistic regression with logit linkage.
//    
//    Input, NumericVector& (double) z, working response for
//    the logistic regression with logit linkage.
//
{
  int i,j,k;
  double nrow = X.rows();
  double ncol = X.cols();
  double prob_min = pow(10, -5);
  double prob_max = 1.0 - prob_min;
  
  std::fill(eta.begin(), eta.end(), 0);

  for(j = 0; j < ncol; j++){ 
    if(B[j] != 0){
      for(i = 0; i < nrow; i++){
        eta[i] += B[j]*X(i,j);
      }
    }
  }

  for(k = 0; k < nrow; k++){
    mu[k] = 1.0/(1.0 + exp_approx(-eta[k]));
    
    if(mu[k] < prob_min){
      mu[k] = prob_min;
    }
    else if(mu[k] > prob_max){
      mu[k] = prob_max;
    }
    
    w[k] = mu[k]*(1.0-mu[k]);
    z[k] = log(mu[k]/(1.0 - mu[k])) + (y[k] - mu[k])*(1.0/w[k]);
  }
}


// [[Rcpp::export]]
List CPP_clownet(
    const NumericMatrix& X, 
    const IntegerVector& y, 
    const NumericVector& escala,
    const NumericVector& centro,
    NumericVector lambdas = NA_REAL,
    bool cv = false,
    int n_lambda = 50, 
    bool trace = false,
    int maxit = 10)
//****************************************************************************
//
//  Purpose:
//
//    CPP_clownet is the workhorse function, it calls all others and return the result
//    for the R Code.
//
//  Discussion:
//    
//
//  Parameters:
//
//    Input, NumericMatrix& (double **) X, design matrix of the regression.
//
//    Input, IntegerVector& (int) y, response vector (0's and 1's).
//
//    Input, NumericVector& (double) escala, could be used if we wanted to deescale the coefficients.
//
//    Input, NumericVector& (double) centro, could be used if we wanted to decenter the coefficients.
//
//    Input, NumericVector& (double) lambdas, lambda sequence provided by the user
//
//    Input, bool cv, if true then the lambda sequence provided is followed.
//    
//    Input, int n_lambda, number of lambdas for the solution.
//
//    Input, bool trace, if true then the lambda is printed as the computation run.
//
//    Input, int maxit, number of lambdas for the solution.
//
//    Output, List containing two elements: the matrix of the coefficients B (not sparse in
//    storage) and the lambda sequence that generated it.
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
  int it_lasso;
  Rcpp::Clock clock;
  
  NumericVector eta(nrow);
  NumericVector mu(nrow);
  NumericVector z(nrow);
  NumericVector w(nrow);
  NumericVector seq_lambda(n_lambda + 1);
  NumericVector betas_atuais(ncol);
  NumericVector betas_antigos(ncol);
  NumericVector erros_lasso(ncol);
  NumericVector grad(ncol);
  NumericVector out(nrow);
  NumericMatrix betas(n_lambda + 1, ncol);
    

  if(cv == false){
    CPP_calc_lambda(lambda, X, y);
    
    lambda_min = log((1/n_lambda)*lambda);
    
    if(lambda_min == R_NegInf){ //C++ has lower precision than R
      lambda_min = -10;
    }
      
    lambda_dist = (log(lambda) - lambda_min)/n_lambda;
    
    seq_lambda[0] = exp(lambda_min); 
    
    for(int i = 1; i < seq_lambda.length(); i++){
      seq_lambda[i] = exp(log(seq_lambda[i-1]) + lambda_dist);
    }
    
    seq_lambda[n_lambda] = lambda; 
  }
  else{
    seq_lambda = lambdas;
  }
  
  //Intercept
  betas_atuais[0] = log(mean(y)/(1- mean(y)));


  for(int j = n_lambda; j >= 0; j--){
      Rcpp::checkUserInterrupt();
    
      if(trace == true){
        Rcout << "Lambda: " << seq_lambda[j] << ";\n";
      }
      
      erro_lasso = 1;
      it_lasso = 1;

      s = seq_lambda[j];
      
      while((erro_lasso > tol) & (it_lasso < maxit)){
        
        calc_eta(X, y, betas_atuais, eta, mu, z, w);

        betas_antigos = clone(betas_atuais);
        
        for(int k = 0; k < ncol; k++){
            val = step_size*calc_grad(X, z, w, betas_atuais, eta, k);
            aux = betas_atuais[k] - val; //calc_hess(X, z, w, betas_atuais, eta, k);

            if(k == 0){
              betas_atuais[k] = aux; 
            }
            else {
              betas_atuais[k] = soft_threshold(aux, s);
            }
            
            erros_lasso[k] = fabs(betas_atuais[k] - betas_antigos[k]);
        }
        
        erro_lasso = sum(erros_lasso);
        
        if(erro_lasso < tol){
          break;
        }
        
        it_lasso += 1;
      }
      

      betas(j, _) = betas_atuais;

  }
  // End it all;
  
  return(List::create(Named("betas") = betas,
                      Named("lambdas") = seq_lambda));
}

// [[Rcpp::export]]
//****************************************************************************
//
//  Purpose:
//
//    CPP_calc_deviance calculates the deviance for each observation.
//
//  Discussion:
//    
//
//  Parameters:
//
//    Input, IntegerVector& (int) y, response vector (0's and 1's).
//
//    Input, NumericMatrix& (double) mu, estimated vector of probabilities.
//
//    Output, desvio&, deviance for each observation.
//
void CPP_calc_deviance(IntegerVector& y,
                       NumericVector& mu,
                       NumericVector& desvio)
{
  int i;
  for(i = 0; i < y.length(); i++){ 
    desvio[i] = y[i]*log(mu[i]) + (1 - y[i])*log(1 - mu[i]);
  }
}
  
  
// [[Rcpp::export]]
//****************************************************************************
//
//  Purpose:
//
//    CPP_calc_deviance_vector calculates the deviance for each observation and sums it.
//
//  Discussion:
//    
//
//  Parameters:
//
//    Input, IntegerVector& (int) y, response vector (0's and 1's).
//
//    Input, NumericMatrix& (double) mu, estimated vector of probabilities.
//
//    Output, double out, the deviance for the entire model.
//
NumericVector CPP_calc_deviance_vector(NumericMatrix& M,
                                       const IntegerVector& y)
{
  int i, j;
  //double nrow = M.rows();
  double ncol = M.cols();
  double val;  
  NumericVector out(ncol);
  double e = pow(10, -5);
  
  for(j = 0; j < ncol; j++){
    val = 0;
    for(i = 0; i < y.length(); i++){ 
      if(fabs(1 - M(i,j)) < e) {
        M(i,j) = 1 - e;
      }
      else if (M(i,j) < e){
        M(i,j) = e;
      }
      val += y[i]*log(M(i,j)) + (1 - y[i])*log(1 - M(i,j));
    }
    out[j] = -2*val;
  }
 
  return(out);
}

// [[Rcpp::export]]
double CPP_calc_lambda_double(const NumericMatrix& X, const IntegerVector& y)
//****************************************************************************
//
//  Purpose:
//
//    CPP_calc_lambda_double is a function that calculates the maximum 
//    value of lambda in which all coefficients are 0.
//
//  Discussion:
//
//  Parameters:
//
//    Input, NumericMatrix& (double **) m, design matrix of the regression.
//
//    Input, IntegerVector& (int *) y, numeric vector of response, expected to
//    be integer (0's and 1's).
//
//    Output, double lambda, value of lambda, to be used later, in CPP_clownet.
//
{
  double nrow = X.rows();
  double ncol = X.cols();
  
  //double lambda;
  NumericVector lambdas_calc(ncol);
  double inner_prod;
  double max = 0;
  double pond = 1/nrow;
  
  for(int j = 1; j < ncol; j++){ 
    for(int i = 0; i < nrow; i++){
      inner_prod += inner_prod + X(i,j)*y[i];
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
NumericMatrix CPP_predict(const NumericMatrix& X,
                          const NumericMatrix& B,
                          const int& n_lambdas,
                          const IntegerVector& y)
//****************************************************************************
//
//  Purpose:
//
//    CPP_calc_lambda_double is a function that calculates the maximum 
//    value of lambda in which all coefficients are 0.
//
//  Discussion:
//
//  Parameters:
//
//    Input, NumericMatrix& (double **) X, design matrix of the regression.
//
//    Input, NumericVector& B, matrix of estimates, each column is for a lambda value.
//
//    Input, int& n_lambdas, number of lambdas used for B.
//
//    Input, IntegerVector& (int *) y, numeric vector of response, expected to
//    be integer (0's and 1's).
//
//    Output, NumericVector out, predicted values for each lambda (for deviance calculation in CV).
//
{
  double nrow = X.rows();
  double ncol = X.cols();
  int i, j, k;
  NumericMatrix out(nrow, n_lambdas);
  
  for(k = 0; k < n_lambdas; k++){
    for(i = 0; i < ncol; i++){
      for(j = 0; j < nrow; j++){
        out(j,k) += B(k,i)*X(j,i);
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
