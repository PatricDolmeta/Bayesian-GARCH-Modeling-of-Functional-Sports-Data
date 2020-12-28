#include <RcppArmadillo.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <list>			   
#include<stdio.h> 
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// Update of the GARCH variance according to equation (7)

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_h(arma::mat zeta, arma::rowvec alpha, double GARCHbeta, int n, arma::vec Si, int max_s){
	arma::mat h_new(n, max_s);
	h_new.zeros();
	for (int i = 0; i < n; i++){ 
		h_new(i,0) = alpha[0];
		for (int s=1; s<Si[i]; s++){
			h_new(i,s) = alpha[0] + alpha[1] * pow(zeta(i, s-1),2) + GARCHbeta * h_new(i,s-1); 
		}			  
	}
	return h_new;
}

// Adaptive scaling within the Adaptive Metropolis Hastings algorithm for alpha GARCH coefficient 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void mh_alpha(arma::rowvec alpha, double GARCHbeta, arma::cube& adapt_cov, arma::mat& Gamma_adapt,
			  double& Sd, arma::mat h, arma::mat zeta, arma::mat mu, double m,  
			  arma::mat sigma0alpha, arma::vec nu0alpha, arma::colvec Si, int max_s, int n, int nn,
			  arma::vec &alpha_new, arma::vec & alpha_p, list<double>& lambda) {
  
  //transform input parameter to be in the real plane abort
  
  arma::vec Gamma = arma::log(alpha.t());
  
  //covariance matrix given by previous iterations
	 	 
  arma::mat Cov_alpha = adapt_cov.slice(nn-1);			   
  
  //symmetric and adaptive propsal distribution 
  arma::vec Gamma_prop =  mvnrnd(Gamma, Cov_alpha);
  
  //original parameter proposal
  arma::vec alpha_prop(2);
  alpha_prop = arma::exp(Gamma_prop); 
  alpha_p = alpha_prop;
  
  arma::mat h_new = update_h(zeta, alpha_prop.t(), GARCHbeta, n, Si, max_s);
  
  //acceptance probability computations
  double lnlambda = 0;
  lambda.push_back(lnlambda);
  
  for (int i = 0; i < n; i++){
	for(int s = 0; s<Si[i]; s++){
	  lnlambda -= 0.5 * log(2 * h_new(i,s));
	  lnlambda -= pow((mu(i,s)-m),2)/(2*h_new(i,s)) ;
	  lnlambda += 0.5 * log(2 * h(i,s));
	  lnlambda += pow((mu(i,s)-m),2)/(2*h(i,s)) ;
	}
	lambda.push_back(lnlambda);
  }
  lnlambda -= 0.5 * as_scalar((alpha_prop - nu0alpha).t() * inv(sigma0alpha) * (alpha_prop - nu0alpha));
  lambda.push_back(lnlambda);
  lnlambda += Gamma_prop[0] + Gamma_prop[1];
  lambda.push_back(lnlambda);
  lnlambda += 0.5 * as_scalar((alpha.t() - nu0alpha).t() * inv(sigma0alpha) * (alpha.t() - nu0alpha));
  lambda.push_back(lnlambda);
  lnlambda -= Gamma[0] + Gamma[1];
  lambda.push_back(lnlambda);
  lnlambda = std::min(lnlambda, 0.0);
		
  double lnu=std::log(R::runif(0.0,1.0));
	  
  if (lnu<lnlambda){
	  alpha_new = alpha_prop;
	  Gamma_adapt.row(nn) = Gamma_prop.t();
						  
	
  } else {
	  alpha_new = alpha.t();
	  Gamma_adapt.row(nn) = Gamma.t();
  }
  
  double w_nn = std::pow(nn-1,-0.7);
  Sd = Sd * std::exp(w_nn * (std::exp(lnlambda) - 0.234));
  
  if(Sd<std::pow(10,-3)){
		Sd=std::pow(10,-3);
  } else if (Sd>std::pow(10,2)){
		Sd=std::pow(10,2);
  }
   
  if(nn < 500){
	adapt_cov.slice(nn) = 0.1 * arma::eye(2,2);
  }else{
	arma::mat aux1(2,2);
	aux1.zeros();
	arma::vec aux2(2);
	aux2.zeros();
	for (int i = 1; i < nn; i++){ 
		aux1 += Gamma_adapt.row(i).t() * (Gamma_adapt.row(i));
		aux2 += Gamma_adapt.row(i).t();
	}
	adapt_cov.slice(nn) = Sd/(nn-2) * (aux1 - 1.0/(nn-1) * aux2 * aux2.t()) + Sd * 0.1 * arma::eye(2,2);  
  }
} 
  
// Adaptive scale Metropolis algorithm for beta GARCH coefficient
 
  
void mh_beta(double GARCHbeta, arma::rowvec alpha, arma::mat h, 
			 arma::mat zeta, arma::mat mu, double m, 
			 double sigma0beta, double nu0beta, int n, 
			 arma::colvec Si, int max_s, int nn,
			 double& beta_new, double& beta_p, 
			 list<double>& lambda, arma::vec& adapt_var){
  
  double theta = log(GARCHbeta);
  double var_beta = adapt_var[nn];
  double theta_prop = theta + R::rnorm(0, pow(var_beta, 1/2));
  
  double lnlambda = 0;
  lambda.push_back(lnlambda);
  
  double beta_prop = std::exp(theta_prop);
  beta_p = beta_prop;  
  arma::mat h_new = update_h(zeta, alpha, beta_prop, n, Si, max_s);  
	  
  for (int i = 0; i < n; i++){ 
	for(int s = 0; s<Si[i]; s++){
	  lnlambda -= 0.5*log(2 * 2 * acos(0.0) * h_new(i,s));
      lnlambda -= pow((mu(i,s)-m),2)/(2*h_new(i,s)) ;
	  lnlambda += 0.5*log(2 * 2 * acos(0.0) * h(i,s));
	  lnlambda += pow((mu(i,s)-m),2)/(2*h(i,s)) ;
	}
	lambda.push_back(lnlambda);
  }
  lnlambda -= 0.5*pow( (beta_prop-nu0beta)/sigma0beta, 2);
  lambda.push_back(lnlambda);
  lnlambda += theta_prop;
  lambda.push_back(lnlambda);
  lnlambda += 0.5*pow( (GARCHbeta-nu0beta)/sigma0beta, 2);
  lambda.push_back(lnlambda);
  lnlambda -= theta;
  lambda.push_back(lnlambda);
  lnlambda = std::min(lnlambda, 0.0);
				
  double lnu=std::log(R::runif(0.0,1.0));
	  
  if (lnu<lnlambda){
	  beta_new = beta_prop;
  } else {
	  beta_new = GARCHbeta;
  }
  
  double w_nn = std::pow(nn,-0.7);
  double var_beta_new = var_beta * std::exp(w_nn * (std::exp(lnlambda) - 0.234));
  
  if(var_beta_new<std::pow(10,-5)){
		var_beta_new=std::pow(10,-5);
  } else if (var_beta_new>std::pow(10,5)){
		var_beta_new=std::pow(10,5);
  }
   
  adapt_var[nn+1] = var_beta_new;
  
} 
 
// GIBBS sampler

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List gibbs_compl_garch(arma::mat& data, arma::mat& O, arma::mat& B, 
								arma::mat tt,arma::colvec Si,
								int nit, int burn, int thin,
								double as, double bs, double df, double ad1, double bd1, 
								double ad2, double bd2, double aw, double bw,
								double nu0alpha1, double sigma0alpha1, double nu0alpha2, 
								double sigma0alpha2, double corr0alpha,
								double nu0beta, double sigma0beta, double m0, double s0, 
								double nu0psi, double sigma0psi,
								double prop, double b0, double b1, double epsilon, int k) 
{
	
  // open txt files for printed iterations							   
  ofstream th;
  th.open("Theta.txt");
  ofstream bet;
  bet.open("ExternalBeta.txt");
  ofstream ps;
  ps.open("Psi.txt");
  ofstream lambda_norm;
  lambda_norm.open("Lambda_Norm.txt");
  ofstream kkk;
  kkk.open("k.txt");
  ofstream mu_print;
  mu_print.open("Mu.txt");
  ofstream Lambda_tot;
  Lambda_tot.open("Lambda.txt");
  ofstream eta;
  eta.open("Eta.txt");
  ofstream m_print;
  m_print.open("m.txt");  
  ofstream h_is;
  h_is.open("h.txt");
  ofstream z_is;
  z_is.open("z.txt");
  ofstream alpha_print;
  alpha_print.open("alpha.txt");
  ofstream beta_print;
  beta_print.open("beta.txt");
  ofstream acc_prob_alpha;
  acc_prob_alpha.open("acc_prob_alpha.txt");
  ofstream acc_prob_beta;
  acc_prob_beta.open("acc_prob_beta.txt");
  ofstream alpha_proposal;
  alpha_proposal.open("alpha_prop.txt");
  ofstream beta_proposal;
  beta_proposal.open("beta_prop.txt");
  ofstream residuals;
  residuals.open("Residuals.txt");  
  
  arma::arma_rng::set_seed(230518);  

  // DATA READING
  
  arma::vec id = data.col(0);														// ID
  arma::colvec tij = data.col(1)/max(data.col(1)); 									// t_ij between 0 and 1, may already be in this form
  arma::colvec y = data.col(2); 													// results
  arma::colvec season = data.col(3);												// Season numbers for every performance 
  int max_s = max(season);															// max_s: longest career in seasons
								
  int N = y.n_elem; 																// N: total number of observation 
  arma::colvec n1 = unique(id);  													// n1: unique identifiers of athletes
  int n = n1.size();																// n: number of unique athletes
  arma::colvec ni = sum(tt,1);														// ni: number of observations for every athlete
  arma::colvec NI = cumsum(ni);														// NI: total observations up to athlete i 
  
  int r = O.n_cols;   									   
																				  
																	   
  int q = B.n_cols;   														// X_t: external covariates matrix

    // Data Centering
  y.subvec(0,sum(ni.subvec(0,0))-1)=y.subvec(0,sum(ni.subvec(0,0))-1)-mean(y.subvec(0,sum(ni.subvec(0,0))-1));
  for (int h=0; h<n-1; h++) {
    y.subvec(sum(ni.subvec(0,h)),sum(ni.subvec(0,h+1))-1)=y.subvec(sum(ni.subvec(0,h)),sum(ni.subvec(0,h+1))-1)-mean(y.subvec(sum(ni.subvec(0,h)),sum(ni.subvec(0,h+1))-1));
  }
  
  // Individual components in the mddel
  arma::colvec y1(N);
  arma::colvec y2(N);
  arma::colvec y3(N);
  arma::colvec y4(N);
  
  // DEFINITIONS
  
  // Functional Part
  

  arma::colvec sig(q);
  arma::mat Sigma(q,q);
  arma::mat Sigmainv(q,q);
  arma::mat Lambda(q,k);
  arma::mat InvEta(k,k);
  arma::mat eta2(n,k);
  arma::mat THETA2(q,n);
  arma::mat t(q, k);
  arma::colvec delta(k);
  arma::colvec tau(k);
  arma::mat Ptht(q,k);
  double ad,bd;
  arma::mat THETAtil(q,n);
  arma::colvec Thetatil2(THETAtil.n_rows);
  arma::colvec yy=y;
  double totsum=0;
  arma::vec temp1(q);
  
  //Seasonal component

  arma::vec m(nit);
  arma::mat h(n, max_s);
  h.zeros();
  arma::mat alpha(nit,2);
  arma::vec GARCHbeta(nit);
  arma::mat zeta(n, max_s);
  zeta.zeros();
  vector <arma::rowvec> g;
  arma::mat mu_mat(n,max_s);
  mu_mat.zeros();
  arma::colvec col_mean_mu(n);  
  arma::mat y_sum(n,max_s);
  double meanm=0,varm=0;
  double meanmu=0,varmu=0;
  double sumSi=sum(Si);
  double sumh=0;
  double sum_mu=0;
  int aux1=0,aux2=0;   
  int accepted_alpha = 0;
  int accepted_beta = 0;
  arma::vec adapt_var_beta(nit);
  adapt_var_beta.zeros();
  arma::cube adapt_Cov_alpha(2,2,nit);
  adapt_Cov_alpha.zeros();
  double Sd = 1.0;
  arma::mat Gamma_adapt(nit,2);
  Gamma_adapt.zeros();  
  
  
  // External regression part
  arma::mat XX_regr(r,r);
  arma::mat Xy_regr(r,1);
  arma::mat beta(nit,r);
  arma::mat Sigma_regr(r,r);
  arma::mat Sigmainv_regr(r,r);
  arma::colvec mu_regr(r);
  arma::colvec res_regr(N);
  double ssr_regr;
  arma::colvec beta0_regr(r);
  arma::mat Sigma0_regr(r,r);
  arma::vec acc(n);
  acc.zeros();
  
  // Error
  arma::colvec psi(nit);
  
  double nunpsi=0,sigmanpsi=0;
  arma::colvec psi2(nit);
    
  // Prediction
  
  arma::mat y_pred(nit,N);
  
  cout << " ok definitions" << endl;
  
  // INITIALIZATION
  
  // Functional Part
  
  for (int i=0; i<q; i++ ) {
    sig(i)=1/as_scalar(arma::randg(1, arma::distr_param(as,1/bs))); 
  }
  Sigma=arma::diagmat(sig);
  Sigma=(Sigma+Sigma.t())/2;
  Sigmainv=arma::inv(Sigma);
  Sigmainv=(Sigmainv+Sigmainv.t())/2;
  Lambda=arma::zeros(q,k); 
  InvEta = arma::eye(k,k);
							  
  eta2 =  arma::mvnrnd(arma::zeros(1,k).t(), InvEta, n).t();
	
  
  arma::mat aux22=Lambda*eta2.t();
  THETA2 = arma::mvnrnd((arma::zeros(1,Lambda.n_rows)).t(), Sigma, eta2.n_rows );  
  for (int j=0; j<eta2.n_rows; j++) {
    THETA2.col(j) = arma::mvnrnd(aux22.col(j),Sigma);
  }
  t = arma::randg( q, k, arma::distr_param(df/2,2/df) );
  delta(0)=as_scalar(arma::randg( 1, arma::distr_param(ad1,bd1)));
  delta.subvec(1,k-1)=arma::randg( k-1, 1, arma::distr_param(ad2,bd2));
  tau=arma::cumprod(delta);
  arma::mat aux_mat=(arma::repmat(tau.t(),q,1));
  Ptht=t%aux_mat;
  
  // Seasonal component part
  
  arma::vec nu0alpha = {nu0alpha1, nu0alpha2};
  arma::mat sigma0alpha = {{sigma0alpha1, corr0alpha},
						   {corr0alpha, sigma0alpha2}};
  sigma0alpha = sigma0alpha;
  arma::vec start_alpha; 
  start_alpha = arma::mvnrnd(nu0alpha, sigma0alpha);
  while(start_alpha[0] <= 0 || start_alpha[1] <= 0){
	  start_alpha = arma::mvnrnd(nu0alpha, sigma0alpha);
  }  
  alpha.row(0) = start_alpha.t();
  double u = R::runif(0,1);
  GARCHbeta[0]=R::qnorm((u+1)/2, nu0beta, sqrt(sigma0beta), 1, 0);
  m[0]=as_scalar(arma::randn(1)*sqrt(s0)+m0);
  
  // Initialization of adaptive variances for AMH
  
  adapt_var_beta[0] = 1.0;
  adapt_Cov_alpha.slice(0) = 0.1 * arma::eye(2,2);
  Gamma_adapt.row(0) = arma::log(start_alpha.t());
  
  // Initialization of h_is, z_is and mu_is
  
  for(int i = 0; i < n; i++){
	h(i,0) = alpha(0,0);
	mu_mat(i,0)=as_scalar(arma::randn(1) * sqrt(h(i,0)) + m[0]);
	zeta(i,0) = mu_mat(i,0) - m[0];
	for(int s = 1; s < Si[i]; s++){
		h(i,s) = alpha(0,0) + alpha(0,1) * pow(zeta(i, s-1),2) + GARCHbeta[0] * h(i,s-1);
		mu_mat(i,s)=as_scalar(arma::randn(1) * sqrt(h(i,s)) + m[0]);
		zeta(i,s) = mu_mat(i,s) - m[0];		
	}		
  }
 
  // g_is: number of shotput per season
  for(int i = 0; i < n; i++) {  														
    arma::rowvec temp = tt.submat(i,0,i,Si[i]-1);
    g.push_back(temp);
  }										  

  // Regressive component 
  beta0_regr.zeros();
  Sigma0_regr=arma::eye(r,r);
  beta.row(0)=arma::mvnrnd(beta0_regr,Sigma0_regr).t();

  
  // Error
  
  psi[0]=1/as_scalar(arma::randg(1, arma::distr_param(nu0psi/2,2/(nu0psi*sigma0psi))));
  psi2[0]=1/as_scalar(arma::randg(1, arma::distr_param(nu0psi/2,2/(nu0psi*sigma0psi))));
  
  // Additive components initialization 
  arma::mat x4 = B.submat(0,0, sum(ni.subvec(0,0))-1,B.n_cols-1);
  y1.subvec(0,sum(ni.subvec(0,0))-1)=x4*THETA2.col(0);
  for (int h=0; h<n-1; h++){
    arma::mat x4 = B.submat(sum(ni.subvec(0,h)),0, sum(ni.subvec(0,(h+1)))-1,B.n_cols-1);
    y1.subvec(sum(ni.subvec(0,h)),sum(ni.subvec(0,h+1))-1)=x4*THETA2.col(h);
  }
  
  for (int i=0; i<g.size();i++) {
    for (int s=0; s<g[i].size(); s++) {      
      for (int j=0; j<g[i][s]; j++) {
        y2[aux1]=mu_mat(i,s);
        aux1++;
      }      
    }
  }
  aux1=0;
  
  y3 = O*beta.row(0).t();		

  cout << " ok initializations" << endl;
  
  
    //GIBBS
  for (int nn=1; nn<nit; nn++) {
    
	// Functional part update
	
    y1=y-y2-y3;
	
	// Update Lammda
    arma::mat Lambda(q,k);
    for (int j=0; j<Lambda.n_rows; j++) {
      arma::mat Vlam1 = arma::diagmat(Ptht.row(j))+sig(j)*(eta2.t()*eta2);
      arma::mat Vlam=arma::inv(Vlam1);
      Vlam=(Vlam+Vlam.t())/2;
      arma::colvec Elam=Vlam*sig(j)*eta2.t()*THETA2.row(j).t();
      Lambda.row(j)=arma::mvnrnd(Elam,Vlam).t();
    }
    k=Lambda.n_cols;
	
	// Update Shinkage		
    for (int j=0; j<Lambda.n_cols; j++){
      for (int i=0; i<Lambda.n_rows; i++) {
        t(i,j) = arma::randg(arma::distr_param(df/2+0.5,1/(df/2+0.5*Lambda(i,j)*Lambda(i,j)*tau(j))));
      }
    }
	// Update Shinkage
    totsum=sum(tau.t()%sum(t%(Lambda%Lambda)));
    ad = ad1 + q*k/2;
    bd = bd1 + 0.5*(1/delta(0))*totsum;
    delta(0)=arma::randg(arma::distr_param(ad,1/bd));
    
	// Update Shinkage
	tau=arma::cumprod(delta);													   
    for (int h=0; h<k; h++){
      ad = ad2 + q*(k-h)/2;
      temp1=(tau.t()%sum(t%(Lambda%Lambda))).t();
      bd=bd2+0.5*(1/delta(h))*sum(temp1.subvec(h,k-1));
      delta(h)=arma::randg(arma::distr_param(ad,1/bd));
      tau=arma::cumprod(delta);
    }
	
	// Update diagonal elemnts of matrix D
    for (int i=0;i<Ptht.n_rows;i++)
      for (int j=0; j<Ptht.n_cols;j++)
        Ptht(i,j)=t(i,j)*tau(j);
    
    // Update Sigma          
    THETAtil=THETA2-Lambda*eta2.t();
    Thetatil2=sum(THETAtil%THETAtil,1);
    
    for (int i=0; i<q; i++) {
      sig(i)=1/as_scalar(arma::randg(1,arma::distr_param(as+n/2,1/(bs+0.5*Thetatil2(i)))));
    }
    Sigma=arma::diagmat(sig);
    Sigma=(Sigma+Sigma.t())/2;
    Sigmainv=arma::inv(Sigma);
    Sigmainv=(Sigmainv+Sigmainv.t())/2;
	 
    // Update Eta
    arma::mat x = B.submat(0,0,sum(ni.subvec(0,0))-1,B.n_cols-1);
    arma::mat inner = ((1/psi[nn-1])*arma::eye(ni(0),ni(0)) + x*Sigma*(x.t()));
    inner=(inner + inner.t())/2;
	arma::mat invinner = arma::inv(inner);
    invinner = (invinner + invinner.t())/2;
    arma::mat etavar = (Lambda.t())*(x.t())*(invinner)*x*Lambda + arma::eye(k,k);
    etavar = (etavar + etavar.t())/2;
	arma::mat invetavar = arma::inv(etavar);
    invetavar = (invetavar + invetavar.t())/2;
	arma::colvec meaneta = (Lambda.t())*(x.t())*invinner*(y1.subvec(0,sum(ni.subvec(0,(0)))-1));           
    eta2.row(0) = arma::mvnrnd(invetavar*meaneta, invetavar).t();
    for (int h =0; h< n-1; h++) {
      arma::mat x = B.submat(sum(ni.subvec(0,h)),0, sum(ni.subvec(0,(h+1)))-1,B.n_cols-1);
      arma::mat inner = ((1/psi[nn-1])*arma::eye(ni(h+1),ni(h+1)) + x*Sigma*(x.t()));
      inner = (inner + inner.t())/2;
      arma::mat invinner = arma::inv(inner);
      invinner = (invinner + invinner.t()) / 2;
      arma::mat etavar = (Lambda.t())*(x.t())*(invinner)*x*Lambda + arma::eye(k,k) ;
      etavar = (etavar + etavar.t())/2;
      arma::mat invetavar = arma::inv(etavar);
      invetavar = (invetavar + invetavar.t())/2;
      arma::colvec meaneta = (Lambda.t())*(x.t())*invinner*(y1.subvec(sum(ni.subvec(0,h)),sum(ni.subvec(0,(h+1)))-1));
      eta2.row(h+1) = arma::mvnrnd(invetavar*meaneta, invetavar).t();
    }
	
	
    // Update Theta
    arma::mat x2 = B.submat(0,0, sum(ni.subvec(0,0))-1,B.n_cols-1);
    arma:: mat covfun = Sigmainv + psi[nn-1]*(x2.t())*x2;
    covfun = (covfun + covfun.t()) / 2;
    arma::mat Invcovfun = arma::inv(covfun);
    Invcovfun = (Invcovfun + Invcovfun.t()) / 2;
    arma::mat meanvec= (Invcovfun*(Sigmainv*(Lambda*eta2.row(0).t()) + psi[nn-1]*(x2.t())*(y1.subvec(0,sum(ni.subvec(0,0))-1)))).t();
    THETA2.col(0) = arma::mvnrnd(meanvec.t(), Invcovfun);
    for (int h=0; h<n-1; h++) {
      arma::mat x2 = B.submat(sum(ni.subvec(0,h)),0, sum(ni.subvec(0,(h+1)))-1,B.n_cols-1);
      arma:: mat covfun = Sigmainv + psi[nn-1]*(x2.t())*x2;
      covfun = (covfun + covfun.t()) / 2;
      arma::mat Invcovfun = arma::inv(covfun);
      Invcovfun = (Invcovfun + Invcovfun.t()) / 2;
      arma::mat meanvec = (Invcovfun*(Sigmainv*(Lambda*eta2.row(h+1).t()) + psi[nn-1]*(x2.t())*(y1.subvec(sum(ni.subvec(0,h)),sum(ni.subvec(0,h+1))-1)))).t();
      THETA2.col(h+1) = arma::mvnrnd(meanvec.t(), Invcovfun);
    }
	
	// Update y1
	
	IntegerVector ni_cum(n+1);
    ni_cum[0]=0;
    for(int i1=1; i1<(n+1); i1++){
      ni_cum[i1]=ni_cum[i1-1]+ni[i1-1];
    }
    
    for(int h=0; h<(n);h++){
      for(int i=ni_cum[h];i<(ni_cum[h+1]);i++){
        y1[i]=dot(B.row(i),THETA2.col(h)) ;
      } 
    }	
	
    double prob;
    double uu;
    NumericVector lind(Lambda.n_cols);
    NumericVector vec(Lambda.n_cols);
    double num;
    
    // probability of adapting
    prob = 1/exp(b0 + b1*nn);               
    uu=as_scalar(arma::randu(1));

    // proportion of elements in each column less than eps in magnitude
    for(int ll=0; ll<Lambda.n_cols; ll++){
      lind[ll] = sum(abs(Lambda.col(ll)) < epsilon)/q  ;   
    }
    
    // location and number of redundant columns
    vec = lind >=prop; num = sum(vec);

    if( uu < min(prob,0.001*(nn>=5000)+(nn<5000) ) ){
      //cout << "adaptive k" << k << endl;
		if(nn > 20 & num == 0 & (sum(lind < 0.995)==(lind.size())) ){
			
			k = k+1;
		  
			Lambda.insert_cols(Lambda.n_cols,arma::zeros( q ));
			eta2.insert_cols(eta2.n_cols,arma::randn(n));
			t.insert_cols(t.n_cols, arma::randg(q,arma::distr_param(df/2,df/2)));
			delta.insert_rows(delta.n_rows,arma::randg(1,arma::distr_param(ad2,1/bd2)));
			tau=arma::cumprod(delta);
			arma::mat Ptht_new(q,k);
			for (int i = 0; i < q; i++){
				for (int j = 0; j < k; j++){
				Ptht_new(i,j) = t(i,j)*tau(j);
				}
			}
			Ptht = Ptht_new;
			//cout << "middle k" << k << endl;
			
		} 
		else if(num>0){
			int k_tmp;
			k_tmp = (k - num);
            if(k_tmp < 1){k_tmp=1;}

            arma::mat Lambda_temp(q,k_tmp);
            arma::mat t_temp(q, k_tmp);
            arma::mat eta2_temp(n,k_tmp);
            arma::colvec delta_temp(k_tmp);
			
            double cnt_temp = 0;
            for(int kk = 0; kk < k; kk++){
				if(vec[kk]==0){
                Lambda_temp.col(cnt_temp) = Lambda.col(kk);
                t_temp.col(cnt_temp) = t.col(kk);
                eta2_temp.col(cnt_temp) = eta2.col(kk);
                delta_temp(cnt_temp) = delta(kk);
                cnt_temp = cnt_temp+1;
				}
            }
			
			k = k_tmp;
			
            Lambda = Lambda_temp;
            t = t_temp;
            eta2 = eta2_temp;
            delta = delta_temp;
            tau = arma::cumprod(delta);
            arma::mat Ptht_temp(q,k);
            for (int i = 0;i < q;i++){
				for (int j = 0; j < k;j++){
					Ptht_temp(i,j) =  t(i,j)*tau(j);
                }
            }
			Ptht = Ptht_temp;
        }
		//cout << "updated k" << k << endl;
	}
   		   
    // Seasonal component update
	
    y2=y-y1-y3; 
	
	// computation sum over observations in each season
    for (int i=0; i<n;i++) {
      for (int s=0; s<Si[i]; s++) {      
        aux2=aux2+int(g[i][s]);
        if (aux2==aux1) y_sum(i,s)=0;
        else y_sum(i,s)=sum(y2.subvec(aux1,(aux2-1)));
        aux1=aux2;
      }
    } 
    aux1=0,aux2=0;
	
	// Update mu(i,s) GARCH
    for (int i=0; i<n; i++) {															
      for (int s=0; s<Si[i]; s++) {
	    meanmu=(y_sum(i,s)/(psi[nn-1])+m[nn-1]/h(i,s));
        varmu=(g[i][s]/(psi[nn-1])+1/h(i,s));
        varmu=pow(varmu,-1);
        mu_mat(i,s)=as_scalar(arma::randn(1)*sqrt(varmu)+(varmu*meanmu));
      }
      
    }
    
    // Update m		 
	sumh = 0;
	sum_mu = 0;		
    for (int i=0; i<n; i++){
      for (int s=0; s<Si[i]; s++){
          sumh+=(1./h(i,s));
		  sum_mu += mu_mat(i,s)/h(i,s);
	  }				   
    }
    
	varm = 1/(sumh + 1/s0);
	meanm = varm * (sum_mu + m0/s0);
    
    m[nn]=as_scalar(arma::randn(1)*sqrt(varm)+(meanm));

	//Update zeta
	for (int i = 0; i < n; i++){ 
		for (int s=0; s<Si[i]; s++){
			zeta(i,s) = mu_mat(i,s)-m[nn]; 
		}			  
	}
	
	//Update alpha with MH
	arma::vec alpha_new(2);
	alpha_new.zeros();
	
	arma::vec alpha_p(2);
	alpha_p.zeros();
	
	list<double> lambda;
	
	mh_alpha(alpha.row(nn-1), GARCHbeta[nn-1],adapt_Cov_alpha, Gamma_adapt,
           Sd, h, zeta, mu_mat, m[nn], sigma0alpha, nu0alpha, 
           Si, max_s, n, nn, alpha_new, alpha_p, lambda);
	
	alpha_proposal << alpha_p.t() <<endl;	
	acc_prob_alpha << "acc prob is" ; 
	
	if (alpha_new[0] != alpha.row(nn-1)[0] || alpha_new[1] != alpha.row(nn-1)[1]){
	  alpha.row(nn) = alpha_new.t();
	  accepted_alpha += 1;
	} else {
	  alpha.row(nn) = alpha.row(nn-1);
	}
	
	for( const auto& v : lambda ) acc_prob_alpha << v << " " ;
	acc_prob_alpha << " end" << endl; 

	//Update GARCHbeta with MH
	double GARCHbeta_new;
	double beta_p;
	
	list<double> lambda2;
	
	 mh_beta(GARCHbeta[nn-1], alpha.row(nn), update_h(zeta, alpha.row(nn), GARCHbeta[nn-1], n, Si, max_s),
		  zeta, mu_mat, m[nn], sigma0beta, nu0beta, n, Si, max_s, nn,
		  GARCHbeta_new, beta_p, lambda2, adapt_var_beta);
	
	beta_proposal << beta_p <<endl;	

	acc_prob_beta << "acc prob is" ; 
	
	if (GARCHbeta_new != GARCHbeta[nn-1]) {
	  GARCHbeta[nn] = GARCHbeta_new;
	  accepted_beta += 1;
	} else {
	  GARCHbeta[nn] = GARCHbeta[nn-1];
	}
	
	for( const auto& v : lambda2 ) acc_prob_beta << v << " " ;
	acc_prob_beta << " end" << endl; 

	//Update h_is
	h = update_h(zeta, alpha.row(nn), GARCHbeta[nn], n, Si, max_s);
    
    
    // Update y2
    for (int i=0; i<g.size();i++) {
      for (int s=0; s<g[i].size(); s++) {      
        for (int j=0; j<g[i][s]; j++) {
          y2[aux1]=mu_mat(i,s);
          aux1++;
        }      
      }
    }
    aux1=0;
    
    
    // Update regression part
    y3=y-y1-y2;
    
    XX_regr = O.t()*O;
    Xy_regr = O.t()*y3;
    Sigma_regr=(arma::inv(Sigma0_regr)+(XX_regr/psi[nn-1]));
    Sigmainv_regr=arma::inv(Sigma_regr);
    mu_regr=Sigmainv_regr*(Sigma0_regr*beta0_regr+Xy_regr/psi[nn-1]);
    beta.row(nn) =  arma::mvnrnd(mu_regr,Sigmainv_regr).t();
    y3 = O*beta.row(nn).t();	
	
    // Update psi
    res_regr  = y - y1 - y2 - y3; //
    ssr_regr=0;
    for (int j=0; j<N; j++) {
      ssr_regr+=pow(res_regr(j),2);
    }
    
    nunpsi= nu0psi+N;
    sigmanpsi= nu0psi*sigma0psi+ssr_regr;
    psi[nn]=1/as_scalar(arma::randg(1, arma::distr_param(nunpsi/2,2/sigmanpsi)));
    
    if ((nn>burn)&(nn%thin==0)){
		th << vectorise(THETA2).t() <<endl;
		bet << beta.row(nn) << endl;					  
		for (int j = 0; j < Lambda.n_cols; j++)
			lambda_norm << arma::norm(Lambda.col(j))<<" ";
		lambda_norm << endl;
		ps << psi[nn] << endl;
		kkk << k << endl;
		mu_print << vectorise(mu_mat,1);
									   
										  
		Lambda_tot << vectorise(Lambda).t() << endl;
		eta << vectorise(eta2).t() <<endl;
		m_print << m[nn] <<endl;									
						  
		h_is << vectorise(h, 1) <<endl;
		z_is << vectorise(zeta, 1) <<endl;
		alpha_print << alpha.row(nn) <<endl;
		beta_print << GARCHbeta(nn) <<endl;
		residuals << res_regr.t() << endl;			
	}
	cout << "iteration" << nn << endl; 
  
  } 
   // END GIBBS
  
  // Outout print
  
  th.close();
  bet.close();
  ps.close();
  lambda_norm.close();
  kkk.close();
  mu_print.close();
  Lambda_tot.close();
  eta.close();
  m_print.close();
  h_is.close();
  z_is.close();
  alpha_print.close();
  beta_print.close();
  acc_prob_alpha.close();
  acc_prob_beta.close();
  alpha_proposal.close();
  beta_proposal.close();
residuals.close();  		  
  return List::create(); //
}

