#include <RcppEigen.h>
#include<cmath>
//[[Rcpp::depends(RcppEigen)]]

//[[Rcpp::export]]
double get_rf_fitness_cpp(Eigen::VectorXi& p, Eigen::MatrixXd& peaks, double wl){
  if(p.minCoeff()<1) return -1000.;
  int nrow=peaks.rows();
  int ncol=peaks.cols();
  double scalef = 1.0/(3.14159265*sqrt((double)nrow));
  double fitness = 0;
  for(int i = 0; i<nrow;i++){
    double d = 0;
    for(int j = 0; j<ncol; j++){
      d+=pow((double)p(j)-peaks(i,j),2);
    }
    d=sqrt(d);
    fitness+=sin(d/wl)*scalef;
  }
  
  return fitness;
}

//[[Rcpp::export]]
Eigen::VectorXi peak_step(Eigen::VectorXi& p, Eigen::MatrixXd& peaks, double wl){

  Eigen::VectorXi p_tmp = p;
  Eigen::VectorXi p_best = p;
  double f_best=get_rf_fitness_cpp(p,peaks,wl); 
  
  for(int i = 0; i<p_tmp.size(); i++){
    p_tmp = p;
    p_tmp(i)+=1;
    double f_tmp = get_rf_fitness_cpp(p_tmp,peaks,wl); 
    if(f_tmp>f_best) p_best = p_tmp;
    f_best = std::max(f_tmp,f_best);
    p_tmp = p;
    p_tmp(i)-=1;
    f_tmp = get_rf_fitness_cpp(p_tmp,peaks,wl); 
    if(f_tmp>f_best) p_best = p_tmp;
    f_best = std::max(f_tmp,f_best);
  }
  
  return p_best;
}
//[[Rcpp::export]]
Eigen::VectorXi find_peak(Eigen::VectorXi p, Eigen::MatrixXd peaks, double wl){
  Eigen::VectorXi ptmp = peak_step(p,peaks,wl);
  while(ptmp!=p){
    p=ptmp;
    ptmp = peak_step(p,peaks,wl);
  }
  return p;
}



