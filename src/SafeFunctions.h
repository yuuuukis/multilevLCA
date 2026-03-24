#ifndef SAFEFUNCTIONS_H
#define SAFEFUNCTIONS_H

double probcheck(double prob);
double LogOfSum(arma::vec vLog);
arma::mat GammaCheck(arma::mat mGamma, int iJ);
arma::vec OmegaCheck(arma::vec vOmega, int iK);
arma::vec PostCheck(arma::vec vPost, int iK);
double LogDensityCheck(double dLLK);
double DensityCheck(double dDensity);
double MixtDensityScale(arma::vec vOmega, arma::vec vD_log, int iJ);
// double zerocheck(double value);

#endif
