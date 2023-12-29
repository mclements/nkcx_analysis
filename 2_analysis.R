## Implement an extension of Day and Walter (1984)
## Mark Clements

library(minqa) # bobyqa
library(parallel) # mclapply
library(Rcpp) # sourceCpp
library(xtable) # xtable
logit = binomial()$linkfun
expit = binomial()$linkinv
pow=function(a,b) a^b
      
## Gauss-Kronrod -- log-logistic for f1 and otherwise exponentials
sourceCpp(code="
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <vector>
#include <cmath>
using namespace boost::math::quadrature;
class ScreeningModel {
public:
size_t fulln, n;
std::vector<double> ti, tj;
double a1, b1, rate2, rate3, pcure1, beta2, beta3, betaTx, betaUnder, error, error2, error3;
ScreeningModel(size_t fulln, std::vector<double> ti, double a1, double b1, double rate2, double rate3, double pcure1, double beta2, double beta3, double betaTx, double betaUnder) : fulln(fulln), ti(ti), a1(a1), b1(b1), rate2(rate2), rate3(rate3), pcure1(pcure1), beta2(beta2), beta3(beta3), betaTx(betaTx), betaUnder(betaUnder) { }
// double f1(double s) { return rate1*std::exp(-rate1*s); }
// double S1(double s) { return std::exp(-rate1*s); }
double S1(double x) { return pcure1 + (1-pcure1)/(1+pow(x/b1,a1)); }
double f1(double x) { double a2=a1-1; return (1-pcure1)*(a1/b1)*pow(x/b1,a2)/pow(1+pow(x/b1,a1),2); }
double f2(double s) { return rate2*std::exp(-rate2*s); }
double S2(double s) { return std::exp(-rate2*s); }
double f3(double s) { return rate3*std::exp(-rate3*s); }
double S3(double s) { return std::exp(-rate3*s); }
// *NB*: use tj and n in the calculations!
void setup(double t) {
  tj.resize(0);
  tj.push_back(0.0);
  for (size_t i=0; i<fulln; i++) {
    if (ti[i]<t) tj.push_back(ti[i]);
    else break;
  }
  tj.push_back(t);
  n = tj.size()-2;
}
double W(double s, bool reset = true) {
  if (reset) setup(s);
  return S1(s);
}
double X(double t, bool reset = true) {
  if (reset) setup(t);
  double value = 0.0;
  for (size_t i=0; i<=n; i++) {
    auto fn = [&](double x) { return f1(x)*S2(t-x)*pow(beta2,n-i); };
    value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
  }
  return value;
}
double Y(double t, bool reset=true) {
  if (reset) setup(t);
  double value = 0.0, error, error2;
  for (size_t i=0; i<=n; i++) {
    for (size_t j=i; j<=n; j++) {
      auto fn = [&](double s) { 
        auto inner = [&](double u) { return f1(s)*f2(u-s)*S3(t-u)*pow(beta2,j-i)*pow(beta3,n-j)*pow(1-betaUnder,n-j); };
        return gauss_kronrod<double, 15>::integrate(inner, std::max(s,tj[j]), tj[j+1], 5, 1e-6, &error2);
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
    }
  }
  return value;
}
double I(double t, bool reset=true) {
  if (reset) setup(t);
  double value = 0.0;
  for (size_t i=0; i<=n; i++) {
    for (size_t j=i; j<=n; j++) {
      auto fn = [&](double s) { 
        auto inner = [&](double u) { return f1(s)*f2(u-s)*f3(t-u)*pow(beta2,j-i)*pow(beta3,n-j)*pow(1-betaUnder,n-j); };
        return gauss_kronrod<double, 15>::integrate(inner, std::max(s,tj[j]), tj[j+1], 5, 1e-6, &error2);
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
    }
  }
  return value;
}
double Z(double t, bool reset=true) {
  if (reset) setup(t);
  double value = 0.0;
  for (size_t i=0; i<=n; i++) {
    for (size_t j=i; j<=n; j++) {
    for (size_t k=j; k<=n; k++) {
      auto fn = [&](double s) { 
        auto fn2 = [&](double u) { 
          auto fn3 = [&](double v) {return f1(s)*f2(u-s)*f3(v-u)*pow(beta2,j-i)*pow(beta3,k-j); };
          return gauss_kronrod<double, 15>::integrate(fn3, std::max(u,tj[k]), tj[k+1], 5, 1e-6, &error3); 
        };
        return gauss_kronrod<double, 15>::integrate(fn2, std::max(s,tj[j]), tj[j+1], 5, 1e-6, &error2); 
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
    }
    }
  }
  if (n>1) {
  for (size_t i=0; i<n; i++) {
    for (size_t j=i; j<n; j++) {
    for (size_t k=j+1; k<=n; k++) {
      auto fn = [&](double s) { 
        auto fn2 = [&](double u) { return f1(s)*f2(u-s)*S3(tj[k]-u)*pow(beta2,j-i)*pow(beta3,k-j-1)*(1-beta3); };
        return gauss_kronrod<double, 15>::integrate(fn2, std::max(s,tj[j]), tj[j+1], 5, 1e-6, &error2); 
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
    }
    }
  }
  }
  return value;
}
double Wstar(double t, bool reset=true) {
  if (reset) setup(t);
  if (n==0) return 0.0;
  double value = 0.0;
  for (size_t i=0; i<n; i++) {
    for (size_t j=i+1; j<=n; j++) {
      auto fn = [&](double s) {  return f1(s)*S2(tj[j]-s)*pow(beta2,j-i-1)*(1-beta2)*(1-betaTx)*S1(t-tj[j]);
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
    }
  }
  return value;
}
double Wstarj(double t, size_t j, bool reset=true) {
  if (reset) setup(t);
  if (n==0 || j<1 || j>n) return 0.0;
  double value = 0.0;
  for (size_t i=0; i<j; i++) {
      auto fn = [&](double s) {  return f1(s)*S2(tj[j]-s)*pow(beta2,j-i-1)*(1-beta2)*(1-betaTx)*S1(t-tj[j]);
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
  }
  return value;
}
double Xstar(double t, bool reset=true) {
  if (reset) setup(t);
  if (n==0) return 0.0;
  double value = 0.0;
  for (size_t i=0; i<n; i++) {
    for (size_t j=i+1; j<=n; j++) {
      auto fn = [&](double s) {  
        double expr=0.0;
        for (size_t k=j; k<=n; k++) {
          auto fn2 = [&](double u) { return f1(u-tj[j])*S2(t-u)*pow(beta2,n-k); };
          expr += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
        }
        return f1(s)*S2(tj[j]-s)*pow(beta2,j-i-1)*(1-beta2)*(betaTx*pow(beta2,n-j)*S2(t-tj[j]) +
          (1-betaTx)*expr);
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
    }
  }
  return value;
}
double Xstarj(double t, size_t j, bool reset=true) {
  if (reset) setup(t);
  if (n<1 || j<1 || j>n) return 0.0;
  double value = 0.0;
  for (size_t i=0; i<j; i++) {
      auto fn = [&](double s) {  
        double expr=0.0;
        for (size_t k=j; k<=n; k++) {
          auto fn2 = [&](double u) { return f1(u-tj[j])*S2(t-u)*pow(beta2,n-k); };
          expr += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
        }
        return f1(s)*S2(tj[j]-s)*pow(beta2,j-i-1)*(1-beta2)*(betaTx*pow(beta2,n-j)*S2(t-tj[j]) +
          (1-betaTx)*expr);
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
  }
  return value;
}
double Ustarstar(double t, bool reset=true) {
  if (reset) setup(t);
  if (n<2) return 0.0;
  double value = 0.0;
  for (size_t i=0; i<=n-2; i++) {
    for (size_t j=i+1; j<=n-1; j++) {
      auto fn = [&](double s) {
        double expr1=0.0, expr2=0.0;
        for (size_t k=j+1; k<=n; k++) {
          expr1 += pow(beta2,k-j-1)*S2(tj[k]-tj[j]);
        }
        for (size_t k=j; k<n; k++) {
          for (size_t m=k+1; m<=n; m++) {
            auto fn2 = [&](double u) { return f1(u-tj[j])*S2(tj[m]-u)*pow(beta2,m-k-1); };
            expr2 += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
          }
        }
        return f1(s)*S2(tj[j]-s)*pow(beta2,j-i-1)*pow(1-beta2,2)*(betaTx*expr1 + (1-betaTx)*expr2);
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
    }
  }
  return value;
}
double Ystar(double t, bool reset=true) {
  if (reset) setup(t);
  if (n==0) return 0.0;
  double value = 0.0;
  for (size_t i=0; i<=n; i++) {
    for (size_t j=i; j<=n; j++) {
      auto fn = [&](double s) {
        double expr1=0.0, expr2=0.0;
        for (size_t k=j; k<=n; k++) {
          auto fn2 = [&](double u) { return f2(u-tj[j])*S3(t-u)*pow(beta2,k-j)*pow(beta3,n-k); };
          expr1 += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
        }
        for (size_t k=j; k<=n; k++) {
          for (size_t m=k; m<=n; m++) {
            auto fn2 = [&](double u) { 
              auto fn3 = [&](double v) { return f1(u-tj[j])*f2(v-u)*S3(t-v)*pow(beta2,m-k)*pow(beta3,n-m); };
              return gauss_kronrod<double, 15>::integrate(fn3, std::max(u,tj[m]), tj[m+1], 5, 1e-6, &error3);
            };
            expr2 += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
          }
        }
        return f1(s)*S2(tj[j]-s)*pow(beta2,j-i-1)*(1-beta2)*(betaTx*expr1 + (1-betaTx)*expr2);
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
    }
  }
  return value;
}
double Ystarj(double t, size_t j, bool reset=true) {
  if (reset) setup(t);
  if (n<1 || j<1 || j>n) return 0.0;
  double value = 0.0;
  for (size_t i=0; i<j; i++) {
      auto fn = [&](double s) {
        double expr1=0.0, expr2=0.0;
        for (size_t k=j; k<=n; k++) {
          auto fn2 = [&](double u) { return f2(u-tj[j])*S3(t-u)*pow(beta2,k-j)*pow(beta3,n-k); };
          expr1 += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
        }
        for (size_t k=j; k<=n; k++) {
          for (size_t m=k; m<=n; m++) {
            auto fn2 = [&](double u) { 
              auto fn3 = [&](double v) { return f1(u-tj[j])*f2(v-u)*S3(t-v)*pow(beta2,m-k)*pow(beta3,n-m); };
              return gauss_kronrod<double, 15>::integrate(fn3, std::max(u,tj[m]), tj[m+1], 5, 1e-6, &error3);
            };
            expr2 += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
          }
        }
        return f1(s)*S2(tj[j]-s)*pow(beta2,j-i-1)*(1-beta2)*(betaTx*expr1 + (1-betaTx)*expr2);
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
  }
  for (size_t i=0; i<j; i++) {
      auto fn = [&](double s) {
        double expr3=0.0;
        for (size_t k=i; k<j; k++) {
          auto fn2 = [&](double u) { return f1(s)*f2(u-s)*S3(t-u)*pow(beta2,k-i)*pow(beta3,n-k)*pow(1-betaUnder,j-k-1)*betaUnder; };
          expr3 += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
        }
        return expr3;
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
  }
  return value;
}
double Zstar(double t, bool reset=true) {
  return 1.0 - (W(t)+X(t)+Y(t)+Z(t)+Wstar(t)+Xstar(t)-Ystar(t)-Ustarstar(t));
}
double Istarj(double t, size_t j, bool reset=true) {
  if (reset) setup(t);
  if (n<1 || j<1 || j>n) return 0.0;
  double value = 0.0;
  for (size_t i=0; i<j; i++) {
      auto fn = [&](double s) {
        double expr1=0.0, expr2=0.0;
        for (size_t k=j; k<=n; k++) {
            auto fn2 = [&](double u) { return f2(u-tj[j])*f3(t-u)*pow(beta2,k-j)*pow(beta3,n-k); };
            expr1 += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
        }
        for (size_t k=j; k<=n; k++) {
            for (size_t m=k; m<=n; m++) {
              auto fn2 = [&](double u) { 
                auto fn3 = [&](double v) { return f1(u-tj[j])*f2(v-u)*f3(t-v)*pow(beta2,m-k)*pow(beta3,n-m); };
                return gauss_kronrod<double, 15>::integrate(fn3, std::max(u,tj[m]), tj[m+1], 5, 1e-6, &error3);
              };
              expr2 += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
            }
        }
        return f1(s)*S2(tj[j]-s)*pow(beta2,j-i-1)*(1-beta2)*(betaTx*expr1 + (1-betaTx)*expr2);
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
  }
  for (size_t i=0; i<j; i++) {
      auto fn = [&](double s) {
        double expr3=0.0;
        for (size_t k=i; k<j; k++) {
          auto fn2 = [&](double u) { return f1(s)*f2(u-s)*f3(t-u)*pow(beta2,k-i)*pow(beta3,n-k)*pow(1-betaUnder,j-k-1)*betaUnder; };
          expr3 += gauss_kronrod<double, 15>::integrate(fn2, tj[k], tj[k+1], 5, 1e-6, &error2);
        }
        return expr3;
      };
      value += gauss_kronrod<double, 15>::integrate(fn, tj[i], tj[i+1], 5, 1e-6, &error);
  }
  return value;
}
};
//[[Rcpp::export]]
std::vector<double> test(size_t n, double t, std::vector<double> tj, double a1=1.0, double b1=1.0, double rate2=0.1, double rate3=0.1, double pcure1=0.9, double beta2=0.2, double beta3=0.05, double betaTx=0.2, double betaUnder=0.05) {
ScreeningModel m(n, tj, a1, b1, rate2, rate3, pcure1, beta2, beta3, betaTx, betaUnder);
double W=m.W(t), X=m.X(t), Y=m.Y(t), Z=m.Z(t), Wstar=m.Wstar(t), Xstar=m.Xstar(t), 
       Ystar=m.Ystar(t), Ustarstar=m.Ustarstar(t); 
return {W,X,Y,Z,Wstar,Xstar,Ustarstar,Ystar,1.0-(W+X+Y+Z+Wstar+Xstar+Ystar+Ustarstar),m.I(t),
m.Istarj(t,1), m.Xstarj(t,1)*(1-beta2), m.Ystarj(t,1)*(1-beta3), m.Wstarj(t,1)+m.Xstarj(t,1)+m.Ystarj(t,1),m.Wstar(t)-m.Wstarj(t,1)-m.Wstarj(t,2),m.Xstar(t)-m.Xstarj(t,1)-m.Xstarj(t,2),m.Ystar(t)-m.Ystarj(t,1)-m.Ystarj(t,2)};
}
double like(size_t n, double t, std::vector<double> tj, int type, size_t j=0, double a1=1.0, double b1=30.0, double rate2=0.1, double rate3=0.1, double pcure1=0.9, double beta2=0.2, double beta3=0.05, double betaTx=0.2, double betaUnder = 0.05) {
ScreeningModel m(n, tj, a1, b1, rate2, rate3, pcure1, beta2, beta3, betaTx, betaUnder);
if (type==-1) return m.W(t);
if (type==-2) return m.X(t);
if (type==-3) return m.Y(t);
if (type==-4) return m.Z(t);
if (type==-5) return m.Wstarj(t,j);
if (type==-6) return m.Xstarj(t,j);
if (type==-7) return m.Ystarj(t,j);
if (type==-8) return m.Zstar(t,j);
if (type==1) return m.W(t)+m.X(t)+m.Y(t);
if (type==2) return m.Wstarj(t,j)+m.Xstarj(t,j)+m.Ystarj(t,j);
if (type==3) return m.Xstarj(t,j)*(1-beta2);
if (type==4) return m.I(t);
if (type==5) return m.Istarj(t,j);
if (type==6) return m.Y(t)*(1-beta3);
if (type==7) return m.Ystarj(t,j)*(1-beta3);
else return -1.0;
}
//[[Rcpp::export]]
std::vector<double> likes(Rcpp::List inputs, double a1=10.0, double b1=30.0, double rate2=0.1, double rate3=0.1, double pcure1=0.9, double beta2=0.2, double beta3=0.05, double betaTx=0.2, double betaUnder=0.05) {
using Rcpp::as;
std::vector<double> out(inputs.size());
std::vector<double> one = {1.0};
for (int i=0; i<inputs.size(); i++) {
Rcpp::List input = inputs(i);
std::vector<double> ti = as<size_t>(input(\"n\"))>0 ? as<std::vector<double>>(input(\"ti\")) : one;
out[i] = like(as<size_t>(input(\"n\")), as<double>(input(\"t\")), ti, as<int>(input(\"type\")), as<size_t>(input(\"j\")), a1, b1, rate2, rate3, pcure1, beta2, beta3, betaTx, betaUnder);
if (i % 10 == 0) R_CheckUserInterrupt();  /* be polite -- did the user hit ctrl-C? */
}
return out;
}")
print(est4 <- test(2,1,c(0.25,0.5)))

## Predictions
noScreening = function(t, type=-1) lapply(t, function(ti)
    list(t = ti, cohort = 1960, type = type, nhsil = 0, 
         ti = NULL, j = NA, tj = NA, j2 = NA, tj2 = NA, n = 0L))
stacked = function(x,ys, cols=1:6, ...) {
    plot(x,rowSums(ys), type="n", ...)
    lower=rep(0,length(x))
    for (i in 1:ncol(ys)) {
        upper = lower+ys[,i]
        polygon(c(x,rev(x)), c(lower,rev(upper)), col=cols[i], border=cols[i])
        lower = upper
    }
}
band = function(x,lower,upper,col="grey") {
    if (length(lower)==1) lower=rep(lower,length(upper))
    if (length(upper)==1) upper=rep(upper,length(lower))
    polygon(c(x,rev(x)),c(lower,rev(upper)),col=col,border=col)
}
pdf("~/Documents/clients/HEAP/nkcx_analysis/stacked.pdf")
par(mfrow=1:2)
ts = seq(15,65,len=301)
healthy = likes(noScreening(ts,type=-1))
HSIL = likes(noScreening(ts,type=-2))
Preclin = likes(noScreening(ts,type=-3))
Dx = likes(noScreening(ts,type=-4))
## Dx = 1-healthy-HSIL-Preclin
stacked(ts,cbind(HSIL,Preclin,Dx),xlab="Age (years)", ylab="Proportion", main="No screening", cols=1:4)
legend("topleft", legend=c("HSIL","Preclin","Dx"),
       fill=1:8)
##
singleScreening = function(t, type=-1, HSIL=FALSE) lapply(t, function(ti)
    list(t = ti, cohort = 1960, type = type, nhsil = 0, 
         ti = 30, j = if(HSIL) 1 else NA, tj = if (HSIL) 30 else NA, j2 = NA, tj2 = NA, n = 1L))
ts = seq(15,65,len=301)
healthy = likes(singleScreening(ts,type=-1))
HSIL = likes(singleScreening(ts,type=-2))
Preclin = likes(singleScreening(ts,type=-3))
Dx = likes(singleScreening(ts,type=-4))
healthyStar = likes(singleScreening(ts,type=-5,HSIL=TRUE))
HSILStar = likes(singleScreening(ts,type=-6,HSIL=TRUE))
PreclinStar = likes(singleScreening(ts,type=-7,HSIL=TRUE))
DxStar = likes(singleScreening(ts,type=-8,HSIL=TRUE))
stacked(ts, cbind(HSIL,Preclin,Dx,healthyStar,HSILStar,PreclinStar,DxStar),
        xlab="Age (years)", ylab="Proportion", cols=1:8, main="Single screen at 30")
legend("topleft", legend=c("HSIL","Preclin","Dx","healthyStar","HSILStar","PreclinStar","DxStar"),
       fill=1:8)
dev.off()

## Plots of HSIL onset based on fits
## log(a1) and log(b1) only
## neglli = function(theta) {
##     value = -sum(wt*log(unlist(mclapply(splitIndices(length(temp2),80), function(set) pmax(1e-100,likes(temp2[set], a1=exp(theta[1]), b1=exp(theta[2]), rate2=-log(1-0.05), rate3=-log(1-0.5))), mc.cores=80))))
##     cat(theta,value,"\n")
##     value
## }
## out = bobyqa(c(1.73642, 3.511851), neglli) # a1 and b1 *are* identifiable:)
outA = list(par=c(1.73642, 3.511851))
exp(outA$par)
## hess = optimHess(out$par, neglli)
hessA = structure(c(28492.3114159028, -4912.39072289318, -4912.39072289318, 
                   278980.309070903), dim = c(2L, 2L))
##
## negllii = function(theta) {
##     value = -wt*log(unlist(mclapply(splitIndices(length(temp2),80), function(set) pmax(1e-100,likes(temp2[set], a1=exp(theta[1]), b1=exp(theta[2]), rate2=-log(1-0.05), rate3=-log(1-0.5))), mc.cores=80)))
## }
calcScores = function(par, fn, eps=1e-5) {
    n = length(par)
    delta = function(i,eps) { out=rep(0,n); out[i]=eps; out}
    sapply(1:n, function(i) (fn(par+delta(i,eps))-fn(par-delta(i,eps)))/(2*eps))
}
## scores = calcScores(out$par, negllii) # matrix
source("screening-2-scores.R")
meat = Reduce("+",apply(scores,1,function(col) outer(col,col),simplify=FALSE))
bread = solve(hessA)
vcovA = bread %*% meat %*% bread
## apply(exp(out$par+outer(1.96*sqrt(diag(bread)),c(0,-1,1))),1,function(x) sprintf("%.2f (%.2f,%.2f)", x[1], x[2], x[3]))
cat(apply(exp(outA$par+outer(1.96*sqrt(diag(vcovA)),c(0,-1,1))),1,function(x) sprintf("%.2f (%.2f,%.2f)", x[1], x[2], x[3])), "\n")
##
pdf("/home/marcle/Documents/clients/HEAP/nkcx_analysis/onset_density.pdf")
par(mfrow=1:2)
grad = function(par, fn, eps=1e-5) {
    n = length(par)
    delta = function(i,eps) { out=rep(0,n); out[i]=eps; out}
    sapply(1:n, function(i) (fn(par+delta(i,eps))-fn(par-delta(i,eps)))/(2*eps))
}
f1 = function(x, a1=10, b1=30, pcure1=0.9) {a2=a1-1; (1-pcure1)*(a1/b1)*pow(x/b1,a2)/pow(1+pow(x/b1,a1),2) }
x=seq(0,60,length=1001)[-1]
g = grad(outA$par, function(par) log(f1(x,a=exp(par[1]),b=exp(par[2]),pcure=0.9)))
ci = exp(log(f1(x,a1=exp(outA$par[1]),b1=exp(outA$par[2]),pcure1=0.9)) +
    outer(sqrt(diag(g %*% vcov %*% t(g))),c(0,-1,1)))
matplot(x, ci, type="n", xlab="Age (years)", ylab="Density for onset of HSIL", main="Model A")
polygon(c(x,rev(x)),c(ci[,2],rev(ci[,3])),col="grey",border="grey")
lines(x, ci[,1])
##
outB = list(theta=c(6.0387605684121, 3.0977749365898, -1.57061328971769, 0.546152677743288, -8.22406189480543, 2.96611675405796, 6.67311456679623))
pars = with(outB, list(a1=exp(theta[1]), b1=exp(theta[2]), pcure1=expit(theta[3]), rate2=exp(theta[4]), rate3=exp(theta[5]), beta2=expit(theta[6]), beta3=expit(theta[7])))
plot(x, f1(x,a1=exp(outB$theta[1]),b1=exp(outB$theta[2]),pcure1=expit(outB$theta[3])),
     type="l", xlab="Age (years)", ylab="Density for onset of HSIL", main="Model B")
dev.off()

## Read in data (local only:()
if (local <- FALSE) {
    source("~/Downloads/temp-20231217.R")
    gc()
    screening1960 = screening[sapply(screening,"[[","cohort")==1960]
    ## Truncate follow-up after the second HSIL
    temp = lapply(screening1960, function(lst) {
        if (lst$nhsil>=2) {
            lst$type=3
            lst$t=lst$tj2
            lst$ti=lst$ti[1:(lst$j2-1)]
        }
        if (length(lst$ti)>0 && tail(lst$ti,1)==lst$t) {
            lst$ti = head(lst$ti,-1)
            ## lst$t = lst$t-1e-10
        }
        lst$n=length(lst$ti)
        lst
    })
    type = sapply(temp, "[[", "type")
    n1 = sum(type==1)
    set.seed(12345)
    index = c(which(type>1), which(type==1)[sample(n1, floor(n1/10))])
    temp2 = temp[index]
    wt = c(1,rep(sum(type>1)), rep(10, sum(type==1)))
    temp2=temp2[order(index)]
    wt=wt[order(index)]

    ## Table: description
    addTotals = function(tab) {
        withColTotals = cbind(tab,Total=apply(tab,1,sum))
        rbind(withColTotals,Total=apply(withColTotals,2,sum))
    }
    xtable(addTotals(table(sapply(screening1960, "[[", "type"), sapply(screening1960, "[[", "nhsil"))))


    ## log(a1) and log(b1) only
    neglli = function(theta) {
        value = -sum(wt*log(unlist(mclapply(splitIndices(length(temp2),5), function(set) pmax(1e-100,likes(temp2[set], a1=exp(theta[1]), b1=exp(theta[2]), rate2=-log(1-0.05), rate3=-log(1-0.5))), mc.cores=5))))
        cat(theta,value,"\n")
        value
    }
    out = bobyqa(c(log(10),log(30)), neglli)
}
