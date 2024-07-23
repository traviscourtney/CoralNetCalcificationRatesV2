######################################################################
##’ @title Area-normalized scaling of Indo-Pacific ReefBudget calcification and bioerosion rates for use with CoralNet
##’ ReefBudget2CoralNet
##’ @author Travis A Courtney, PhD
##’ @contact travis.courtney@upr.edu
##’ @date 2024-05-22
##’ @log Version 2.0
######################################################################

library(triangle)

set.seed(999)

coeff_mean=1.24
coeff_lwr=0.66
coeff_upr=1.98
int_mean=0
int_lwr=0
int_upr=0
length_mean=13
length_lower=5
length_upper=89.5
rugosity_mean=3.3
rugosity_sd=0.8

ReefBudget2CoralNet = function(coeff_mean,coeff_lwr,coeff_upr,int_mean,int_lwr,int_upr,length_mean,length_lower,length_upper,rugosity_mean,rugosity_sd)
{
  nsim=100000
  s=rtriangle(n=nsim,c=length_mean,a=length_lower,b=length_upper)
  c=rtriangle(n=nsim,c=coeff_mean,a=coeff_lwr,b=coeff_upr)
  i=rtriangle(n=nsim,c=int_mean,a=int_lwr,b=int_upr)
  r=rnorm(n=nsim,mean=rugosity_mean,sd=rugosity_sd)
  n=100/s
  g=(n*(c*s*r+i))/10
  calc=round(quantile(g,probs=0.5,na.rm=TRUE),digits=2)
  calc_lower=round(quantile(g,probs=0.25,na.rm=TRUE),digits=2)
  calc_upper=round(quantile(g,probs=0.75,na.rm=TRUE),digits=2)
  calc_rates=data.frame(calc,calc_lower,calc_upper)
  return(calc_rates)
}

ReefBudget2CoralNet(coeff_mean,coeff_lwr,coeff_upr,int_mean,int_lwr,int_upr,length_mean,length_lower,length_upper,rugosity_mean,rugosity_sd)