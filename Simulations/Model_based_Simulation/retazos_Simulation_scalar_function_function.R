X          = rbinom(sub,1,0.5)

epsilon    = matrix(rnorm(sub*len,0,sd_epsilon), ncol=len,nrow=(sub))

design_M = cbind(int=1,X=X)
coeff_M  = cbind(delta1,alpha)
colnames(coeff_M) = c("delta1","alpha")

M    = t(tcrossprod(coeff_M, design_M))+epsilon

eta  = matrix(rnorm(sub*len,0,sd_eta), ncol=len,nrow=(sub))

Y    = t(sapply(1:sub, function(i) create_functional_y(x=X[i], m=M[i,], etai=eta[i,], delta2=delta2, gamma=gamma, beta=beta)))

dta  = data.frame(X=X,M=M,Y=Y)

N= sub; y = Y; x=X; m=M; penalty_ff=c(3,3); len=

fit  <- pffr(y ~ x + ff(m,limits="s<t", integration="riemann", basistype = "s", splinepars=list(k=ifelse(N < 52, N-2, 50),m=c(3,2)))) # Defaults to quartic (m[1]=3) P-splines (bs="ps") with 2nd derivative order penalty (m[2]=2), and at most 50-dimensional basis 

system.time(fit2 <- pffr(y ~ x + ff(m,limits="s<t", integration="riemann", splinepars=list(bs="pss", k=ifelse(N < 52, N-2, 50),m=c(3,2)))) # Defaults to quartic (m[1]=3) P-splines (bs="ps") with 2nd derivative order penalty (m[2]=2), and at most 50-dimensional basis 
)
system.time(
fit3 <- pffr(y ~ x + ff(m,limits="s<t", integration="riemann", splinepars=list(bs="pss",m=list(c(2, 1), c(2,1))))) # Defaults to quartic (m[1]=3) P-splines (bs="ps") with 2nd derivative order penalty (m[2]=2), and at most 50-dimensional basis 
)

system.time(
  fit3 <- pffr(y ~ x + ff(m,limits="s<t", integration="riemann", splinepars=list(bs="pss"))) # Defaults to quartic (m[1]=3) P-splines (bs="ps") with 2nd derivative order penalty (m[2]=2), and at most 50-dimensional basis 
)

plot(fit)
plot(fit2)
plot(fit3)