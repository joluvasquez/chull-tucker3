p <- parse('R/CHull_Ceulemans-Kiers.R')
for (ex in p) if (is.call(ex) && identical(ex[[1]],as.name('<-')) && is.call(ex[[3]]) && identical(ex[[3]][[1]],as.name('function'))) eval(ex)
set.seed(432)
s <- generar_tucker3(7,6,5,3,2,2,0)
manual <- array(0,c(7,6,5))
for(i in 1:7) for(j in 1:6) for(k in 1:5) for(a in 1:3) for(b in 1:2) for(c in 1:2) manual[i,j,k] <- manual[i,j,k]+s$A[i,a]*s$B[j,b]*s$C[k,c]*s$G[a,b,c]
cat('Independent reconstruction max error:',max(abs(manual-s$M)),'\n')
tab <- ajustes_tucker3(s$X,5)
cat('Noiseless true fit:',subset(tab,P==3 & Q==2 & R==2)$fit,'\n')
# Verify the hull against base R geometry on random strictly increasing fits.
for(z in 1:100) {
 x <- cumsum(runif(20,.1,2)); y <- cumsum(runif(20)); y <- y/max(y)
 t <- data.frame(P=1:20,Q=1,R=1,fp=x,fit=y)
 actual <- CHull(t)$hull$P
 # Points lying on or above every chord are not a simple test; use chull cycle.
 v <- chull(x,y); start <- which(v==1); cycle <- c(v[start:length(v)],if(start>1) v[1:(start-1)]); a <- cycle[1:which(cycle==20)]
 b <- c(1,rev(cycle[which(cycle==20):length(cycle)]))
 height <- function(idx) mean(y[idx]-(y[1]+(y[20]-y[1])*(x[idx]-x[1])/(x[20]-x[1])))
 expected <- if(height(a)>=height(b)) a else b
 stopifnot(identical(as.integer(actual),sort(as.integer(expected))))
}
cat('100 independent hull geometry tests: OK\n')


stopifnot(max(abs(manual-s$M))<1e-12)
stopifnot(abs(subset(tab,P==3 & Q==2 & R==2)$fit-1)<1e-12)
for(eps in c(0,1e-15)) {
  plateau <- data.frame(P=1:4,Q=1,R=1,fp=1:4,S=1:4,fit=c(.1,.6,1,1-eps))
  for(metric in c('fp','S')) stopifnot(CHull(plateau,metric)$seleccionado$P==3)
}
set.seed(321)
sim <- generar_tucker3(27,27,27,3,2,2,sqrt(.6/.4))
stopifnot(nrow(ajustes_tucker3(sim$X))==341)
stopifnot(abs(sim$error_nominal-sim$error_energia_separada)<1e-10)
stopifnot(abs(sim$SS_X-sim$SS_senal-sim$SS_ruido-sim$termino_cruzado)<1e-10)
cat('Todas las pruebas satisfactorias.\n')
