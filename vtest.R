library(metafor)


dat <- read.csv("rmatest.csv")
dat <- read.csv("cohentest.csv")

### construct dataset and var-cov matrix of the correlations
#dat <- dat %>% filter(Env=="Temp")

dat <- dat[which(!is.na(dat$Z_cor)),]
dat <- dat %>% filter(Z_cor>-1)
dat <- dat %>% filter(Z_cor<1)

tmp <- rcalc(Z_cor ~ var_1 + var_2 | effect_ID, ni=Length, data=dat)
V <- tmp$V
V<-as.matrix(V)
#dat <- tmp$dat


### multivariate random-effects model

theRegressionModel <- ~Lat + Long

theResults <- rma.mv(yi = Z_cor,
                  V = V,
                  mods = ~ Lat + Long + cor_ID,
                  random = ~ 1 | effect_ID,
                  #struct = "UN",
                  data = dat,
                  control = list(optimizer = "nlminb", maxiter=1000))
