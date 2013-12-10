####################################################
# To recreate the logo of HSF from the MEME matrix #
####################################################

#install if necessary
#source("http://bioconductor.org/biocLite.R")
#biocLite("seqLogo")
library(seqLogo)

aa <- c(0.62,  0.  ,  0.96,  0.86,  0.46,  0.62,  0.  ,  0.96,  0.86,  0.46,  0.62,  0.  ,  0.96,  0.86,  0.46)
cc <- c(0.28,  0.  ,  0.04,  0.02,  0.02,  0.28,  0.  ,  0.04,  0.02,  0.02,  0.28,  0.  ,  0.04,  0.02,  0.02)
gg <- c(0.08,  1.  ,  0.  ,  0.1 ,  0.28,  0.08,  1.  ,  0.  ,  0.1,   0.28,  0.08,  1.  ,  0.  ,  0.1 ,  0.28)
tt <- c(0.02,  0.  ,  0.  ,  0.02,  0.24,  0.02,  0.  ,  0.  ,  0.02,  0.24,  0.02,  0.  ,  0.  ,  0.02,  0.24)

ttt <- rev(aa)
ggg <- rev(cc)
ccc <- rev(gg)
aaa <- rev(tt)

df <- data.frame(aaa,ccc,ggg,ttt)

# Divide the frequency by the row sum i.e. proportions
proportion <- function(x){
   rs <- sum(x);
   return(x / rs);
}

pwm <- apply(df, 1, proportion)
pwm <- makePWM(pwm)

png("logo.png", width=500, height=250)
seqLogo(pwm)
dev.off()
