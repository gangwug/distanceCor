## permutation-test of distance correlation coefficient
rm(list=ls())
#install.packages("energy")
library(energy)
run_start=proc.time()
dcorPermu <- function(v1,v2,permutation.num=100)
{
    t = dcor(v1,v2);                                   #t = true correlation value 
    f = c();                                           #f = recipient vector
    set.seed(1000);                                    #then for one expression profile, with different times analysis, the output results will be constant
    for(i in 1:permutation.num)                        # Script to run 10,000 and store in a vector
    {
      temp= sample(v2, length(v1), replace=TRUE); 
      f[i] = dcor(v1, temp);
    }
    prob = ecdf(f);                                    # Ranking the true value in the sea of false values
    p_val = prob(t);                                   # Non-parametric p-value estimate
    p_val = 1-p_val;
    return (c(t,p_val));
}

run_start=proc.time();
## ten levels (eg., dosing levels)
dose.array <- c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 1000);
ary.expM <- read.delim("exampleData.txt", stringsAsFactors = FALSE);
dary.Labname <- c("RdCor","RdCorPvalue","RmeanExp","RmaxExp",
                "CdCor","CdCorPvalue","CmeanExp","CmaxExp");
genesym <- as.character(ary.expM$ID)
ary.rhoM <- as.matrix(ary.expM[,-1])
dary.rhoP <- NULL
for (i in 1:nrow(ary.rhoM))             
{
  tep.rhoP <- NULL
  for (j in 1:2)
  {
    st <- (j-1)*11+1;
    ed <- (j-1)*11+11;
    tep.exp <- ary.rhoM[i,st:ed];
    tep.sd <- sd(tep.exp);
    if (tep.sd == 0)
    {  tep.rhoP <- c(tep.rhoP,NaN,1,mean(tep.exp),max(tep.exp));  }                                 ##if constantly express, set P=1
    if (tep.sd > 0)
    {
      tepa <- dcorPermu(v1=dose.array,v2=tep.exp,permutation.num=100);                             #since few data in this file comparing with RNA-Seq, more permutation times are used
      tep.rhoP <- c(tep.rhoP,tepa,mean(tep.exp),max(tep.exp));
    } 
  }
  dary.rhoP <- rbind(dary.rhoP,tep.rhoP);
}
dimnames(dary.rhoP) <- list("r"=genesym,"c"=dary.Labname);
dary.rhoP <- as.data.frame(dary.rhoP)
dary.rhoP$ID <- row.names(dary.rhoP)
write.table(dary.rhoP[,c("ID", dary.Labname)], file="example.dCorPermuTest.txt", quote=FALSE, sep="\t", row.names = FALSE);
print ( c("Time used:", (proc.time()-run_start)) )

