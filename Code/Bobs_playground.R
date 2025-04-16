boxplot(outs$ETR.21p[outs$setTleaf==25], outs$ETR.0p[outs$setTleaf==25])


plot(c(median(outs$ETR.21p[outs$setTleaf==25]), median(outs$ETR.0p[outs$setTleaf==25])), type="l",
     ylim=c(40,90), col='blue', lwd=4)
lines(c(median(outs$ETR.21p[outs$setTleaf==30]), median(outs$ETR.0p[outs$setTleaf==30])), 
      col='orange', lwd=4)
lines(c(median(outs$ETR.21p[outs$setTleaf==35]), median(outs$ETR.0p[outs$setTleaf==35])), 
      col='red', lwd=4)




par(mar=c(5,5,1,1))
plot(outs$JO2, outs$pr.real, pch=21, bg=c('blue','orange','red')[as.factor(outs$setTleaf)])

plot(outs$JO2.percent, outs$pr.percent, 
     pch=21, bg=c('blue','orange','red')[as.factor(outs$setTleaf)],
     cex.lab=2, cex=2)

cor.test(outs$JO2.percent, outs$pr.percent)


plot(outs$JO2.percent, outs$pr.percent)
abline(0,1)



plot(c(median(outs$ETR.21p[outs$setTleaf==25]),
       median(outs$ETR.0p[outs$setTleaf==25])),
     type="l", ylim=c(45, 85), col="blue", lwd=4)
lines(c(median(outs$ETR.21p[outs$setTleaf==30]),
       median(outs$ETR.0p[outs$setTleaf==30])),
      col="orange", lwd=4)
lines(c(median(outs$ETR.21p[outs$setTleaf==35]),
        median(outs$ETR.0p[outs$setTleaf==35])),
      col="red", lwd=4)



par(mfrow=c(3,3))

for(sp in seq_along(unique(outs$sp))){
  focsp <- unique(outs$sp)[sp]
  plot(c(median(outs$ETR.21p[outs$setTleaf==25 & outs$sp==focsp]), 
         median(outs$ETR.0p[outs$setTleaf==25 & outs$sp==focsp])), 
       type="l", ylim=c(0, 150), col=NA, lwd=4, main=focsp)
  lines(c(median(outs$ETR.21p[outs$setTleaf==25 & outs$sp==focsp]), 
          median(outs$ETR.0p[outs$setTleaf==25 & outs$sp==focsp])), 
        type="l", ylim=range(c(outs$ETR.21p, outs$ETR.0p)), col='blue', lwd=2)
  lines(c(median(outs$ETR.21p[outs$setTleaf==30 & outs$sp==focsp]), 
          median(outs$ETR.0p[outs$setTleaf==30 & outs$sp==focsp])), 
        type="l", ylim=range(c(outs$ETR.21p, outs$ETR.0p)), col='orange', lwd=2)
  lines(c(median(outs$ETR.21p[outs$setTleaf==35 & outs$sp==focsp]), 
          median(outs$ETR.0p[outs$setTleaf==35 & outs$sp==focsp])), 
        type="l", ylim=range(c(outs$ETR.21p, outs$ETR.0p)), col='red', lwd=2)
}




par(mfrow=c(3,3))

for(sp in seq_along(unique(outs$sp))){
  focsp <- unique(outs$sp)[sp]
  
  plot(median(outs$pr.real[outs$setTleaf==25 & outs$sp==focsp]),
       median(outs$ETR.21p[outs$setTleaf==25 & outs$sp==focsp]) - 
           median(outs$ETR.0p[outs$setTleaf==25 & outs$sp==focsp]), 
       col="blue", main=focsp, ylim=c(-10,50), xlim=c(0,13), pch=16, cex=2)
  points(median(outs$pr.real[outs$setTleaf==30 & outs$sp==focsp]),
       median(outs$ETR.21p[outs$setTleaf==30 & outs$sp==focsp]) - 
         median(outs$ETR.0p[outs$setTleaf==30 & outs$sp==focsp]), 
       col="orange", main=focsp, pch=16, cex=2)
  points(median(outs$pr.real[outs$setTleaf==35 & outs$sp==focsp]),
         median(outs$ETR.21p[outs$setTleaf==35 & outs$sp==focsp]) - 
           median(outs$ETR.0p[outs$setTleaf==35 & outs$sp==focsp]), 
         col="red", main=focsp, pch=16, cex=2)
}



