
load("../../../gmma/release/nonreleased_output/gmma_result.rda")

# Multi-mutant data frame
use_mask = mutant$gmma == "use"
mi_use = which(use_mask)
max_n_mut = max(mutant[mi_use,"N_sub"])
d = data.frame(table(mutant[mi_use,"N_sub"]))
colnames(d) = c("n_mut","obs")
rownames(d) = d$n_mut
d$n_mut = as.numeric(levels(d$n_mut)[d$n_mut])  # numbers in stead of factors here
d$frac_lib = d$obs/sum(d$obs)
t = table(mutant[which(mutant$signal > .5 & use_mask), 'N_sub'])   #assume mutant$active is either 1 or 0
d[names(t),'active'] = t
d$frac_act = d$active/d$obs
d$signal_avg = sapply(d$n_mut, function(n) { mean(mutant[which(mutant$N_sub == n & use_mask),"signal"]) })
d$signal_med = sapply(d$n_mut, function(n) { median(mutant[which(mutant$N_sub == n & use_mask),"signal"]) })


# quartz(width=6, height=7.5)
# tiff("library_composition_300dpi.tiff", width=6, heigh=4, units="cm", res=300, pointsize=7, compression="lzw")
jpeg(file="library_composition_300dpi.jpg", width=8, height=6, units="cm", res=300, pointsize=8, quality=90)

par(mar=c(4,4,1,2)+.1)
plot(d$n_mut, d$frac_lib, ylim=c(0,1), col="white", xlab="", ylab="")
title(xlab="Substitutions per variant", ylab="Fraction", line=2.4)
points(d$n_mut, d$frac_lib, pch=16, lwd=1.5, type="o", col=2)
points(d$n_mut, d$frac_act, pch=15, lwd=1.5, type="o", col=4)
legend("topright",c("Library composition","Fraction complementing"), pch=c(16,15), lty=1, lwd=1.5, col=c(2,4), cex=0.9)
dev.off()



require(Matrix)
require(igraph)
load("../../../gmma/release/nonreleased_output/gmma_graph.rda")

breaks = seq(0,1200,10)
h  = hist(degree(graph_clean)[which(V(graph_clean)$type)], breaks=breaks, plot=F)

# Plot substitution node degree distribution
# tiff("library_complexity_300dpi.tiff", width=6, heigh=4, units="cm", res=300, pointsize=7, compression="lzw")
jpeg(file="library_complexity_300dpi.jpg", width=8, height=6, units="cm", res=300, pointsize=8, quality=90)

par(mar=c(4,4,1,2)+.1)
hi = which(h$counts > 0)
plot(h$mids[hi], h$counts[hi], col="white", xlab="", ylab="", xlim=c(1,1200), log="y", yaxp=c(1,100,1))
title(xlab="Variants per substitution", ylab="Counts [substitutions]", line=2.4)
# hpi = which(hp$counts > 0)
# lines(hp$mids[hpi], hp$counts[hpi], type="o", pch=18, col=4)
lines(h$mids[hi], h$counts[hi], type="o", pch=20, col=2)
# legend("topright", c("Substitution nodes","Projected substitution graph"), pch=c(20,18), lty=1, col=c(2,4))

dev.off()
