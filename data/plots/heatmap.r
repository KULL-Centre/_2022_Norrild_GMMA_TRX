
load("../../../gmma/release/gmma_result.rda")

nmut = length(mutant[,1])
nsubst = length(subst[,1])

# What2plot, gap not in heatmap
iss = which(substr(subst_sig$taa,1,1) != '*')
subst_sig_nogap = subst_sig[iss,]

# Convert to kJ/mol
cal2j = 8.314/1.985
subst_sig_nogap$ddG_glob = subst_sig_nogap$ddG_glob * cal2j
subst$ddG_glob = subst$ddG_glob * cal2j

# Order of amino acids
aa2int = seq(20)
# names(aa2int) = strsplit("ACDEFGHIKLMNPQRSTVWY","")[[1]]
names(aa2int) = strsplit("CDEKRHNQAGSTVMLIFYWP","")[[1]]

grad_range_min = -5
grad_range_max = 15
ns = abs(grad_range_min)*100
nd = abs(grad_range_max)*100
col_not_in_lib = "white"
col_native     = "black"
col_uncertain  = "gray80"
# -----------------------------------------------------------------------------------------------
# Variant :   destabilizing           neutral               stabilizing
      col_destab = "#ff5a6b";  col_neutral = "#ffff00";    col_stab = "#079700"
      # col_destab = "#47a600";  col_neutral = "#ffff00";    col_stab = "#e2002e"
      # col_destab = "#1b9d00";  col_neutral = "#f4e700";    col_stab = "#d4000d"
      # col_destab = "cyan";  col_neutral = "yellow";    col_stab = "red"
# Position: stable/conserved      neutral/tolerant     unstable/engineering potential
# If bith gradients pass through brown-black the difference between slightly stab and destab disapears
# white as tolerant works well because the small effects are easier to distinguis and I get a saturation gradient fom
#     tolerant to non- tolerant if both stab and destab are saturated colors
# -----------------------------------------------------------------------------------------------

col_grad = c(colorRampPalette( c(col_stab, col_neutral), space="rgb")(ns), colorRampPalette( c(col_neutral,col_destab), space="rgb")(nd))
col_breaks = seq(grad_range_min, grad_range_max, length.out=ns+nd+1)

first_resn = 49  # 48 does not have subst
last_resn = 95   # 96 and 97 does not have subst
wt_seq = as.character(residue[first_resn:last_resn,"wt"])
nres = length(wt_seq)

# quartz(width=5, height=10)
tiff("heatmap_300dpi.tiff", width=7.3, heigh=12, units="cm", res=300, pointsize=9, compression="lzw")
layout(matrix(c(1,1,3,2), ncol=2, ), width=c(4,1.1), height=c(3,2))
par(mar=c(2,4,2,1)+.1)

m = matrix(1, ncol=nres, nrow=20)
image(m, xaxt="n", yaxt="n", col=col_not_in_lib, ylim=c(1+.7/nres,-.7/nres), xlim=c(-.03,1.03))
axis(1, seq(0, 20, length=20)/20, names(aa2int), cex.axis=.8, las=1, gap.axis=0)
axis(3, seq(0, 20, length=20)/20, names(aa2int), cex.axis=.8, las=1, gap.axis=0)
res_lab = paste(wt_seq, seq(first_resn,first_resn-1+length(wt_seq)), sep="")

axis(2, seq(0, nres, length=nres)/nres, res_lab, cex.axis=.7, las=2, gap.axis=0)
# mask = rep_len(c(TRUE,FALSE), length.out=nres)
# axis(2,         (seq(0, nres, length=nres)/nres)[mask], labels=F, tcl=-.7, lwd=.6)
# axis(2,         (seq(0, nres, length=nres)/nres)[!mask], labels=F, tcl=-2.2, lwd=.6)
# axis(2,         (seq(0, nres, length=nres)/nres)[mask], res_lab[mask], cex.axis=.8, las=2, gap.axis=0, tick=F)
# axis(2, line=1.4, (seq(0, nres, length=nres)/nres)[!mask], res_lab[!mask], cex.axis=.8, las=2, gap.axis=0, tick=F)

# Mark all substitutions in library
m[] = NA
m[cbind(aa2int[substr(subst$taa,1,1)], subst$resi-first_resn+1)] = 1
image(m, col=col_uncertain, add=T)

# Mark native
m[] = NA
m[cbind(aa2int[wt_seq], seq(nres))] = 1
image(m, col=col_native, add=T)

# Mark stabilizing outside range
m[] = NA
ssi = which(subst_sig_nogap$ddG_glob < grad_range_min)
if (length(ssi) > 0) {
    m[cbind(aa2int[substr(subst_sig_nogap[ssi,'taa'],1,1)], subst_sig_nogap[ssi,'resi']-first_resn+1)] = 1
    image(m, col=col_stab, add=T)
}

# Mark destabilizing outside range
m[] = NA
si = which((subst$ddG_glob > grad_range_max & !is.na(subst$rank)) | (subst$ddG_glob > grad_range_max & subst$eff=="destab"))
if (length(si) > 0) {
    m[cbind(aa2int[substr(subst[si,'taa'],1,1)], subst[si,'resi']-first_resn+1)] = 1
    image(m, col=col_destab, add=T)
}

# Gradient color of low uncertainty subst
m[] = NA
m[cbind(aa2int[substr(subst_sig_nogap$taa,1,1)], subst_sig_nogap$resi-first_resn+1)] = subst_sig_nogap$ddG_glob
image(m, zlim=c(grad_range_min,grad_range_max), col=col_grad, breaks=col_breaks, add=T)

# quartz.save("heatmap.png", type="png")

# quartz(width=1.3, height=5)
par(mar=c(3,1,1,3)+.1)
image(t(col_breaks), zlim=c(grad_range_min,grad_range_max), col=col_grad, breaks=col_breaks, xaxt="n", yaxt="n")
# n = grad_range_max - (grad_range_min-1)
n = 5
axis(4, seq(0, n, length=n)/n, seq(grad_range_min,grad_range_max,length.out=n), las=2)
mtext("kJ/mol", 1, 1, cex=.8)
# quartz.save("heatmap_scalebar.png", type="png")
dev.off()

