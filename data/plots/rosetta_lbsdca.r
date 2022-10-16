load("../../../gmma/release/nonreleased_output/gmma_result.rda")

options(width=180)
pcol = c("signal","active","inactive","init_m","init_ddG","gmma","ddG_glob","stderr_meas","stderr_subfit_est","rank","eff","rosetta","dca")

nmut = length(mutant[,1])
nsubst = length(subst[,1])

si = which(! is.na(subst$rank))

# Convert to kJ/mol
cal2j = 8.314/1.985
subst_sig$ddG_glob = subst_sig$ddG_glob * cal2j
subst$ddG_glob = subst$ddG_glob * cal2j


# Rosetta, Oct. 20, 2020
ros = read.csv("rosetta_edf106_ddg.csv")
wt_aa = substr(ros$variant, 1, 1)
to_aa = substr(ros$variant, nchar(ros$variant), nchar(ros$variant))
resi = substr(ros$variant, 2, nchar(ros$variant)-1)
# numbering is GMMA+1
rownames(ros) = paste0(wt_aa, as.numeric(resi)-1, to_aa)

subst$rosetta = ros[rownames(subst), "Rosetta_ddg_score"]

ros_pear_all = cor(subst$ddG_glob, subst$rosetta, method="pearson", use="complete.obs")
ros_pear = cor(subst[si,"ddG_glob"], subst[si,"rosetta"], method="pearson", use="complete.obs")
ros_spear = cor(subst[si,"ddG_glob"], subst[si,"rosetta"], method="spearman", use="complete.obs")
print(sprintf("Rosetta correlation: %.3f pearson, %.3f spearman",ros_pear,ros_spear))

i = which(! is.na(subst$rank))
i_top_gmma = i[order(subst[i,"ddG_glob"])]

# Plot Rosetta
quartz(height=6, width=6)
plot(subst$ddG_glob, subst$rosetta, xlab="GMMA stability [kJ/mol]", ylab="Rosetta [REU]", pch=16)
points(subst[si,"ddG_glob"], subst[si,"rosetta"], pch=16, cex=.6, col=2)
legend("topleft", c(sprintf("Pearson %.3f",ros_pear_all),sprintf("Pearson %.3f",ros_pear)), pch=16, col=c(1,2))
quartz.save("cor_rosetta.png", type="png")

print("Rosetta top 15 (among subst with accurate GMMA stability)")
i = which(! is.na(subst$rank))
i_top_ros = i[order(subst[i,"rosetta"])]
print(subst[i_top_ros[1:15],c("signal","active","inactive","ddG_glob","rank","rosetta")])

print("Summary of Rosetta scores (subst. w. accurate GMMA est.)")
summary(subst[i,"rosetta"])

print(sprintf("Number of GMMA stabilizing within Rosetta top 15: %d",sum(subst[i_top_ros[1:15],"ddG_glob"] < 0.0)))

for (mut in c("V83L","P92S","P92T","E54V","T57I","L87F","L55V")) {
    rank = which(rownames(subst[i_top_ros,])==mut)
    print(sprintf("Rosetta %s rank %d and score %.2f",mut,rank,subst[mut,"rosetta"]))
}


# lbsDCA
dca = read.table("trx_edf106_msa_variants_scores_prismdb.dat", header=T)
rownames(dca) = dca$variant

subst$dca = dca[rownames(subst), "score"]

dca_pear_all = cor(subst$ddG_glob, subst$dca, method="pearson", use="complete.obs")
dca_pear = cor(subst[si,"ddG_glob"], subst[si,"dca"], method="pearson", use="complete.obs")
dca_spear = cor(subst[si,"ddG_glob"], subst[si,"dca"], method="spearman", use="complete.obs")
print(sprintf("lbsDCA correlation: %.3f pearson, %.3f spearman",dca_pear,dca_spear))

# plot lbsDCA
quartz(height=6, width=6)
plot(subst$ddG_glob, subst$dca, xlab="GMMA stability [kJ/mol]", ylab="lbsDCA [score]", pch=16)
points(subst[si,"ddG_glob"], subst[si,"dca"], pch=16, cex=.6, col=2)
legend("topleft", c(sprintf("Pearson %.3f",dca_pear_all),sprintf("Pearson %.3f",dca_pear)), pch=16, col=c(1,2))
quartz.save("cor_lbsDCA.png", type="png")

print("lbsDCA top 15 (among subst with accurate GMMA stability)")
i = which(! is.na(subst$rank))
i_top_dca = i[order(subst[i,"dca"])]
print(subst[i_top_dca[1:15],pcol])

print("Summary of Rosetta scores (subst. w. accurate GMMA est.)")
summary(subst[i,"dca"])

print(sprintf("Number of GMMA stabilizing within lbsDCA top 15: %d",sum(subst[i_top_dca[1:15],"ddG_glob"] < 0.0)))

for (mut in c("V83L","P92S","P92T","E54V","T57I","L87F","L55V")) {
    rank = which(rownames(subst[i_top_dca,pcol])==mut)
    print(sprintf("lbsDCA %s rank %d, score %.4f",mut,rank,subst[mut,"dca"]))
}


# exclude strand M(K)51-D56 that is in contact with NT fusion to CPOP.
si2 = which(subst$resi > 55 & subst$resi != 79 & ! is.na(subst$rank))
ros_pear2 = cor(subst[si2,"ddG_glob"], subst[si2,"rosetta"], method="pearson", use="complete.obs")
print(sprintf("Pearson when excluding %d of %d variants with resi<=55 or resi==79: %.2f",length(si2),length(si),ros_pear2))
dca_pear2 = cor(subst[si2,"ddG_glob"], subst[si2,"dca"], method="pearson", use="complete.obs")
print(sprintf("Pearson when excluding %d of %d variants with resi<=55 or resi==79: %.2f",length(si2),length(si),dca_pear2))

quartz(height=6, width=6)
plot(subst[si,"ddG_glob"], subst[si,"rosetta"], xlab="GMMA stability [kJ/mol]", ylab="Rosetta [REU]", pch=16, main="Without terminals")
points(subst[si2,"ddG_glob"], subst[si2,"rosetta"], pch=16, cex=.6, col=2)
legend("topleft", c(sprintf("Pearson %.3f",ros_pear),sprintf("Pearson %.3f, resi>56",ros_pear2)), pch=16, col=c(1,2))

quartz(height=6, width=6)
plot(subst[si,"ddG_glob"], subst[si,"dca"], xlab="GMMA stability [kJ/mol]", ylab="lbsDCA [score]", pch=16, main="Without terminals")
points(subst[si2,"ddG_glob"], subst[si2,"dca"], pch=16, cex=.6, col=2)
legend("topleft", c(sprintf("Pearson %.3f",dca_pear),sprintf("Pearson %.3f, resi>56",dca_pear2)), pch=16, col=c(1,2))


##########################
# Paper plot
##########################
require("TeachingDemos")
si = which(! is.na(subst$rank))
# expect that accurate GMMA estimates have stability and conservation scores assigned
stopifnot( sum( is.na( subst[si,"rosetta"] ) ) == 0 )
stopifnot( sum( is.na( subst[si,"dca"] ) ) == 0 )

# top 15 according to different methods
top15 = rownames(subst_sig[1:15,])
ros_sii = order(subst[si,"rosetta"], decreasing=F)
rtop15 = rownames(subst[si[ros_sii],])[1:15]
dca_sii = order(subst[si,"dca"], decreasing=F)
dtop15 = rownames(subst[si[dca_sii],])[1:15]

# value between rank 15 and 16 that separates the top15 region
gmma_cut = mean(subst_sig[c(15,16),"ddG_glob"])
ros_cut = mean(subst[si[ros_sii[c(15,16)]],"rosetta"])
dca_cut = mean(subst[si[dca_sii[c(15,16)]],"dca"])

# quartz(height=10, width=5)
tiff("silico_1200dpi.tiff", width=8, heigh=16, units="cm", res=1200, pointsize=9, compression="lzw")
# jpeg("silico_600dpi.jpg", width=8, heigh=16, units="cm", res=600, pointsize=9, quality=90)
par(mfcol=c(2,1), mar=c(5,4,1,2)+.1)

plot(subst[si,"ddG_glob"], subst[si,"rosetta"], xlab="GMMA stability [kJ/mol]", ylab="Rosetta [REU]", pch=16)
# points(subst[top15,"ddG_glob"], subst[top15,"rosetta"], pch=16, cex=.6, col=2)
abline(v=gmma_cut, h=ros_cut, lty=2)
mark=c("T57I", "E54Y", "E54V")
shadowtext(subst[mark,"ddG_glob"], subst[mark,"rosetta"], mark, pos=c(2,3,4), cex=.8, xpd=T, offset=.3, col=1, bg="white", font=2, r=0.1)
# shadowtext(subst[mut,'ddG_glob'], subst[mut,'signal'], mut_name, pos=pos, cex=.7, col=1, bg="white", font=2, r=0.1)
# legend("topleft", c("Substitution effects","GMMA top 15"), pch=16, col=c(1,2))
text(-23, 9.3, "a", font=2, cex=1.2, xpd=T)

plot(subst[si,"ddG_glob"], subst[si,"dca"], xlab="GMMA stability [kJ/mol]", ylab="lbsDCA [score]", pch=16)
# points(subst[top15,"ddG_glob"], subst[top15,"dca"], pch=16, cex=.6, col=2)
abline(v=gmma_cut, h=dca_cut, lty=2)
mark=c("M51R", "M51K", "M51T", "L55V")
shadowtext(subst[mark,"ddG_glob"], subst[mark,"dca"], mark, pos=c(3,4,4,4), cex=.8, xpd=T, offset=.3, col=1, bg="white", font=2, r=0.1)
# legend("topleft", c(sprintf("Pearson %.3f",dca_pear_all),sprintf("Pearson %.3f",dca_pear)), pch=16, col=c(1,2))
text(-23, 0.275, "b", font=2, cex=1.2, xpd=T)

dev.off()
