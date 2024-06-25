###
#
###
unit_object <- read.csv("analyses/efficiency/unit_object_level_2v13v16_Unit_10pct.tab", sep = '\t')

unit_object<-subset(unit_object, (unit_object$DGY1657_ms_pt_median != 0 &
                                    unit_object$DGY1657_rpf_pt_median != 0 &
                                    unit_object$DGY1726_ms_pt_median != 0 &
                                    unit_object$DGY1726_rpf_pt_median != 0 &
                                    unit_object$DGY1735_ms_pt_median != 0 &
                                    unit_object$DGY1735_rpf_pt_median != 0 &
                                    unit_object$DGY1741_ms_pt_median != 0 &
                                    unit_object$DGY1741_rpf_pt_median != 0 &
                                    unit_object$DGY1743_ms_pt_median != 0 &
                                    unit_object$DGY1743_rpf_pt_median != 0))

#DGY1726
unit_object$anc_teff<-log2(unit_object$DGY1657_ms_pt_median/unit_object$DGY1657_rpf_pt_median)
unit_object$evo_teff<-log2(unit_object$DGY1726_ms_pt_median/unit_object$DGY1726_rpf_pt_median)

pdf('figures/_DGY1726_ms_rpf_l2_peff.pdf')

plot((unit_object$anc_teff), (unit_object$evo_teff), col=rgb(0,0,0,0.1), pch = 16)
fit <- lm(unit_object$evo_teff ~ unit_object$anc_teff)
abline(fit)
summary(fit)
unit_object$sres_DGY1726<-rstandard(fit)

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1726) >= 2 & unit_object$DGY1726_ms_rpf_fdr <= 0.05) &
                              (unit_object$DGY1726_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(1,0,1,0.5), pch = 16)
print(c("Both",length(sig$X)))

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1726) >= 2 & unit_object$DGY1726_ms_rpf_fdr > 0.05) &
                              (unit_object$DGY1726_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(0,1,0,0.5), pch = 16)
print(c("Sres",length(sig$X)))

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1726) < 2 & unit_object$DGY1726_ms_rpf_fdr <= 0.05) & 
                              (unit_object$DGY1726_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(0,0,1,0.5), pch = 16)
print(c("FDR",length(sig$X)))

dev.off()

#DGY1735
unit_object$anc_teff<-log2(unit_object$DGY1657_ms_pt_median/unit_object$DGY1657_rpf_pt_median)
unit_object$evo_teff<-log2(unit_object$DGY1735_ms_pt_median/unit_object$DGY1735_rpf_pt_median)

pdf('figures/_DGY1735_ms_rpf_l2_peff.pdf')

plot((unit_object$anc_teff), (unit_object$evo_teff), col=rgb(0,0,0,0.1), pch = 16)
fit <- lm(unit_object$evo_teff ~ unit_object$anc_teff)
abline(fit)
summary(fit)
unit_object$sres_DGY1735<-rstandard(fit)

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1735) >= 2 & unit_object$DGY1735_ms_rpf_fdr <= 0.05) &
                              (unit_object$DGY1735_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(1,0,1,0.5), pch = 16)
print(c("Both",length(sig$X)))

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1735) >= 2 & unit_object$DGY1735_ms_rpf_fdr > 0.05) &
                              (unit_object$DGY1735_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(0,1,0,0.5), pch = 16)
print(c("Sres",length(sig$X)))

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1735) < 2 & unit_object$DGY1735_ms_rpf_fdr <= 0.05) & 
                              (unit_object$DGY1735_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(0,0,1,0.5), pch = 16)
print(c("FDR",length(sig$X)))

dev.off()

#DGY1741
unit_object$anc_teff<-log2(unit_object$DGY1657_ms_pt_median/unit_object$DGY1657_rpf_pt_median)
unit_object$evo_teff<-log2(unit_object$DGY1741_ms_pt_median/unit_object$DGY1741_rpf_pt_median)

pdf('figures/_DGY1741_ms_rpf_l2_peff.pdf')

plot((unit_object$anc_teff), (unit_object$evo_teff), col=rgb(0,0,0,0.1), pch = 16)
fit <- lm(unit_object$evo_teff ~ unit_object$anc_teff)
abline(fit)
summary(fit)
unit_object$sres_DGY1741<-rstandard(fit)

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1741) >= 2 & unit_object$DGY1741_ms_rpf_fdr <= 0.05) &
                              (unit_object$DGY1741_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(1,0,1,0.5), pch = 16)
print(c("Both",length(sig$X)))

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1741) >= 2 & unit_object$DGY1741_ms_rpf_fdr > 0.05) &
                              (unit_object$DGY1741_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(0,1,0,0.5), pch = 16)
print(c("Sres",length(sig$X)))

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1741) < 2 & unit_object$DGY1741_ms_rpf_fdr <= 0.05) & 
                              (unit_object$DGY1741_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(0,0,1,0.5), pch = 16)
print(c("FDR",length(sig$X)))

dev.off()

#DGY1743
unit_object$anc_teff<-log2(unit_object$DGY1657_ms_pt_median/unit_object$DGY1657_rpf_pt_median)
unit_object$evo_teff<-log2(unit_object$DGY1743_ms_pt_median/unit_object$DGY1743_rpf_pt_median)

pdf('figures/_DGY1743_ms_rpf_l2_peff.pdf')

plot((unit_object$anc_teff), (unit_object$evo_teff), col=rgb(0,0,0,0.1), pch = 16)
fit <- lm(unit_object$evo_teff ~ unit_object$anc_teff)
abline(fit)
summary(fit)
unit_object$sres_DGY1743<-rstandard(fit)

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1743) >= 2 & unit_object$DGY1743_ms_rpf_fdr <= 0.05) &
                              (unit_object$DGY1743_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(1,0,1,0.5), pch = 16)
print(c("Both",length(sig$X)))

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1743) >= 2 & unit_object$DGY1743_ms_rpf_fdr > 0.05) &
                              (unit_object$DGY1743_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(0,1,0,0.5), pch = 16)
print(c("Sres",length(sig$X)))

sig <- subset(unit_object, ((abs(unit_object$sres_DGY1743) < 2 & unit_object$DGY1743_ms_rpf_fdr <= 0.05) & 
                              (unit_object$DGY1743_ms_rpf_signal_to_noise == "True")))
points((sig$anc_teff), (sig$evo_teff), col=rgb(0,0,1,0.5), pch = 16)
print(c("FDR",length(sig$X)))

dev.off()

#
write.csv(unit_object,"analyses/efficiency/unit_object_level2_ms_rpf_16_Unit_Sres_fdr.csv", row.names = FALSE)


#

