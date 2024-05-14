#####
# homemade function to compute standardized iH metric for a list of snps that do not have a proper estimation of ancestry
# after having computed iHs whole genome
# KS - 2023-10-17
#####


ihh2ihs_ks<-function (res_ihh, mis_ihh, freqbin = 0.025, minmaf = 0.05) 
{
  # 1 - first redo the summary_class from whole genome data
    freq_class = seq(minmaf, 1 - minmaf, freqbin)
    summary_class = matrix(0, length(freq_class) - 1, 5)
    colnames(summary_class) = c("Freq Lim Inf", "Freq Lim Sup", "Size", "Mean iHH", 
                                "SD iHH")
    #ihs = log(res_ihh[, 4]/res_ihh[, 5])
    #res_ihh$uniHS[res_ihh$uniHS == "Inf" | res_ihh$uniHS == "-Inf"] = NA
    for (c in 1:(length(freq_class) - 1)) {
      lim_inf = freq_class[c]
      lim_sup = freq_class[c + 1]
      mrk_sel = (res_ihh[, 3] >= lim_inf & res_ihh[, 3] < 
                   lim_sup)
      tmp_ihs = res_ihh$uniHS[mrk_sel]
      tmp_nmrk = sum(mrk_sel)
      if (tmp_nmrk < 10) {
        warning(paste("Size of Allele Frequency Class: ", 
                      lim_inf, "-", lim_sup, " <10: You should probably increase freqbin\n", 
                      sep = ""))
      }
      tmp_mean = mean(tmp_ihs, na.rm = T)
      tmp_sd = sd(tmp_ihs, na.rm = T)
      summary_class[c, 1] = lim_inf
      summary_class[c, 2] = lim_sup
      summary_class[c, 3] = tmp_nmrk
      summary_class[c, 4] = tmp_mean
      summary_class[c, 5] = tmp_sd
    }
  # 2 - compute uniHs and iHs for new data
    mis_ihh$uniHS  = log(mis_ihh[, 4]/mis_ihh[, 5])
    mis_ihh$uniHS[mis_ihh$uniHS == "Inf" | mis_ihh$uniHS == "-Inf"] = NA
    mis_ihh$iHS <- NA
    mis_ihh$Pvalue <- NA 
    for (c in 1:dim(summary_class)[1]) {
      mrk_sel = (mis_ihh[, 3] >= summary_class[c,1]  & mis_ihh[, 3] < summary_class[c,2])
      tmp_ihs = mis_ihh$uniHS[mrk_sel]
      tmp_nmrk = sum(mrk_sel)
      if (tmp_nmrk > 0) {
        mis_ihh$iHS[mrk_sel] = (mis_ihh$uniHS[mrk_sel] - summary_class[c,4])/summary_class[c,5]
      }
    }
  # 3 - compute pvalue
    tmp_pval = -1 * log10(1 - 2 * abs(pnorm(c(mis_ihh$iHS,res_ihh$iHS)) - 0.5))
    tmp_pval2 = tmp_pval
    tmp_pval2[tmp_pval2 == "Inf"] = NA
    tmp_pval[tmp_pval == "Inf"] = max(tmp_pval2, na.rm = TRUE) + 1
    mis_ihh$Pvalue<-tmp_pval[1:dim(mis_ihh)[1]]
    return(mis_ihh)
}
