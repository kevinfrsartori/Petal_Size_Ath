##########
#
# Custom calc_ehh function from rehh package
#
##########


my_calc_ehh<-function (haplohh, mrk, limhaplo = 2, limehh = 0.05, plotehh = TRUE, 
                       main_leg = "EHH plot") 
{
  if (mrk < 1 | mrk > haplohh@nsnp) {
    stop(paste("Focal snp index must be between", 1, "and", 
               haplohh@nsnp))
  }
  if (limhaplo < 2) {
    stop("limhaplo must be >1")
  }
  if (limehh < 0 | limehh > 1) {
    stop("limehh must be between 0 and 1")
  }
  ehh <- matrix(0, nrow = haplohh@nsnp, ncol = 2)
  nhaplo_eval <- matrix(0, nrow = haplohh@nsnp, ncol = 2)
  ihh <- rep(0, 2)
  res.ehh <- .C("r_ehh", Rdata = as.integer(haplohh@haplo), 
                number_SNPs = as.integer(haplohh@nsnp), number_chromosomes = as.integer(haplohh@nhap), 
                focal_SNP = as.integer(mrk), map = as.double(haplohh@position), 
                number_haplotypes = as.integer(nhaplo_eval), EHH = as.double(ehh), 
                IHH = as.double(ihh), min_number_haplotypes = as.integer(limhaplo), 
                min_EHH = as.double(limehh))
  ehh = matrix(res.ehh$EHH, 2, haplohh@nsnp, byrow = T)
  nhaplo_eval = matrix(res.ehh$number_haplotypes, 2, haplohh@nsnp, 
                       byrow = T)
  #rownames(ehh) = rownames(nhaplo_eval) = c("Anc. Allele","Der. Allele")
  A1<-map[mrk,4]#
  A2<-map[mrk,5]#
  rownames(ehh) = rownames(nhaplo_eval) = c(A1,A2)#
  colnames(ehh) = colnames(nhaplo_eval) = haplohh@snp.name
  ihh = res.ehh$IHH
  #names(ihh) = c("Anc. Allele", "Der. Allele")
  names(ihh) = c(A1,A2)
  if (plotehh) {
    sel_reg <- (colSums(nhaplo_eval) > 0)
    if (sum(sel_reg) > 0) {
      matplot(haplohh@position[sel_reg]/1000, t(ehh[, sel_reg]), 
              col = c("red", "blue"), lty = 1, type = "l", 
              main = paste(main_leg,"- Chrm",map$V2[mrk] ,"Pos",map$V3[mrk]  ), bty = "n", xlab = "Position (K bp)", 
              ylab = "EHH",las=1)
      abline(v = haplohh@position[mrk]/1000, lty = 2)
      #smartlegend("left", "top", c("Anc. Allele", "Der. Allele"),col = c("red", "blue"), bty = "n", lty = 1)
      legend("left", "top", c(A1,A2),col = c("red", "blue"), bty = "n", lty = 1)
    }
  }
  return(list(ehh = ehh, nhaplo_eval = nhaplo_eval, freq_all1 = nhaplo_eval[1, 
                                                                            mrk]/sum(nhaplo_eval[, mrk]), ihh = ihh))
}