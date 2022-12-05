information <- readRDS("data/VNTR_n1797.rds")
































#' Tracking the origin of monkeypox virus (MPXV) by using Variable Number
#' Tandem Repeats (VNTR) markers
#'
#' The function \code{VNTRtracker} calculates the MPXV strain-strain genetic
#' distance to track the origin of a given MPXV by using VNTRs.
#'
#' The four VNTR markers nt133095[A25R - A26L] (T), nt150560[A47R - A49R]
#' (TATGATGGA), nt173267[B14R - B15L] (AT), and nt179074[B18R - B19R]
#' (ATATACATT) are proposed to track the original of MPXV strains. Based on
#' these four VNTRs, we calculate three distances (ploymorphism information
#' content, entropy, and length of repetitive unit (see Yang et al. (2022))
#' between a query MPXV strain and the strains in our database.
#'
#' @param r A vector of copies of VNTR of MPXV.
#' @param file_name Filename of the output.
#' @return The output contains multiple columns, including accession ID,
#' country name, collection date, clade, lineage, and three distances mentioned
#' above.
#' @references Yang, H.-C. et al (2022) Monkeypox genome contains variable
#' number tandem repeats enabling accurate virus tracking.
#' @examples
#'
#' r=c(27,7,24,16)
#'
#' ## The output is a csv file with the file name given by the parameter "file_name".
#' ## In this file, each row represents a strain in our database, and column variables contain
#' ## accession ID, country name, collection date, clade, lineage, and three distances (see Yang et al. (2022))
#' ## between a query MPXV strain and a strains in our database.
#' ## The rank of strains are ordered by ploymorphism information content distance (see Yang et al. (2022)) from low to high.
#' VNTRtracker(r, file_name="Tracking_output")
#'
#' @importFrom utils write.csv
#' @export VNTRtracker
VNTRtracker=function(r, file_name="Tracking_output"){
  L=c(1,9,2,9)#nt 133095 150554 173287 179074

  nts=c("nt133095_T", "nt150554_TATGATGGA", "nt173267_AT", "nt179074_ATATACATT")
  da_r=information[,nts]

  maxabs=apply(rbind(da_r,r),2,function(x) (max(x,na.rm=T)-min(x,na.rm=T)))
  miss_points=apply(da_r,1,function(x) (abs(x-r)/maxabs))
  STR_dist=apply(da_r,2,function(x) table(x)/sum(table(x)))

  weight_L=L/sum(L)

  weight=1-sapply(STR_dist,function(x) sum(x^2)+sum(x^2)^2-sum(x^4))
  weight=weight/sum(weight)

  weight_entropy=sapply(STR_dist,function(x) sum(x*log(x))/log(1/length(x)))
  weight_entropy=weight_entropy/sum(weight_entropy)

  Dis=apply(miss_points,2,function(x) {sum(x*weight,na.rm=T)})
  Dis_L=apply(miss_points,2,function(x) {sum(x*weight_L,na.rm=T)})
  Dis_entropy=apply(miss_points,2,function(x) {sum(x*weight_entropy,na.rm=T)})

  sel=order(Dis)
  Distance_PIC=Dis[sel]
  Distance_L=Dis_L[sel]
  Distance_entropy=Dis_entropy[sel]

  output=cbind(information[order(Dis),1:5],Distance_PIC,Distance_L,Distance_entropy)

  dir.create("output",showWarnings = F)
  dir.create("output/VNTRTracker",showWarnings = F)
  write.csv(output,file=paste0("output/VNTRTracker/",file_name,".csv"),row.names = F)
  print(paste0("Tracking VNTR is done and the results are saved in ", getwd(), "/output/VNTRTracker"))
}
