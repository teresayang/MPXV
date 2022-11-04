r=c(27,7,24,16) #TAIWAN 1
r=c(20,7,30,16) #TAIWAN 2
r=c(19,7,24,15) #å·´è¥¿



#' %% ~~function to do ... ~~ Tracking the origin of monkeypox virus (MPXV) by
#' using Variable Number Tandem Repeats (VNTR) markers
#' 
#' %% ~~ A concise (1-5 lines) description of what the function does. ~~ The
#' function \code{VNTRtracker} calculates the MPXV strain-strain genetic
#' distance to track the origin of a given MPXV by using VNTRs.
#' 
#' %% ~~ If necessary, more details than the description above ~~ The four VNTR
#' markers nt133095[A25R<e2><88>’A26L] (T), nt150560[A47R<e2><88>’A49R]
#' (TATGATGGA), nt173267[B14R<e2><88>’B15L] (AT), and nt179074[B18R<e2><88>’B19R]
#' (ATATACATT) are proposed to track the original of MPXV strains. Based on
#' these four VNTRs, we calculate three distances (ploymorphism information
#' content, entropy, and length of repetitive unit (see Yang et al. (2022))
#' between a query MPXV strain and the strains in our database.
#' 
#' @param r A vector of copies of VNTR of MPXV.
#' @param out_dir Directory name of the output files.
#' @param file_name Filename of the output.
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...  The output contains multiple columns, including accession
#' ID, country name, lineage, release date, clade, lineage, collection date,
#' and three distances mentioned above.
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~ Yang, H.-C.
#' et al (2022) \emph{Monkeypox genome contains variable number tandem repeats
#' enabling accurate virus tracking}.
#' @examples
#' 
#'   r=c(27,7,24,16)
#' 
#'   VNTRtracker(r,
#'               out_dir="C:/",
#'               file_name="Tracking_output")
#' 
VNTRtracker=function(r, L=NULL, out_dir="/", file_name="Tracking_output"){
	if(is.null(L)){L=c(1,9,2,9)}#nt 133095 150554 173287 179074
	da=read.csv("J:/monkeypox/VNTR_listbaseonly_TRUE.csv",check.names=F)
	information=read.csv("J:/monkeypox/STR/Metadata_forCSS_n1798.csv",check.names=F)
	information=information[match(da$ID,information$Accession),]
	nts=c("nt133095_T", "nt150554_TATGATGGA", "nt173267_AT", "nt179074_ATATACATT")
	da_r=da[,nts]

	maxabs=apply(rbind(da_r,r),2,function(x) (max(x)-min(x)))
	miss_points=apply(da_r,1,function(x) (abs(x-r)/maxabs))
	STR_dist=apply(da_r,2,function(x) table(x)/sum(table(x)))

	weight_L=L/sum(L)

	weight=1-sapply(STR_dist,function(x) sum(x^2)+sum(x^2)^2-sum(x^4))
	weight=weight/sum(weight)

	weight_entropy=sapply(STR_dist,function(x) sum(x*log(x))/log(1/length(x)))
	weight_entropy=weight_entropy/sum(weight_entropy)

	Dis=apply(miss_points,2,function(x) {mean(x*weight)})
	Dis_L=apply(miss_points,2,function(x) {mean(x*weight_L)})
	Dis_entropy=apply(miss_points,2,function(x) {mean(x*weight_entropy)})

	sel=order(Dis)
	Distance_PIC=Dis[sel]
	Distance_L=Dis_L[sel]
	Distance_entropy=Dis_entropy[sel]

	output=cbind(information[order(Dis),],Distance_PIC,Distance_L,Distance_entropy)

	write.csv(output,file=paste(out_dir,file_name,".csv",sep=""))
}













