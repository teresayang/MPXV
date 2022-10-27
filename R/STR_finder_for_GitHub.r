r=c(27,7,24,16) #TAIWAN 1
r=c(20,7,30,16) #TAIWAN 2
r=c(19,7,24,15) #巴西

STR_finder=function(r, L=NULL, out_dir="J:/monkeypox/STR_finder/", file_name="matching_output_n_1798"){ 
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

