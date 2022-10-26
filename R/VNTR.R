
library(seqinr)
library(stringr)
library(Biostrings)
library(stringi)

ref_fa <- read.fasta("D:/ctyang/Monkeypox/MA001.fasta",as.string = T)



VNTR_sub <- function(data, STR=STR, match_s=match_s, mismatch_s=mismatch_s,
                     regionStart=regionStart, regionEnd=regionEnd,baseonly = baseonly,VNTRoutput=VNTRoutput){


  str <- paste(rep(STR,100),collapse = "")
  match <- rep(NA,length(data))
  mismatch <- rep(NA,length(data))
  score <- rep(NA,length(data))
  r <- rep(NA,length(data))
  str_ind <- rep(NA,length(data))
  indel <- rep(NA,length(data))
  STR_align <- list()
  l_score <- rep(NA,length(data))
  r_score <- rep(NA,length(data))
  probe_add_l <- rep(NA,length(data))
  probe_add_r <- rep(NA,length(data))
  start_pos <- rep(NA,length(data))
  r_start <- rep(NA,length(data))
  r_end <- rep(NA,length(data))
  ATCG_p <- rep(NA,length(data))

  s_fa <- DNAString(as.character(substr(ref_fa,regionStart-150,regionEnd+150)))
  s2 <- DNAString(str)
  mat <- nucleotideSubstitutionMatrix(match = match_s, mismatch = mismatch_s, baseOnly = F)
  mat[!mat%in%c(match_s,mismatch_s)] <- match_s
  reflAlign <- pairwiseAlignment(s_fa, s2, substitutionMatrix = mat,gapOpening = -mismatch_s, gapExtension = 2,type="local")

  str_pos_ref <- str_locate_all(gsub("-","",reflAlign@subject),STR)[[1]]

  gap_l <- rep(0,nchar(reflAlign@subject))
  gap_l[str_locate_all(reflAlign@subject,"-")[[1]][,1]] <- 1
  gap_l<- cumsum(gap_l)
  #
  a <- str_pos_ref[1,1]+gap_l[str_pos_ref[1,1]]
  b <- str_pos_ref[nrow(str_pos_ref),2]+gap_l[str_pos_ref[nrow(str_pos_ref),2]]
  b <- b-str_count(substr(reflAlign@pattern,a,b),"-")

  nt <- str_locate_all(toupper(ref_fa),gsub("-","",as.character(reflAlign@pattern)))[[1]][1]



  seq_ref_l <- substr(s_fa,reflAlign@pattern@range@start-25,reflAlign@pattern@range@start-1)
  seq_ref_r <-  substr(s_fa,reflAlign@pattern@range@start+b,reflAlign@pattern@range@start+b+24)





  for(i in 1:length(data)){

    raw_fa <- toupper(data[i])


    seq_raw <- DNAString(raw_fa)
    mat <- nucleotideSubstitutionMatrix(match = 2, mismatch = -4, baseOnly = F)
    ##
    mat[!mat%in%c(2,-4)] <- 1
    Align_l <- pairwiseAlignment(seq_ref_l, seq_raw ,type="global-local", substitutionMatrix = mat,gapOpening = 4, gapExtension = 1)

    Align_r <- pairwiseAlignment(seq_ref_r, seq_raw ,type="global-local", substitutionMatrix = mat,gapOpening = 4, gapExtension = 1)

    if(str_count(Align_l@subject,"N")>nchar(Align_l@subject)*0.5){
      seq_ref_l_1 <- substr(s_fa,max(1,reflAlign@pattern@range@start-150),reflAlign@pattern@range@start-1)
      Align_l <- pairwiseAlignment(seq_ref_l_1, seq_raw ,type="global-local", substitutionMatrix = mat,gapOpening = 4, gapExtension = 1)
      probe_add_l[i] <- 1
    }

    if(str_count(Align_r@subject,"N")>nchar(Align_r@subject)*0.5){
      seq_ref_r_1 <- substr(s_fa,reflAlign@pattern@range@start+reflAlign@pattern@range@width,min(nchar(s_fa),reflAlign@pattern@range@start+reflAlign@pattern@range@width+149))
      Align_r <-  pairwiseAlignment(seq_ref_r_1, seq_raw ,type="global-local", substitutionMatrix = mat,gapOpening = 4, gapExtension = 1)
      probe_add_r[i] <- 1
    }

    l_score[i] <- Align_l@score
    r_score[i] <- Align_r@score

    if(str_count(Align_r@subject,"N")>100|str_count(Align_l@subject,"N")>100){
      seq_txt <- ""
    }else{
      seq_txt <- substr(raw_fa,Align_l@subject@range@start+Align_l@subject@range@width,
                        Align_r@subject@range@start-1)

      deg <- setdiff(LETTERS,c("A","T","C","G"))
      seq_txt <- gsub(ifelse(baseonly==T,paste0(deg,collapse = "|"),""),"",seq_txt)



      l_seq <- strsplit(as.character(Align_l@subject), split = "")[[1]]
      l_ref <- strsplit(as.character(Align_l@pattern), split = "")[[1]]


      l_seq <- gsub("-","",paste0(l_seq,collapse =  ""))

      r_seq <- strsplit(as.character(Align_r@subject), split = "")[[1]]
      r_ref <- strsplit(as.character(Align_r@pattern), split = "")[[1]]


      r_seq <- gsub("-","",paste0(r_seq,collapse =  ""))


      txt <- seq_txt
      seq_raw <- DNAString(paste0(substr(raw_fa,1,Align_l@subject@range@start+Align_l@subject@range@width-1),txt,substr(raw_fa,Align_r@subject@range@start,nchar(raw_fa))))

    }


    if(nchar(txt)!=0&nchar(seq_txt)>=nchar(STR)){

      if(nchar(STR)==1){
        STR_align[[i]] <- txt
        start_pos[i] <- Align_l@subject@range@start+Align_l@subject@range@width-1+a

        s1_1 <- txt
        s2_1 <- rep(STR,nchar(txt))
        TF <- table(factor(strsplit(s1_1,"")[[1]]==strsplit(s2_1,"")[[1]], levels = c("TRUE","FALSE")))
        match[i] <- TF["TRUE"]
        mismatch[i] <- TF["FALSE"]
        indel[i] <- 0
        r[i] <- str_count(txt,STR)

        r_start[i] <- 1
        r_end[i] <- nchar(txt)

      }else{
        s1 <- DNAString(txt)
        s2 <- DNAString(str)
        mat <- nucleotideSubstitutionMatrix(match = match_s, mismatch = mismatch_s, baseOnly = F)
        mat[!mat%in%c(match_s,mismatch_s)] <- ifelse(baseonly==T,mismatch_s,match_s)


        globalAlign <- pairwiseAlignment(s1,s2, substitutionMatrix = mat,gapOpening = -mismatch_s, gapExtension = 2,type="overlap")



        s1_1 <- globalAlign@pattern
        s2_1 <- globalAlign@subject

        str_pos <- str_locate_all(gsub("-","",s2_1),STR)[[1]]


        if(nrow(str_pos)!=0){
          gap_l <- rep(0,nchar(s2_1))
          gap_l[str_locate_all(s2_1,"-")[[1]][,1]] <- 1
          gap_l<- cumsum(gap_l)

          a <- str_pos[1,1]+gap_l[str_pos[1,1]]
          b <- str_pos[nrow(str_pos),2]+gap_l[str_pos[nrow(str_pos),2]]

          ATCG_p[i] <- paste(sapply(1:nrow(str_pos), function(x) str_count(substr(s1_1,str_pos[x,1]+gap_l[str_pos[x,1]],str_pos[x,2]+gap_l[str_pos[x,2]]),"A|T|C|G|-"))-sapply(1:nrow(str_pos), function(x) str_count(substr(s2_1,str_pos[x,1]+gap_l[str_pos[x,1]],str_pos[x,2]+gap_l[str_pos[x,2]]),"-")),collapse = "-")

          s1_1 <- substr(s1_1,a,b)
          s2_1 <- substr(s2_1,a,b)



          STR_align[[i]] <- txt
          start_pos[i] <- Align_l@subject@range@start+Align_l@subject@range@width-1+a
          TF <- table(factor(strsplit(s1_1,"")[[1]]==strsplit(s2_1,"")[[1]], levels = c("TRUE","FALSE")))
          str_loc <- str_locate_all(gsub("-","",as.character(globalAlign@subject)),STR)[[1]]
          r_start[i] <- str_locate(txt,gsub("-","",s1_1))[1]
          r_end[i] <- str_locate(txt,gsub("-","",s1_1))[2]

          match[i] <- TF["TRUE"]
          mismatch[i] <- TF["FALSE"]
          indel[i] <- str_count(s2_1,"-")+str_count(s1_1,"-")
          r[i]<-nrow(str_loc)



        }else{

          r[i]<-0
          STR_align[[i]] <- ""
        }

      }



    }else{
      r[i]<-0
      STR_align[[i]] <- ""
    }


  }


  names(STR_align) <- names(data)

  out <- data.frame(ID=names(data),r=r,match=match,mismatch=mismatch,
                    indel=indel,score=match*match_s+mismatch*mismatch_s,STR=str_ind,
                    l_score=l_score,r_score=r_score,
                    probe_add_l=probe_add_l,probe_add_r=probe_add_r,start_pos=start_pos,
                    r_start=r_start,r_end=r_end,ATCG_p=ATCG_p)
  out$start_pos[out$r==0] <- NA

  if(VNTRoutput==T){
    dir.create("VNTR",showWarnings = F)
    write.csv(out,paste0("VNTR/VNTR_nt",paste(nt,STR,"baseonly",baseonly,sep = "_"),".csv"),row.names = F)
    seqinr::write.fasta(STR_align,names(STR_align), paste0("VNTR/VNTR_nt",paste(nt,STR,"baseonly",baseonly,sep = "_"),".fas"))
  }

  list(paste0("VNTR_nt",paste(nt,STR,"baseonly",baseonly,sep = "_")),out,STR_align)
}

STR_finder=function(r, L=NULL, brief=T, out_dir="J:/monkeypox/STR_finder/", file_name="matching_output_n_1407"){
  if(is.null(L)){L=c(1,9,2,9)}#nt 133095 150554 173287 179074

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
  Distance=Dis[sel]
  Distance_L=Dis_L[sel]
  Distance_entro=Dis_entropy[sel]
  if(brief){
    VNTRoutput=cbind(da[order(Dis),],Distance,Distance_L,Distance_entro)
  }else{
    VNTRoutput=cbind(da_r[order(Dis),],Distance,Distance_L,Distance_entro)
  }
  colnames(VNTRoutput)[colnames(VNTRoutput)%in%c("Distance","Distance_L","Distance_entro")] <- c("Distance (PIC)", "Distance (Length)", "Distance (Entropy)")
  write.csv(VNTRoutput,file=paste(out_dir,file_name,".csv",sep=""),row.names = F)
}
# data <- read.fasta("G:/Monkeypox/n633/data/MPXV_20220726_n633_combined.fasta",as.string = T)
# read.fasta("G:/Monkeypox/n633/data/MPXV_20220726_n633_combined.fasta",as.string = T)




VNTR <- function(data, STR=STR, match_s=match_s, mismatch_s=mismatch_s,
                 regionStart=regionStart, regionEnd=regionEnd,baseonly = T,VNTRoutput=F,finder=F,brief=T){
  if(sum(is.na(as.numeric(c(regionStart,regionEnd))))!=0){
    stop("regionStart or regionEnd should be numeric.")
  }
  if(!all.equal(length(STR),length(regionStart),length(regionEnd ))){
    stop("")
  }
  if(sum(is.na(as.numeric(c(match_s,mismatch_s))))!=0){
    stop("match_s or mismatch_s should be numeric.")
  }
  if(!baseonly%in%c("TRUE","FALSE")){
    stop("baseonly should be TRUE or FALSE. If TRUE, only uses the letters in base alphabet i.e. A,C,G,T.")
  }

  if(!VNTRoutput%in%c("TRUE","FALSE")){
    stop("VNTRoutput should be TRUE or FALSE.")
  }
  out <- list()
  invisible(sapply(1:length(STR), function(x)out[[x]] <<-VNTR_sub(data, STR=STR[x], match_s=match_s, mismatch_s=mismatch_s,regionStart=regionStart[x], regionEnd=regionEnd[x],baseonly = baseonly,VNTRoutput=VNTRoutput)))
  names(out) <- sapply(1:length(out), function(x)out[[x]][[1]])


  dt <- data.frame(ID=names(data), sapply(1:length(out), function(x) out[[x]][[2]][,2]))
  colnames(dt)[-1] <- sapply(1:length(names(out)), function(x)paste(strsplit(names(out)[x],"_")[[1]][2:3],collapse  = "_") )
  if(VNTRoutput==T){
    write.csv(dt,paste0("VNTR/VNTR_list",paste("baseonly",baseonly,sep = "_"),".csv"),row.names = F)
  }
  if(finder==T){
    load("G:/Monkeypox/VNTR_program/VNTR_database_n1407.RData")
    nts = colnames(dt)[-1]
    da_r = da[,nts]
    L <- sapply(1:length(STR), function(x)nchar(STR[x]))
    dir.create("VNTR/finder",showWarnings = F)
    invisible(sapply(1:nrow(dt), function(x)STR_finder(r = unlist(dt[x,2:5]),L=L,brief=brief,out_dir=getwd(),file_name= paste0("VNTR/finder","/",dt$ID[x]))))
  }
}

