flanking_seq_match <- function(seq,pattern,start,end,max.mismatch,fixed,ext){

  df <- sapply(1:length(end),function(x){
    m <- matchPattern(pattern,seq[start[x]:end[x]],
                      max.mismatch=max.mismatch, min.mismatch=0,
                      with.indels=TRUE, fixed=fixed,
                      algorithm="auto")
    if(length(m@ranges@start)==0){
      c(0,0)
    }else{
      cbind(start[x]+m@ranges@start-1,m@ranges@width)
    }


  })


  if(is.matrix(df)==TRUE){
    df <-  as.data.frame(t(df))
  }else{
    df <- as.data.frame(do.call("rbind",df))
  }

  colnames(df) <- c("start","width")
  df$deg <- 0

  df <- df[df$width!=0,]


  if(ext==1){
    df_N <- sapply(1:length(end),function(x){
      m <- matchPattern(paste(rep("N",nchar(pattern)),collapse = ""),seq[start[x]:end[x]],
                        max.mismatch=max.mismatch, min.mismatch=0,
                        with.indels=TRUE, fixed=T,
                        algorithm="auto")
      if(length(m@ranges@start)==0){
        c(0,0)
      }else{
        cbind(start[x]+m@ranges@start-1,m@ranges@width)
      }

    })

    if(is.matrix(df_N)==TRUE){
      df_N <-  as.data.frame(t(df_N))
    }else{
      df_N <- as.data.frame(do.call("rbind",df_N))
    }

    colnames(df_N) <- c("start","width")

    df <- df[!df$start%in%df_N$start,]
    if(nrow(df)!=0){
      df$deg <- sapply(1:nrow(df),function(x)    str_count(as.character(substr(seq,df$start[x],df$start[x]+df$width[x]-1)),paste0(setdiff(LETTERS,c("A","T","C","G")),collapse = "|")) )
      df <- df[df$width!=df$deg,]
      df <- df[(df$deg/df$width)==min(df$deg/df$width),]
    }
  }



  df
}
ref_fa <- seqinr::read.fasta("inst/extdata/MA001.fasta",as.string = T)

VNTR_sub <- function(data, vntr=vntr,
                     regionStart=regionStart, regionEnd=regionEnd,
                     match_s=match_s, mismatch_s=mismatch_s,
                     baseonly = baseonly,fseqExtend = fseqExtend ,VNTRoutput=VNTRoutput,
                     len_fl = len_fl, mis_prop = mis_prop ){

  max.mismatch = ceiling(len_fl*mis_prop)
  max.mismatch.ext = ceiling(fseqExtend*mis_prop)

  vntr_t <- paste(rep(vntr,100),collapse = "")
  match <- rep(NA,length(data))
  mismatch <- rep(NA,length(data))
  score <- rep(NA,length(data))
  r <- rep(NA,length(data))
  vntr_ind <- rep(NA,length(data))
  indel <- rep(NA,length(data))
  vntr_align <- list()
  probe_add_l <- rep(0,length(data))
  probe_add_r <- rep(0,length(data))
  start_pos <- rep(NA,length(data))
  r_start <- rep(NA,length(data))
  r_end <- rep(NA,length(data))
  ATCG_p <- rep(NA,length(data))
  count_n <- rep(NA,length(data))
  count_deg <- rep(NA,length(data))
  count_len <- rep(NA,length(data))
  vntr_match_norm <- rep(NA,length(data))

  s_fa <- DNAString(as.character(substr(ref_fa,regionStart-fseqExtend,regionEnd+fseqExtend)))
  s2 <- DNAString(vntr_t)
  mat <- nucleotideSubstitutionMatrix(match = match_s, mismatch = mismatch_s, baseOnly = baseonly)
  mat[!mat%in%c(match_s,mismatch_s)] <- match_s
  reflAlign <- pairwiseAlignment(s_fa, s2, substitutionMatrix = mat,gapOpening = -mismatch_s, gapExtension = 2,type="local")

  vntr_pos_ref <- str_locate_all(gsub("-","",reflAlign@subject),vntr)[[1]]

  gap_l <- rep(0,nchar(reflAlign@subject))
  gap_l[str_locate_all(as.character(reflAlign@subject),"-")[[1]][,1]] <- 1
  gap_l<- cumsum(gap_l)
  #
  a <- vntr_pos_ref[1,1]+gap_l[vntr_pos_ref[1,1]]
  b <- vntr_pos_ref[nrow(vntr_pos_ref),2]+gap_l[vntr_pos_ref[nrow(vntr_pos_ref),2]]
  b <- b-str_count(substr(reflAlign@pattern,a,b),"-")

  nt <- str_locate_all(toupper(ref_fa),gsub("-","",as.character(reflAlign@pattern)))[[1]][1]



  seq_ref_l <- substr(s_fa,reflAlign@pattern@range@start-len_fl,reflAlign@pattern@range@start-1)
  seq_ref_r <-  substr(s_fa,reflAlign@pattern@range@start+b,reflAlign@pattern@range@start+b+(len_fl-1))





  for(i in 1:length(data)){

    raw_fa <- toupper(data[i])


    seq_raw <- DNAString(raw_fa)

    st <- seq(1,nchar(seq_raw),20000-len_fl)
    en <- c(st[-1]-1+len_fl,nchar( seq_raw))

    tmp_mismatch=0
    id = 0
    while(id==0){
      fl_l <- flanking_seq_match(seq=seq_raw,pattern=seq_ref_l,start=st,end=en,max.mismatch=max.mismatch,fixed = T,ext=0)


      if(tmp_mismatch==max.mismatch){
        id = 1
      }else{
        id = nrow(fl_l)
      }
      tmp_mismatch = tmp_mismatch+1
    }

    tmp_mismatch=0
    id = 0
    while(id==0){
      fl_r <- flanking_seq_match(seq=seq_raw,pattern=seq_ref_r,start=st,end=en,max.mismatch=max.mismatch,fixed = T,ext=0)


      if(tmp_mismatch==max.mismatch){
        id = 1
      }else{
        id = nrow(fl_r)
      }
      tmp_mismatch = tmp_mismatch+1
    }

    # fl_l <- flanking_seq_match(seq=seq_raw,pattern=seq_ref_l,start=st,end=en,max.mismatch=max.mismatch,fixed = T,ext=0)
    # fl_r <- flanking_seq_match(seq=seq_raw,pattern=seq_ref_r,start=st,end=en,max.mismatch=max.mismatch,fixed = T,ext=0)

    st_ext <- seq(1,nchar(seq_raw),20000-fseqExtend)
    en_ext <- c(st_ext[-1]-1+fseqExtend,nchar( seq_raw))

    if(nrow(fl_l)!=1){
      seq_ref_l_1 <- substr(s_fa,max(1,reflAlign@pattern@range@start-fseqExtend),reflAlign@pattern@range@start-1)
      tmp_mismatch.ext=0
      id = 0
      while(id==0){
        fl_l <- flanking_seq_match(seq=seq_raw,pattern=seq_ref_l_1,start=st_ext,end=en_ext,max.mismatch=tmp_mismatch.ext,fixed = F,ext=1)
        if(tmp_mismatch.ext==max.mismatch.ext){
          id = 1
        }else{
          id = nrow(fl_l)
        }
        tmp_mismatch.ext = tmp_mismatch.ext+1
      }
      probe_add_l[i] <- 1
    }



    if(nrow(fl_r)!=1){
      seq_ref_r_1 <- substr(s_fa,reflAlign@pattern@range@start+reflAlign@pattern@range@width,min(nchar(s_fa),reflAlign@pattern@range@start+reflAlign@pattern@range@width+(fseqExtend-1)))
      tmp_mismatch.ext=0
      id = 0
      while(id==0){
        fl_r <- flanking_seq_match(seq=seq_raw,pattern=seq_ref_r_1,start=st_ext,end=en_ext,max.mismatch=tmp_mismatch.ext,fixed = F,ext=1)
        if(tmp_mismatch.ext==max.mismatch.ext){
          id = 1
        }else{
          id = nrow(fl_r)
        }
        tmp_mismatch.ext = tmp_mismatch.ext+1
      }
      probe_add_r[i] <- 1
    }


    if(nrow(fl_l)!=1|nrow(fl_r)!=1){
      seq_txt <- ""
      txt <- ""
    }else{
      seq_txt <- substr(raw_fa,fl_l$start+fl_l$width,fl_r$start-1)

      deg <- setdiff(LETTERS,c("A","T","C","G"))

      count_n[i] <- str_count(seq_txt,"N")
      count_deg[i] <- str_count(seq_txt,paste0(setdiff(LETTERS,c("A","T","C","G","N")),collapse = "|"))
      count_len[i] <- nchar(seq_txt)

      seq_txt <- gsub(ifelse(baseonly==T,paste0(deg,collapse = "|"),""),"",seq_txt)

      txt <- seq_txt

    }


    if(nchar(txt)!=0&nchar(seq_txt)>=nchar(vntr)&nchar(seq_txt)<=nchar(vntr)*100){

      if(nchar(vntr)==1){
        vntr_align[[i]] <- txt
        start_pos[i] <- fl_l$start+fl_l$width

        s1_1 <- txt
        s2_1 <- rep(vntr,nchar(txt))
        TF <- table(factor(strsplit(s1_1,"")[[1]]==strsplit(s2_1,"")[[1]], levels = c("TRUE","FALSE")))
        vntr_match_norm[i] <- TF["TRUE"]/sum(TF)
        match[i] <- TF["TRUE"]
        mismatch[i] <- TF["FALSE"]
        indel[i] <- 0
        r[i] <- str_count(txt,vntr)

        r_start[i] <- 1
        r_end[i] <- nchar(txt)

      }else{
        s1 <- DNAString(txt)
        s2 <- DNAString(vntr_t)

        mat <- nucleotideSubstitutionMatrix(match = match_s, mismatch = mismatch_s, baseOnly = F)
        mat[!mat%in%c(match_s,mismatch_s)] <- ifelse(baseonly==T,mismatch_s,match_s)


        globalAlign <- pairwiseAlignment(s1,s2, substitutionMatrix = mat,gapOpening = -mismatch_s, gapExtension = 2,type="overlap")



        s1_1 <- globalAlign@pattern
        s2_1 <- globalAlign@subject

        vntr_pos <- str_locate_all(gsub("-","",s2_1),vntr)[[1]]


        if(nrow(vntr_pos)!=0){
          gap_l <- rep(0,nchar(s2_1))
          gap_l[str_locate_all(as.character(s2_1),"-")[[1]][,1]] <- 1
          gap_l<- cumsum(gap_l)

          a <- vntr_pos[1,1]+gap_l[vntr_pos[1,1]]
          b <- vntr_pos[nrow(vntr_pos),2]+gap_l[vntr_pos[nrow(vntr_pos),2]]

          ATCG_p[i] <- paste(sapply(1:nrow(vntr_pos), function(x) str_count(substr(s1_1,vntr_pos[x,1]+gap_l[vntr_pos[x,1]],vntr_pos[x,2]+gap_l[vntr_pos[x,2]]),"A|T|C|G|-"))-sapply(1:nrow(vntr_pos), function(x) str_count(substr(s2_1,vntr_pos[x,1]+gap_l[vntr_pos[x,1]],vntr_pos[x,2]+gap_l[vntr_pos[x,2]]),"-")),collapse = "-")

          s1_1 <- substr(s1_1,a,b)
          s2_1 <- substr(s2_1,a,b)



          vntr_align[[i]] <- txt
          start_pos[i] <- fl_l$start+fl_l$width
          TF <- table(factor(strsplit(s1_1,"")[[1]]==strsplit(s2_1,"")[[1]], levels = c("TRUE","FALSE")))
          vntr_match_norm[i] <- TF["TRUE"]/sum(TF)
          vntr_loc <- str_locate_all(gsub("-","",as.character(globalAlign@subject)),vntr)[[1]]
          r_start[i] <- str_locate(txt,gsub("-","",s1_1))[1]
          r_end[i] <- str_locate(txt,gsub("-","",s1_1))[2]

          match[i] <- TF["TRUE"]
          mismatch[i] <- TF["FALSE"]
          indel[i] <- str_count(s2_1,"-")+str_count(s1_1,"-")
          r[i]<-nrow(vntr_loc)



        }else{

          r[i]<-0
          vntr_align[[i]] <- ""
        }

      }



    }else{
      r[i]<- ifelse(nrow(fl_l)==1&&nrow(fl_r)==1,
                    ifelse(nchar(seq_txt)>nchar(vntr)*100,NA,0),NA)
      vntr_align[[i]] <- ""
    }


  }


  names(vntr_align) <- names(data)

  #===============


  out <- data.frame(ID=names(data),r=r,match=match,mismatch=mismatch,
                    indel=indel,score=match*match_s+mismatch*mismatch_s,
                    start_pos=start_pos,
                    probe_add_l=probe_add_l,probe_add_r=probe_add_r,
                    count_n=count_n,count_deg=count_deg,count_len=count_len,
                    vntr_match_norm=vntr_match_norm)
  out$start_pos[out$r==0] <- NA
  if(length(out$ID[which(out$score<0)])>0){
    vntr_align[out$ID[which(out$score<0)]]<- ""
  }
  out[which(out$score<0|is.na(out$r)),-1] <- NA

  if(VNTRoutput==T){
    dir.create("output",showWarnings = F)
    dir.create("output/VNTRCaller",showWarnings = F)
    write.csv(out,paste0("output/VNTRCaller/VNTR_listbaseonly_",baseonly,"_flanking",len_fl,"_nt",nt,"_",vntr,".csv"),row.names = F)

    seqinr::write.fasta(vntr_align,names(vntr_align), paste0("output/VNTRCaller/VNTR_listbaseonly_",baseonly,"_flanking",len_fl,"_nt",nt,"_",vntr,".fas"))
  }

  list(paste0("VNTR_nt",paste(nt,vntr,"baseonly",baseonly,sep = "_")),out,vntr_align)


}


























































#' Calling Variable Number Tandem Repeats (VNTR) for the genome sequence of
#' monkeypox virus (MPXV) %% ~~function to do ... ~~
#' 
#' The function \code{VNTRcaller} estimates the copy numbers of VNTRs. %% ~~ A
#' concise (1-5 lines) description of what the function does. ~~
#' 
#' 
#' @param data MPXV sequences from a file in FASTA format
#' @param vntr sequence of the repetitive unit
#' @param match_s matching weight
#' @param mismatch_s mismatching penalty
#' @param regionStart starting position of a VNTR region
#' @param regionEnd ending position of a VNTR region
#' @param baseonly logical. If TRUE, only the letters of the nucleotide bases
#' (i.e. A,C,G,T) are considered. If FALSE, the degenerate codes are also
#' considered.
#' @param fseqExtend the length of extended flanking sequence
#' @param VNTRoutput logical. If TRUE, the output is written to .csv files.
#' @param tracker logical. If TRUE, call function \code{\link{VNTRtracker}}.
#' @return \item{ID}{name of MPXV sequences} \item{r}{copy of tandem repeats}
#' \item{match}{number of matches} \item{mismatch}{number of mismatches}
#' \item{indel}{number of insertions and deletions (indels)}
#' \item{score}{alignment score of a VNTR region for a query strain}
#' \item{start_pos}{starting position of the VNTR region for a query strain}
#' \item{count_n}{number of N characters} \item{count_deg}{number of non-ATCGN
#' characters} \item{count_len}{VNTR region length}
#' \item{vntr_match_norm}{proportion of number of N characters to VNTR region
#' length} %% ...
#' @references %% ~put references to the literature/web site here ~ Yang, H.-C.
#' et al (2022) Monkeypox genome contains variable number tandem repeats
#' enabling accurate virus tracking.
#' 
#' Pages H, Aboyoun P, Gentleman R, DebRoy S (2022) \emph{Biostrings: Efficient
#' manipulation of biological strings}.  R package version 2.64.0.
#' %%\href{https://bioconductor.org/packages/Biostringshttps://bioconductor.org/packages/Biostrings}.
#' @examples
#' 
#' ## Read MPXV example sequences from a file in FASTA format.
#' MPXVseq <- read.fasta(system.file("extdata/MPXV_seq_example.fasta.gz", package = "MPXV"), as.string = T)
#' 
#' ## VNTR
#' vntr <- c("T","TATGATGGA","AT","ATATACATT")
#' regionStart <- c(132436,150542,173240,178413)
#' regionEnd <- c(133216,151501,173320,179244)
#' 
#' ## parameter settings
#' match_s = 2
#' mismatch_s = -5
#' baseonly = TRUE
#' fseqExtend = 150
#' VNTRoutput = FALSE
#' tracker = FALSE
#' 
#' ## computes the copy of the variable number tandem repeats
#' out <- VNTRcaller(data = MPXVseq, vntr = vntr,
#'                   regionStart = regionStart, regionEnd = regionEnd,
#'                   match_s = match_s, mismatch_s=mismatch_s,
#'                   baseonly = baseonly, fseqExtend = fseqExtend,
#'                   VNTRoutput = VNTRoutput,
#'                   tracker = tracker)
#' 
#' ## For more details about the output of VNTRcaller, please vitsit https://github.com/teresayang/MPXV_VNTR.
#' 
#' 
#' @export VNTRcaller
VNTRcaller <- function(data, vntr=vntr, match_s=match_s, mismatch_s=mismatch_s,
                       regionStart=regionStart, regionEnd=regionEnd,baseonly = T,
                       fseqExtend = fseqExtend,VNTRoutput=F,tracker=F,
                       len_fl = 25, mis_prop = 0.1){
  if(sum(is.na(as.numeric(c(regionStart,regionEnd))))!=0){
    stop("regionStart or regionEnd should be numeric.")
  }
  if(!all.equal(length(vntr),length(regionStart),length(regionEnd ))){
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
  invisible(sapply(1:length(vntr), function(x) out[[x]] <<-VNTR_sub(data, vntr=vntr[x],regionStart=regionStart[x], regionEnd=regionEnd[x], match_s=match_s, mismatch_s=mismatch_s,baseonly = baseonly,fseqExtend=fseqExtend,VNTRoutput=VNTRoutput,len_fl = len_fl, mis_prop = mis_prop)))
  names(out) <- sapply(1:length(out), function(x)out[[x]][[1]])


  if(length(data)==1){
    dt <- c(names(data), sapply(1:length(out), function(x) out[[x]][[2]][,2]))
    names(dt) <- c("ID",sapply(1:length(names(out)), function(x)paste(strsplit(names(out)[x],"_")[[1]][2:3],collapse  = "_") ))
    dt <- data.frame(t(dt))
  }else{
    dt <- data.frame(ID=names(data), sapply(1:length(out), function(x) out[[x]][[2]][,2]))
    colnames(dt)[-1] <- sapply(1:length(names(out)), function(x)paste(strsplit(names(out)[x],"_")[[1]][2:3],collapse  = "_") )
  }
  if(VNTRoutput==T){
    write.csv(dt,paste0("output/VNTRCaller/VNTR_list",paste("baseonly",baseonly,"flanking",len_fl,sep = "_"),".csv"),row.names = F)
    print(paste0("Calling VNTR is done and the results are saved in ", getwd(), "/output/VNTRCaller"))
  }
  if(tracker==T){

    if(all(colnames(dt)%in%c("ID","nt133095_T", "nt150554_TATGATGGA", "nt173267_AT", "nt179074_ATATACATT"))){
      dt <- dt[,c("ID","nt133095_T", "nt150554_TATGATGGA", "nt173267_AT", "nt179074_ATATACATT")]
    }else{
      db <- c("ID","nt133095_T", "nt150554_TATGATGGA", "nt173267_AT", "nt179074_ATATACATT")[!c("ID","nt133095_T", "nt150554_TATGATGGA", "nt173267_AT", "nt179074_ATATACATT")%in%c(colnames(dt))]
      m <- matrix(NA, ncol=length(db),nrow=nrow(dt))
      colnames(m) <- db
      dt <- cbind(dt,m)
      dt <- dt[,c("ID","nt133095_T", "nt150554_TATGATGGA", "nt173267_AT", "nt179074_ATATACATT")]
    }
    tryCatch(
      {
        invisible(sapply(1:nrow(dt), function(x)VNTRtracker(r = as.numeric(unlist(dt[x,2:5])),file_name= dt$ID[x])))
      },
      error=function(error_message) {
        message(error_message)

      }
    )
  }
  out

}
