setwd("C:/Users/user/Desktop/new")
library(Biostrings)
# ======== load hits ========
hits <- read.table('K_hits.recovered.last.bed',stringsAsFactors=FALSE, sep='\t',header=FALSE);
hits$flag <- 1;
inds <- which(hits$V5 < 0.25);
hits$flag[inds] <- 0;

for (i in 1:length(inds))
{
  if (length(which(hits$V4 %in% hits$V4[inds[i]] & abs(hits$V2 - hits$V2[inds[i]]) < 30)) > 1)
    hits$flag[inds[i]] <- 1
}
myhits <- hits$V8[which(hits$flag == 1)]

tmp <- read.table('K_hits_for_conservation.tab',stringsAsFactors=FALSE, sep='\t',header=TRUE)
#myhits <- hits$V8[which(hits$V5 >= 0.25)]
tmp2 <- tmp[which(tmp$genomic_location %in% myhits),]
tmp2$bed <- paste(tmp2$genomic_location,tmp2$hgnc,tmp2$refseq,tmp2$position,sep='__')
seqs <- list(tmp2$seq_50bp)
mins <- rep(-1,length(seqs[[1]]))
maxs <- rep(-1,length(seqs[[1]]))
tmp2$len <- rep(-1,length(tmp2$seq_50bp))

# ======== initialise file counts ========
counts_7 = rep(-1,50);counts_r1_7 = rep(-1,50);counts_r2_7 = rep(-1,50);counts_r3_7 = rep(-1,50)
counts_long = rep(-1,50);counts_r1_long = rep(-1,50);counts_r2_long = rep(-1,50);counts_r3_long = rep(-1,50)
counts_bulges = rep(-1,50);counts_r1_bulges = rep(-1,50);counts_r2_bulges = rep(-1,50);counts_r3_bulges = rep(-1,50);counts_GG = rep(-1,50);counts_40 = rep(-1,50);counts_other = rep(-1,50);
counts_pqs = rep(-1,50);counts_r1_pqs = rep(-1,50);counts_r2_pqs = rep(-1,50);counts_r3_pqs = rep(-1,50)
counts_fraction = rep(-1,50);counts_r1_fraction = rep(-1,50);counts_r2_fraction = rep(-1,50);counts_r3_fraction = rep(-1,50)
counts_5 = rep(-1,50);counts_r1_5 = rep(-1,50);counts_r2_5 = rep(-1,50);counts_r3_5 = rep(-1,50)
counts_3 = rep(-1,50);counts_r1_3 = rep(-1,50);counts_r2_3 = rep(-1,50);counts_r3_3 = rep(-1,50)
counts_A = vector("list", 50);counts_r1_A = vector("list", 50);counts_r2_A = vector("list", 50);counts_r3_A = vector("list", 50)
counts_T = vector("list", 50);counts_r1_T = vector("list", 50);counts_r2_T = vector("list", 50);counts_r3_T = vector("list", 50)
counts_C = vector("list", 50);counts_r1_C = vector("list", 50);counts_r2_C = vector("list", 50);counts_r3_C = vector("list", 50)
counts_G = vector("list", 50);counts_r1_G = vector("list", 50);counts_r2_G = vector("list", 50);counts_r3_G = vector("list", 50)
counts_AT = vector("list", 50);counts_r1_AT = vector("list", 50);counts_r2_AT = vector("list", 50);counts_r3_AT = vector("list", 50)
counts_CG = vector("list", 50);counts_r1_CG = vector("list", 50);counts_r2_CG = vector("list", 50);counts_r3_CG = vector("list", 50)
counts_skewGC = vector("list", 50);counts_r1_skewGC = vector("list", 50);counts_r2_skewGC = vector("list", 50);counts_r3_skewGC = vector("list", 50)
counts_skewAT = vector("list", 50);counts_r1_skewAT = vector("list", 50);counts_r2_skewAT = vector("list", 50);counts_r3_skewAT = vector("list", 50)

count_f = 1

for (s1 in 1:length(seqs)) {  
  
  seq <- toupper(seqs[[s1]])
  
  # ======== initialise motif counts ========
  
  countsA = rep(0, length(seq));countsT = rep(0, length(seq));countsC = rep(0, length(seq));countsG = rep(0, length(seq))
  countsAT = rep(0, length(seq));countsCG = rep(0, length(seq));skewGC = rep(0, length(seq));skewAT = rep(0, length(seq))
  
  count_A7 = 0;count_C7 = 0;count_G7 = 0;count_T7 = 0;count7 = 0;count5 = 0;count3 = 0;countC7 = 0;count12 = 0;count21 = 0;count_GG = 0;count_GGATG_1 = 0;count_GGATG_2 = 0;count_GGATG_3 = 0;count_GGATG_4 = 0
  count_GGATG_1_2 = 0;count_GGATG_1_3 = 0;count_GGATG_1_4 = 0;count_GGATG_2_3 = 0;count_GGATG_2_4 = 0;count_GGATG_3_4 = 0
  count_GGATG_1_2_3 = 0;count_GGATG_1_2_4 = 0;count_GGATG_2_3_4 = 0;count_GGATG_1_3_4 = 0;
  count_GGATG_1_1 = 0;count_GGATG_2_2 = 0;count_GGATG_3_3 = 0;count_GGATG_4_4 = 0
  count_GGATG_1n_1n = 0;count_GGATG_2n_2n = 0;count_GGATG_3n_3n = 0;count_GGATG_4n_4n = 0;count_GGATG_1n_2n = 0;count_GGATG_1n_3n = 0;count_GGATG_1n_4n = 0;count_GGATG_2n_3n = 0;count_GGATG_2n_4n = 0;count_GGATG_3n_4n = 0;count_GGATG_1n_2n_3n = 0;count_GGATG_1n_3n_4n = 0;count_GGATG_1n_2n_4n = 0;count_GGATG_2n_3n_4n = 0;
  
  inds_A7 = rep(1,0);inds_C7 = rep(1,0);inds_G7 = rep(1,0);inds_T7 = rep(1,0);inds_12 = rep(1,0);inds_7 = rep(1,0);inds_5 = rep(1,0);inds_3 = rep(1,0);inds_C7 = rep(1,0);inds_21 = rep(1,0);inds_GG = rep(1,0);inds_GGATG_1 = rep(1,0);inds_GGATG_2 = rep(1,0);inds_GGATG_3 = rep(1,0);inds_GGATG_4 = rep(1,0)
  inds_GGATG_1_2 = rep(1,0);inds_GGATG_1_3 = rep(1,0);inds_GGATG_1_4 = rep(1,0);inds_GGATG_2_3 = rep(1,0);inds_GGATG_2_4 = rep(1,0);inds_GGATG_3_4 = rep(1,0)
  inds_GGATG_1_2_3 = rep(1,0);inds_GGATG_1_2_4 = rep(1,0);inds_GGATG_2_3_4 = rep(1,0);inds_GGATG_1_3_4 = rep(1,0)
  inds_GGATG_1_1 = rep(1,0);inds_GGATG_2_2 = rep(1,0); inds_GGATG_3_3 = rep(1,0); inds_GGATG_4_4 = rep(1,0);
  inds_GGATG_1n_1n = rep(1,0);inds_GGATG_2n_2n = rep(1,0);inds_GGATG_3n_3n = rep(1,0);inds_GGATG_4n_4n = rep(1,0); inds_GGATG_1n_2n = rep(1,0);inds_GGATG_1n_3n = rep(1,0);inds_GGATG_1n_4n = rep(1,0);inds_GGATG_2n_3n = rep(1,0);inds_GGATG_2n_4n = rep(1,0);inds_GGATG_3n_4n = rep(1,0);inds_GGATG_1n_2n_3n = rep(1,0);inds_GGATG_1n_2n_4n = rep(1,0);inds_GGATG_1n_3n_4n = rep(1,0);inds_GGATG_2n_3n_4n = rep(1,0)
  
  
  for (i in 1:length(seq))
  {
    min_max_len <- 1000
    flag <- TRUE
    aseq <- seq[i]
    
    
    tmp <- strsplit(aseq,'')[[1]]
    aseq_len <- length(tmp)
    tmp2$len[i] <- length(tmp)
    len <- length(tmp)/100
    countsT[i] <- length(which(tmp == 'T'))/len
    countsA[i] <- length(which(tmp == 'A'))/len
    countsC[i] <- length(which(tmp == 'C'))/len
    countsG[i] <- length(which(tmp == 'G'))/len
    countsCG[i] <- (countsC[i] + countsG[i])/len
    countsAT[i] <- (countsA[i] + countsT[i])/len
    skewGC[i] <- (countsG[i] - countsC[i])/(countsG[i] + countsC[i])
    skewAT[i] <- (countsA[i] - countsT[i])/(countsA[i] + countsT[i])
    
    if(i%%1000 == 0)
      print(i)
    if (length(grep('A{7,}?',aseq))>0) 
    {count_A7 = count_A7 + 1 ; inds_A7= c(inds_A7,i);}
    if (length(grep('T{7,}?',aseq))>0) 
    {count_T7 = count_T7 + 1 ; inds_T7= c(inds_T7,i);}
    if (length(grep('C{7,}?',aseq))>0) 
    {count_C7 = count_C7 + 1 ; inds_C7= c(inds_C7,i);}
    if (length(grep('G{7,}?',aseq))>0) 
    {count_G7 = count_G7 + 1 ; inds_G7= c(inds_G7,i);}
    if (length(grep('(G{3,}?[ATGC]{1,3}?){3,}?G{3,}?',aseq))>0) 
    {count3 = count3 + 1 ; inds_3= c(inds_3,i); flag<-FALSE}
    if (length(grep('(G{3,}?[ATGC]{1,5}?){3,}?G{3,}?',aseq))>0) 
    {count5 = count5 + 1 ; inds_5= c(inds_5,i); flag<-FALSE}
    if (length(grep('(G{3,}?[ATGC]{1,7}?){3,}?G{3,}?',aseq))>0) 
    {
      count7 = count7 + 1 ; inds_7= c(inds_7,i); flag<-FALSE
      aseq1 <-as.character(reverse(DNAString(aseq)))
      minmax <- gregexpr('(G{3,}?[ATGC]{1,7}?){3}?G{3,}?',aseq1)[[1]]
      ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; 
      amin1 <- minmax[ind_short]
      amax1 <- attr(minmax,"match.length")[ind_short]+amin1-1
      amax <- aseq_len - amin1 + 1
      amin <- aseq_len - amax1 + 1
      # if (amax - amin < min_max_len)
      if (amax > maxs[i])
      {
        mins[i] <- amin
        maxs[i] <- amax
        min_max_len <- amax - amin
      }
    }
    if (length(grep('(C{3,}?[ATGC]{1,7}?){3,}?C{3,}?',aseq))>0) 
    {countC7 = countC7 + 1 ; inds_C7= c(inds_C7,i); }
    if (length(grep('(G{3,}?[ATGC]{1,12}){3,}?G{3,}?',aseq))>0) 
    {
      count12 = count12 + 1 ; inds_12= c(inds_12,i); flag<-FALSE
      aseq1 <-as.character(reverse(DNAString(aseq)))
      minmax <- gregexpr('(G{3,}?[ATGC]{1,12}?){3}?G{3,}?',aseq1)[[1]]
      ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
      amin1 <- minmax[ind_short]
      amax1 <- attr(minmax,"match.length")[ind_short]+amin1-1
      amax <- aseq_len - amin1 + 1
      amin <- aseq_len - amax1 + 1
      if (amax > maxs[i])
      {
        count12 = count12 + 1 ; inds_12= c(inds_12,i); flag<-FALSE
        mins[i] <- amin
        maxs[i] <- amax 
        min_max_len <- amax - amin
      }
    }
    if (length(grep('G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{13,21}G{3,}?[ATGC]{1,7}?G{3,}?',aseq))>0) 
    {
      #count21 = count21 + 1 ; inds_21= c(inds_21,i); flag<-FALSE
      minmax <- gregexpr('G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{13,21}G{3,}?[ATGC]{1,7}?G{3,}?',aseq)[[1]]
      ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1];
      amin <- minmax[ind_short]
      amax <- attr(minmax,"match.length")[ind_short]+amin-1
      if (amax > maxs[i])
      {
        count21 = count21 + 1 ; inds_21= c(inds_21,i); flag<-FALSE
        mins[i] <- amin
        maxs[i] <- amax 
        min_max_len <- amax - amin
      }
    }
    if (flag)
    {
      if (length(grep('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]{1,9}?G|G[AT]{1,9}?GG)',aseq))>0){
        count_GGATG_4 = count_GGATG_4 + 1    ; inds_GGATG_1 = c(inds_GGATG_1,i)
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]{1,9}?G|G[AT]{1,9}?GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]{1,9}?G|G[AT]{1,9}?GG)[ATGC]{1,9}?G{3,}?',aseq))>0){
        count_GGATG_3 = count_GGATG_3 + 1    ; inds_GGATG_2 = c(inds_GGATG_2,i); 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]{1,9}?G|G[AT]{1,9}?GG)[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1];
        amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?(GG[AT]{1,9}?G|G[AT]{1,9}?GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq))>0){
        count_GGATG_2 = count_GGATG_2 + 1    ; inds_GGATG_3 = c(inds_GGATG_3,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?(GG[AT]{1,9}?G|G[AT]{1,9}?GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1];
        amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                    
      }
      if (length(grep('(GG[AT]{1,9}?G|G[AT]{1,9}?GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq))>0) {
        count_GGATG_1 = count_GGATG_1 + 1; inds_GGATG_4 = c(inds_GGATG_4,i)
        minmax <- gregexpr('(GG[AT]{1,9}?G|G[AT]{1,9}?GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      
      if (length(grep('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[C]{1,1}?G|G[C]{1,1}?GG)',aseq))>0){
        count_GGATG_4 = count_GGATG_4 + 1    ; inds_GGATG_1 = c(inds_GGATG_1,i)
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[C]{1,1}?G|G[C]{1,1}?GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[C]{1,1}?G|G[C]{1,1}?GG)[ATGC]{1,9}?G{3,}?',aseq))>0){
        count_GGATG_3 = count_GGATG_3 + 1    ; inds_GGATG_2 = c(inds_GGATG_2,i); 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[CC]{1,1}?G|G[C]{1,1}?GG)[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1];
        amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?(GG[C]{1,1}?G|G[C]{1,1}?GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq))>0){
        count_GGATG_2 = count_GGATG_2 + 1    ; inds_GGATG_3 = c(inds_GGATG_3,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?(GG[C]{1,1}?G|G[C]{1,1}?GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1];
        amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                    
      }
      if (length(grep('(GG[C]{1,1}?G|G[C]{1,1}?GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq))>0) {
        count_GGATG_1 = count_GGATG_1 + 1; inds_GGATG_4 = c(inds_GGATG_4,i)
        minmax <- gregexpr('(GG[C]{1,1}?G|G[C]{1,1}?GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      
      
      if (length(grep('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq))>0) {
        count_GGATG_1_2 = count_GGATG_1_2 + 1; inds_GGATG_1_2 = c(inds_GGATG_1_2,i) 
        minmax <- gregexpr('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?',aseq))>0) {
        count_GGATG_1_3 = count_GGATG_1_3 + 1; inds_GGATG_1_3= c(inds_GGATG_1_3,i) 
        minmax <- gregexpr('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq))>0) {
        count_GGATG_1_4 = count_GGATG_1_4 + 1; inds_GGATG_1_4= c(inds_GGATG_1_4,i) 
        minmax <- gregexpr('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?',aseq))>0){
        count_GGATG_2_3 = count_GGATG_2_3 + 1; inds_GGATG_2_3 = c(inds_GGATG_2_3 ,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq))>0){
        count_GGATG_2_4 = count_GGATG_2_4 + 1; inds_GGATG_2_4 = c(inds_GGATG_2_4 ,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq))>0) {
        count_GGATG_3_4 = count_GGATG_3_4 + 1; inds_GGATG_3_4 = c(inds_GGATG_3_4 ,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?',aseq))>0) {
        count_GGATG_1_2_3 = count_GGATG_1_2_3 + 1; inds_GGATG_1_2_3 = c(inds_GGATG_1_2_3 ,i) 
        minmax <- gregexpr('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq))>0) {
        count_GGATG_1_2_4 = count_GGATG_1_2_4 + 1; inds_GGATG_1_2_4 = c(inds_GGATG_1_2_4 ,i) 
        minmax <- gregexpr('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq))>0) {
        count_GGATG_1_3_4 = count_GGATG_1_3_4 + 1; inds_GGATG_1_3_4 = c(inds_GGATG_1_3_4 ,i) 
        minmax <- gregexpr('(GG[AT]G|G[AT]GG)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq))>0) {
        count_GGATG_2_3_4 = count_GGATG_2_3_4 + 1; inds_GGATG_2_3_4 = c(inds_GGATG_2_3_4 ,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)[ATGC]{1,9}?(GG[AT]G|G[AT]GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(G[ATC]G[ATC]G)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq))>0) {
        count_GGATG_1_1 = count_GGATG_1_1 + 1; inds_GGATG_1_1 = c(inds_GGATG_1_1 ,i) 
        minmax <- gregexpr('(G[AT]G[AT]G)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?(G[ATC]G[ATC]G)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq))>0) {
        count_GGATG_2_2 = count_GGATG_2_2 + 1; inds_GGATG_2_2 = c(inds_GGATG_2_2,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?(G[AT]G[AT]G)[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(G[ATC]G[ATC]G)[ATGC]{1,9}?G{3,}?',aseq))>0) {
        count_GGATG_3_3 = count_GGATG_3_3 + 1; inds_GGATG_3_3 = c(inds_GGATG_3_3,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(G[AT]G[AT]G)[ATGC]{1,9}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(G[ATC]G[ATC]G)',aseq))>0) {
        count_GGATG_4_4 = count_GGATG_4_4 + 1; inds_GGATG_4_4 = c(inds_GGATG_4_4,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?G{3,}?[ATGC]{1,9}?(G[AT]G[AT]G)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(G[AT]{2,5}?G[AT]{2,5}?G)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?',aseq))>0) {
        count_GGATG_1n_1n = count_GGATG_1n_1n + 1; inds_GGATG_1n_1n = c(inds_GGATG_1n_1n,i) 
        minmax <- gregexpr('(G[AT]{2,5}?G[AT]{2,5}?G)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,7}?(G[AT]{2,5}?G[AT]{2,5}?G)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?',aseq))>0) {
        count_GGATG_2n_2n = count_GGATG_2n_2n + 1; inds_GGATG_2n_2n = c(inds_GGATG_2n_2n,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,7}?(G[AT]{2,5}?G[AT]{2,5}?G)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(G[AT]{2,5}?G[AT]{2,5}?G)[ATGC]{1,7}?G{3,}?',aseq))>0) {
        count_GGATG_3n_3n = count_GGATG_3n_3n + 1; inds_GGATG_3n_3n = c(inds_GGATG_3n_3n,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(G[AT]{2,5}?G[AT]{2,5}?G)[ATGC]{1,7}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(G[AT]{2,5}?G[AT]{2,5}?G)',aseq))>0) {
        count_GGATG_4n_4n = count_GGATG_4n_4n + 1; inds_GGATG_4n_4n = c(inds_GGATG_4n_4n,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(G[AT]{2,5}?G[AT]{2,5}?G)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?',aseq))>0) {
        count_GGATG_1n_2n = count_GGATG_1n_2n + 1; inds_GGATG_1n_2n = c(inds_GGATG_1n_2n,i) 
        minmax <- gregexpr('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?',aseq))>0) {
        count_GGATG_2n_3n = count_GGATG_2n_3n + 1; inds_GGATG_2n_3n = c(inds_GGATG_2n_3n,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?',aseq))>0) {
        count_GGATG_1n_3n = count_GGATG_1n_3n + 1; inds_GGATG_1n_3n = c(inds_GGATG_1n_3n,i) 
        minmax <- gregexpr('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq))>0) {
        count_GGATG_1n_4n = count_GGATG_1n_4n + 1; inds_GGATG_1n_4n = c(inds_GGATG_1n_4n,i) 
        minmax <- gregexpr('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq))>0) {
        count_GGATG_2n_4n = count_GGATG_2n_4n + 1; inds_GGATG_2n_4n = c(inds_GGATG_2n_4n,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq))>0) {
        count_GGATG_3n_4n = count_GGATG_3n_4n + 1; inds_GGATG_3n_4n = c(inds_GGATG_3n_4n,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq))>0) {
        count_GGATG_2n_4n = count_GGATG_2n_4n + 1; inds_GGATG_2n_4n = c(inds_GGATG_2n_4n,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?',aseq))>0) {
        count_GGATG_1n_2n_3n = count_GGATG_1n_2n_3n + 1; inds_GGATG_1n_2n_3n = c(inds_GGATG_1n_2n_3n,i) 
        minmax <- gregexpr('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq))>0) {
        count_GGATG_1n_2n_4n = count_GGATG_1n_2n_4n + 1; inds_GGATG_1n_2n_4n = c(inds_GGATG_1n_2n_4n,i) 
        minmax <- gregexpr('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq))>0) {
        count_GGATG_1n_3n_4n = count_GGATG_1n_3n_4n + 1; inds_GGATG_1n_3n_4n = c(inds_GGATG_1n_3n_4n,i) 
        minmax <- gregexpr('(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if (length(grep('G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq))>0) {
        count_GGATG_2n_3n_4n = count_GGATG_2n_3n_4n + 1; inds_GGATG_2n_3n_4n = c(inds_GGATG_2n_3n_4n,i) 
        minmax <- gregexpr('G{3,}?[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)[ATGC]{1,7}?(GG[AT]{2,5}?G|G[AT]{2,5}?GG)',aseq)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1]; amin <- minmax[ind_short]
        amax <- attr(minmax,"match.length")[ind_short]+amin-1
        if (amax > maxs[i])
        {
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }                  
      }
      if(length(grep('(G{2,}?[ATGC]{1,9}?){3,}?G{2,}?',aseq))>0) {  #&& mins[i] == -1) {
        #count_GG= count_GG+ 1; inds_GG = c(inds_GG,i)
        
        aseq1 <-as.character(reverse(DNAString(aseq)))
        minmax <- gregexpr('(G{2,}?[ATGC]{1,9}?){3}?G{2,}?',aseq1)[[1]]
        ind_short <- which(min(attr(minmax,"match.length")) == attr(minmax,"match.length"))[1];
        amin <- minmax[ind_short]
        amin1 <- minmax[ind_short]
        amax1 <- attr(minmax,"match.length")[ind_short]+amin1-1
        amax <- aseq_len - amin1 + 1
        amin <- aseq_len - amax1 + 1
        if (amax > maxs[i])
        {
          count_GG= count_GG+ 1; inds_GG = c(inds_GG,i)
          mins[i] <- amin
          maxs[i] <- amax 
          min_max_len <- amax - amin
        }
      }
      
    }     # end if flag for the bulges
    len_end <- aseq_len - maxs[i]
    if (len_end <= 3)
    {
      if (length(which(tmp[(maxs[i]+1):aseq_len] == rep('G',len_end))) == len_end)
        maxs[i] <- aseq_len
    }
  } # end for i  
  
  inds_long <- unique(setdiff(c(inds_12,inds_21, inds_7, inds_5, inds_3), c(inds_7, inds_5, inds_3)))
  inds_bulges <- unique(c(inds_GGATG_1, inds_GGATG_2, inds_GGATG_3, inds_GGATG_4, inds_GGATG_1_1, inds_GGATG_2_2, inds_GGATG_3_3, inds_GGATG_4_4, inds_GGATG_1_2, inds_GGATG_1_3, inds_GGATG_1_4, inds_GGATG_2_3, inds_GGATG_2_4, inds_GGATG_3_4, inds_GGATG_1_2_3, inds_GGATG_1_2_4, inds_GGATG_1_3_4, inds_GGATG_2_3_4, inds_GGATG_1n_1n, inds_GGATG_2n_2n, inds_GGATG_3n_3n, inds_GGATG_4n_4n, inds_GGATG_1n_2n, inds_GGATG_1n_3n, inds_GGATG_1n_4n, inds_GGATG_2n_3n, inds_GGATG_2n_4n, inds_GGATG_3n_4n, inds_GGATG_1n_2n_4n, inds_GGATG_1n_2n_3n, inds_GGATG_1n_3n_4n, inds_GGATG_2n_3n_4n))
  #inds_GG_pure <- setdiff(inds_GG,c(inds_long,inds_bulges, inds_7))
  inds_GG_pure <- setdiff(inds_GG,c(inds_long,inds_7))
  inds_bulges <- setdiff(inds_bulges, inds_GG_pure)
  inds_40 <- which(countsG >= 40)
  inds_all <- unique(c(inds_long,inds_bulges, inds_GG, inds_7, inds_40))
  inds_40_pure <- setdiff(inds_40,c(inds_GG,inds_long,inds_bulges, inds_7))
  inds_all_not <- setdiff(1:length(seq), inds_all)  
  
  counts_7[count_f] = length(inds_7)-length(inds_5);
  counts_5[count_f] = length(inds_5)-length(inds_3);
  counts_3[count_f] = length(inds_3)
  counts_long[count_f] = length(inds_long)
  counts_bulges[count_f] = length(inds_bulges)
  counts_GG[count_f] = length(inds_GG_pure) 
  counts_40[count_f] = length(inds_40_pure)
  counts_other[count_f] = length(seq) - length(inds_all) 
  
} # end for s1


tmp2$class <- rep('NA',length(tmp2$position))
tmp2$class[inds_7] <- 'PQ'
tmp2$class[inds_long] <- 'long'
tmp2$class[inds_bulges] <- 'bulge'
tmp2$class[inds_GG_pure] <- '2-tetrad'
tmp2$class[inds_40_pure] <- 'G_40_pct'
tmp2$start_constraint <- mins
tmp2$end_constraint <- maxs 

tmp2_notoher <- tmp2[-inds_all_not,]         
tmp2_notoher$start_constraint <- tmp2_notoher$position-tmp2_notoher$len+tmp2_notoher$start_constraint
tmp2_notoher$end_constraint <- tmp2_notoher$position-tmp2_notoher$len+tmp2_notoher$end_constraint

write.table(tmp2,'K_hits_classes_table.max_max.tab',col.names = colnames(tmp2),row.names = FALSE,quote=FALSE,sep='\t')
write.table(tmp2_notoher,'K_hits_noothers_classes_table.max_max.tab',col.names = colnames(tmp2_notoher),row.names = FALSE,quote=FALSE,sep='\t')


write.table(seq[inds_3],'KPDS_hits.025.inds_3.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_5],'KPDS_hits.025.inds_5.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_7],'KPDS_hits.025.inds_7.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_long],'KPDS_hits.025.inds_long.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_bulges],'KPDS_hits.025.inds_bulges.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_GG_pure],'KPDS_hits.025.inds_GG.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_40_pure],'KPDS_hits.025.inds_40.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_all_not],'KPDS_hits.025.inds_other.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)

write.table(tmp2$bed[inds_3],'KPDS_hits.025.inds_3.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_5],'KPDS_hits.025.inds_5.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_7],'KPDS_hits.025.inds_7.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_long],'KPDS_hits.025.inds_long.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_bulges],'KPDS_hits.025.inds_bulges.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_GG_pure],'KPDS_hits.025.inds_GG.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_40_pure],'KPDS_hits.025.inds_40.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_all_not],'KPDS_hits.025.inds_other.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)

write.table(seq[inds_3],'K_hits.all.inds_3.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_5],'K_hits.all.inds_5.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_7],'K_hits.all.inds_7.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_long],'K_hits.all.inds_long.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_bulges],'K_hits.all.inds_bulges.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_GG_pure],'K_hits.all.inds_GG.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_40_pure],'K_hits.all.inds_40.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(seq[inds_all_not],'K_hits.all.inds_other.seq',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)

write.table(tmp2$bed[inds_3],'K_hits.all.inds_3.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_5],'K_hits.all.inds_5.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_7],'K_hits.all.inds_7.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_long],'K_hits.all.inds_long.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_bulges],'K_hits.all.inds_bulges.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_GG_pure],'K_hits.all.inds_GG.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_40_pure],'K_hits.all.inds_40.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)
write.table(tmp2$bed[inds_all_not],'K_hits.all.inds_other.ind',col.names=FALSE,row.names=FALSE, sep='\t',quote=FALSE)