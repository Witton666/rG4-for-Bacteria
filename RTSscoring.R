setwd('rG4/')

cond1 <- 'KPDS'
# cond1 <- 'K' 
cond2 <- 'Li'

  fname<-'chr1'

  chr1 <- read.table('chr1.txt', stringsAsFactors=FALSE, sep=" ", header=FALSE)
  
  achr <- unique(chr1$V1)
  astrand <- unique(chr1$V6)
  
  
  # ====== filter parameters to remove exon borders and detect step ======
  filter_len      <- 10
  filter_n_border <- 1
  filter_k_cov    <- 6
  filter_th_step  <- 0.8/4


 
  adata <- chr1[,c(11,12,13)]#K
  adata2 <- chr1[,c(8,9,10)]#Li
  

  
  # ====== convolutional filter to detect coverage steps ======
  Li_1     <- filter(adata2[,1],c(rep(1,filter_len),0,rep(-1,filter_len)),'convolution')
  Li_2     <- filter(adata2[,2],c(rep(1,filter_len),0,rep(-1,filter_len)),'convolution')
  Li_3     <- filter(adata2[,3],c(rep(1,filter_len),0,rep(-1,filter_len)),'convolution')
  KPDS_1   <- filter(adata[,1],c(rep(1,filter_len),0,rep(-1,filter_len)),'convolution')
  KPDS_2   <- filter(adata[,2],c(rep(1,filter_len),0,rep(-1,filter_len)),'convolution')
  KPDS_3   <- filter(adata[,3],c(rep(1,filter_len),0,rep(-1,filter_len)),'convolution')

  # ====== convolutional filter to normalise for upstream coverage ======
  norm_Li  <- filter(adata2[,1],c(rep(1,filter_len),0,rep(0,filter_len)),'convolution')
  norm_Li2 <- filter(adata2[,2],c(rep(1,filter_len),0,rep(0,filter_len)),'convolution')
  norm_Li3 <- filter(adata2[,3],c(rep(1,filter_len),0,rep(0,filter_len)),'convolution')
  norm_K   <- filter(adata[,1],c(rep(1,filter_len),0,rep(0,filter_len)),'convolution')
  norm_K2  <- filter(adata[,2],c(rep(1,filter_len),0,rep(0,filter_len)),'convolution')
  norm_K3  <- filter(adata[,3],c(rep(1,filter_len),0,rep(0,filter_len)),'convolution')


  # ====== calculate normalised convolved signal for both conditions to compare ======
  Li_1_r   <- Li_1 / norm_Li;
  Li_2_r   <- Li_2 / norm_Li2;
  Li_3_r   <- Li_3 / norm_Li3;
  KPDS_1_r <- KPDS_1 / norm_K;
  KPDS_2_r <- KPDS_2 / norm_K2;
  KPDS_3_r <- KPDS_3 / norm_K3;
  
  # ====== create data frames and handle zeros======
  Li_df <- as.data.frame(cbind(Li_1_r,Li_2_r,Li_3_r),stringsAsFactors=FALSE)
  Li_df$len <- rowSums((apply(Li_df,2,is.finite)))
  Li_1_r[which(!is.finite(Li_1_r))] <- 0 
  Li_2_r[which(!is.finite(Li_2_r))] <- 0 
  Li_3_r[which(!is.finite(Li_3_r))] <- 0 
  Li_123_r   <- (Li_1_r + Li_2_r + Li_3_r) / Li_df$len
  
  KPDS_df <- as.data.frame(cbind(KPDS_1_r,KPDS_2_r,KPDS_3_r),stringsAsFactors=FALSE)
  KPDS_df$len <- rowSums((apply(KPDS_df,2,is.finite)))  
  KPDS_1_r[which(!is.finite(KPDS_1_r))] <- 0
  KPDS_2_r[which(!is.finite(KPDS_2_r))] <- 0
  KPDS_3_r[which(!is.finite(KPDS_3_r))] <- 0
  KPDS_123_r <- (KPDS_1_r + KPDS_2_r + KPDS_3_r) / KPDS_df$len


  # ====== find local maxima in normalised convolved signals ======

  inds         <- which(Li_123_r >= filter_th_step & chr1$V7 > filter_n_border & chr1$V7 < (chr1$V3-chr1$V2 - filter_n_border) & rowMeans(adata2) >= (filter_k_cov/2))
  maxima       <- which(diff(sign(diff(Li_123_r)))==-2)+1
  inds_maxima  <- intersect(maxima, inds)

  indsK        <- which(KPDS_123_r >= filter_th_step & chr1$V7 > filter_n_border & chr1$V7 < (chr1$V3-chr1$V2 - filter_n_border) & rowMeans(adata2) >= (filter_k_cov/2))
  maximaK      <- which(diff(sign(diff(KPDS_123_r)))==-2)+1
  inds_maximaK <- intersect(maximaK, indsK)


  # --------------------------------------
  # --- score hits (K or KPDS scoring) ---
  # --------------------------------------
  if (length(inds_maximaK) > 0)
  {
  p_lm <- rep(2,length(inds_maximaK))
  for (f in 1:length(p_lm))
  {
    cov_inds1 <- which(adata2[inds_maximaK[f]+filter_len,] >= filter_k_cov)
    cov_inds2 <- which(adata[inds_maximaK[f]+filter_len,] >= filter_k_cov)
    v1 <- c(Li_1_r[inds_maximaK[f]], Li_2_r[inds_maximaK[f]], Li_3_r[inds_maximaK[f]])[cov_inds1]
    v2 <- c(KPDS_1_r[inds_maximaK[f]], KPDS_2_r[inds_maximaK[f]], KPDS_3_r[inds_maximaK[f]])[cov_inds2]
    v3 <- v1[is.finite(v1)]; v4 <- v2[is.finite(v2)]
    v3 <- v3[which(abs(v3 - median(v3)) <= 0.3)]; v4 <- v4[which(abs(v4 - median(v4)) <= 0.3)]
    if (length(v3) >= 1 && length(v4) >= 1 & mean(v4) > mean(v3))
    {
      if (length(unique(v3)) == 1 && length(unique(v4)) == 1)
      {
        a2 <- adata2[inds_maximaK[f]-filter_len,][cov_inds1]; b2 <- adata2[inds_maximaK[f]+filter_len,][cov_inds1]
        a <- adata[inds_maximaK[f]-filter_len,][cov_inds2];   b <- adata[inds_maximaK[f]+filter_len,][cov_inds2]
        p_lm[f] <- phyper(rowMeans(a),  rowMeans(a) + rowMeans(b), rowMeans(a2) + rowMeans(b2), rowMeans(a) + rowMeans(a2))
        
      }
      else
      {
        group <- as.factor(c(rep('v3',length(v3)),rep('v4',length(v4))))
        weight <- c(v3, v4); lm.D9 <- lm(weight ~ group)
        p_lm[f] <- anova(lm.D9)$'Pr(>F)'[1]
      }
    }
  }

  inds_maximaK_2 <- inds_maximaK[which(p_lm <= 1)]
  if (length(inds_maximaK_2) > 0)
  {
  tmp <- chr1[inds_maximaK_2,]
  chr <- rep(-1,0);start <- rep(-1,0);end <- rep(-1,0);name <- rep(-1,0);score <- rep(-1,0);strand <- rep(-1,0);count_i <- 1
  chr2 <- rep(-1,0);start2 <- rep(-1,0);end2 <- rep(-1,0);name2 <- rep(-1,0);score2 <- rep(-1,0);strand2 <- rep(-1,0);count_i2 <- 1
  hit_id <- rep(-1,0); hit_id2 <- rep(-1,0)
  
  for (t in 1:length(inds_maximaK_2))
  {
            agene <- chr1$V4[inds_maximaK_2[t]]
            inds_sub <- max(1,(inds_maximaK_2[t]-30)):inds_maximaK_2[t]
            genes <- chr1[inds_sub,4]
            inds_min <- head(inds_sub[which(genes == agene)],1)
            
            inds_sub2 <- max(1,(inds_maximaK_2[t]-50)):inds_maximaK_2[t]
	    genes2 <- chr1[inds_sub2,4]
            inds_min2 <- head(inds_sub2[which(genes2 == agene)],1)
            inds_sub2 <- inds_maximaK_2[t]:(inds_maximaK_2[t]+30)
	    genes2 <- chr1[inds_sub2,4]
            inds_max2 <- tail(inds_sub2[which(genes2 == agene)],1)
            
            # get the sequence and split
            
            if (chr1$V2[inds_maximaK_2[t]] != chr1$V2[inds_min] || chr1$V3[inds_maximaK_2[t]] != chr1$V3[inds_min])
            {
              chr   <- c(chr,chr1$V1[inds_min]);
              start <- c(start,chr1$V2[inds_min] + chr1$V7[inds_min] + 1);
              end   <- c(end,chr1$V3[inds_min]);
              name  <- c(name,t);
              score <- c(score,0);
              strand <- c(strand,'+')
              hit_id <- c(hit_id,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
              count_i <- count_i + 1
              chr  <- c(chr,chr1$V1[inds_min]);
	      start <- c(start,chr1$V2[inds_maximaK_2[t]]-1);
	      end   <- c(end,chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]]);
	      name  <- c(name,t)
	      score <- c(score,0)
              strand <- c(strand,'+')
              hit_id <- c(hit_id,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
              count_i <- count_i + 1
            } else {
              chr   <- c(chr,chr1$V1[inds_maximaK_2[t]])
	      start <- c(start,chr1$V2[inds_min] + chr1$V7[inds_min])
	      end   <- c(end,chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]])
	      name  <- c(name,t)
	      score <- c(score,0)
	      strand <- c(strand,'+')
              hit_id <- c(hit_id,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
              count_i <- count_i + 1
            }
            
	    # get the sequence for hairpin and split
            if (chr1$V2[inds_maximaK_2[t]] != chr1$V2[inds_min2] || chr1$V3[inds_maximaK_2[t]] != chr1$V3[inds_min2])
            {
              chr2   <- c(chr2,chr1$V1[inds_min2]);
              start2 <- c(start2,chr1$V2[inds_min2] + chr1$V7[inds_min2] + 1);
              end2   <- c(end2,chr1$V3[inds_min2]);
              name2  <- c(name2,t);
              score2 <- c(score2,0);
              strand2 <- c(strand2,'+')
              hit_id2 <- c(hit_id2,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
              count_i2 <- count_i2 + 1
              chr2  <- c(chr2,chr1$V1[inds_min2]);
	      start2 <- c(start2,chr1$V2[inds_maximaK_2[t]]-1);
	      end2   <- c(end2,chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]]);
	      name2  <- c(name2,t)
	      score2 <- c(score2,0)
              strand2 <- c(strand2,'+')
              hit_id2 <- c(hit_id2,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
              count_i2 <- count_i2 + 1
              if (chr1$V2[inds_maximaK_2[t]] != chr1$V2[inds_max2] || chr1$V3[inds_maximaK_2[t]] != chr1$V3[inds_max2])
              {
                #print(paste(t,"if max and min"))
                chr2   <- c(chr2,chr1$V1[inds_min2]);
                start2 <- c(start2,chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]]);
                end2   <- c(end2,chr1$V3[inds_maximaK_2[t]]);
                name2  <- c(name2,t);
                score2 <- c(score2,0);
                strand2 <- c(strand2,'+')
                hit_id2 <- c(hit_id2,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
                count_i2 <- count_i2 + 1
                chr2  <- c(chr2,chr1$V1[inds_min2]);
       	        start2 <- c(start2,chr1$V2[inds_max2]-1);
 	        end2   <- c(end2,chr1$V2[inds_max2] + chr1$V7[inds_max2]-1);
	        name2  <- c(name2,t)
	        score2 <- c(score2,0)
                strand2 <- c(strand2,'+')
                hit_id2 <- c(hit_id2,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
                count_i2 <- count_i2 + 1
              } else {
                #print(paste(t,"if min and not max"))
                chr2   <- c(chr2,chr1$V1[inds_maximaK_2[t]])
 	        start2 <- c(start2,chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]])
	        end2   <- c(end2,chr1$V2[inds_max2] + chr1$V7[inds_max2])
	        name2  <- c(name2,t)
	        score2 <- c(score2,0)
	        strand2 <- c(strand2,'+')
                hit_id2 <- c(hit_id2,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
                count_i2 <- count_i2 + 1                
              }
            } else {
              if (chr1$V2[inds_maximaK_2[t]] != chr1$V2[inds_max2] || chr1$V3[inds_maximaK_2[t]] != chr1$V3[inds_max2])
              {
                #print(paste(t,"if max and not min"))
	        chr2   <- c(chr2,chr1$V1[inds_min2]);
                start2 <- c(start2,chr1$V2[inds_min2] + chr1$V7[inds_min2]);
                end2   <- c(end2,chr1$V3[inds_maximaK_2[t]]);
                name2  <- c(name2,t);
                score2 <- c(score2,0);
                strand2 <- c(strand2,'+')
                hit_id2 <- c(hit_id2,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
                count_i2 <- count_i2 + 1
                chr2  <- c(chr2,chr1$V1[inds_min2]);
       	        start2 <- c(start2,chr1$V2[inds_max2]-1);
 	        end2   <- c(end2,chr1$V2[inds_max2] + chr1$V7[inds_max2]-1);
	        name2  <- c(name2,t)
	        score2 <- c(score2,0)
                strand2 <- c(strand2,'+')
                hit_id2 <- c(hit_id2,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
                count_i2 <- count_i2 + 1
              } else {
                chr2   <- c(chr2,chr1$V1[inds_maximaK_2[t]])
 	        start2 <- c(start2,chr1$V2[inds_min2] + chr1$V7[inds_min2])
	        end2   <- c(end2,chr1$V2[inds_max2] + chr1$V7[inds_max2])
	        name2  <- c(name2,t)
	        score2 <- c(score2,0)
	        strand2 <- c(strand2,'+')
                hit_id2 <- c(hit_id2,paste(chr1$V1[inds_maximaK_2[t]],chr1$V2[inds_maximaK_2[t]] + chr1$V7[inds_maximaK_2[t]],chr1$V6[inds_maximaK_2[t]], sep='_'))
                count_i2 <- count_i2 + 1                
              }
            }            
            tmp[t,] <- chr1[inds_min,]
  }
  
  tmp$V2 <- tmp$V2 + tmp$V7
  tmp$V3 <- chr1[inds_maximaK_2,2] +  chr1[inds_maximaK_2,7]
  tmp$V5 <- KPDS_123_r[inds_maximaK_2]
  tmp$V20 <- p_lm[which(p_lm <= 1)]
  tmp$V21 <- paste(tmp$V1,tmp$V3,tmp$V6, sep='_')
  #tmp <- tmp[which(p.adjust(tmp$V20,'fdr') <= 0.1),]
  write.table(tmp,paste(fname,'_maxima_step_10_',cond1,'_',cond2,'.lm.new.last.bed',sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
  write.table(tmp,paste('rG4/',fname,'_maxima_step_10_',cond1,'_',cond2,'.lm.new.last.bed',sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
  tmp2 <- as.data.frame(cbind(chr,start,end,name,score,strand,hit_id), stringsAsFactors=FALSE)
  tmp3 <- as.data.frame(cbind(chr2,start2,end2,name2,score2,strand2,hit_id2), stringsAsFactors=FALSE)

  write.table(tmp2,paste(fname,'_maxima_step_10_',cond1,'_',cond2,'.lm.new.last_seq.bed',sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
  write.table(tmp3,paste(fname,'_maxima_step_10_',cond1,'_',cond2,'.lm.new.last_hairpin.bed',sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
  write.table(tmp2,paste('rG4/',fname,'_maxima_step_10_',cond1,'_',cond2,'.lm.new.last_seq.bed',sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
  write.table(tmp3,paste('rG4/',fname,'_maxima_step_10_',cond1,'_',cond2,'.lm.new.last_hairpin.bed',sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
  
  
  }
  }

  # ------------------------------------------
  # --- score false positives (Li scoring) ---
  # ------------------------------------------
  if (length(inds_maxima) > 0)
  {
  p_lm <- rep(2,length(inds_maxima))
  for (f in 1:length(p_lm))
  {
      cov_inds1 <- which(adata2[inds_maxima[f]+filter_len,] >= filter_k_cov)
      cov_inds2 <- which(adata[inds_maxima[f]+filter_len,] >= filter_k_cov)
      v1 <- c(Li_1_r[inds_maxima[f]], Li_2_r[inds_maxima[f]], Li_3_r[inds_maxima[f]])[cov_inds1]
      v2 <- c(KPDS_1_r[inds_maxima[f]], KPDS_2_r[inds_maxima[f]], KPDS_3_r[inds_maxima[f]])[cov_inds2]
      v3 <- v1[is.finite(v1)]
      v4 <- v2[is.finite(v2)]
      v3 <- v3[which(abs(v3 - median(v3)) <= 0.3)]
      v4 <- v4[which(abs(v4 - median(v4)) <= 0.3)]

      if (length(v3) >= 1 && length(v4) >= 1 & mean(v4) < mean(v3))
      {
        if (length(unique(v3)) == 1 && length(unique(v4)) == 1)
        {
          a2 <- adata[inds_maxima[f]-filter_len,][cov_inds2]
          b2 <- adata[inds_maxima[f]+filter_len,][cov_inds2]
          a <- adata2[inds_maxima[f]-filter_len,][cov_inds1]
          b <- adata2[inds_maxima[f]+filter_len,][cov_inds1]
          p_lm[f] <- phyper(rowMeans(a),  rowMeans(a) + rowMeans(b), rowMeans(a2) + rowMeans(b2), rowMeans(a) + rowMeans(a2))
        }
        else
        {
          group <- as.factor(c(rep('v3',length(v3)),rep('v4',length(v4))))
          weight <- c(v3, v4)
          lm.D9 <- lm(weight ~ group)
          p_lm[f] <- anova(lm.D9)$'Pr(>F)'[1]
        }
        
      }
  }
  
  inds_maxima_2 <- inds_maxima[which(p_lm <= 1)]
  if (length(inds_maxima_2) > 0)
  {
  tmp <- chr1[inds_maxima_2,]
  for (t in 1:length(inds_maxima_2))
  {
          agene <- chr1$V4[inds_maxima_2[t]]
          inds_sub <- max(1,(inds_maxima_2[t]-30)):inds_maxima_2[t]
          genes <- chr1[inds_sub,4]
          inds_min <- head(inds_sub[which(genes == agene)],1)
          tmp[t,] <- chr1[inds_min,]
  }
  tmp$V2 <- tmp$V2 + tmp$V7
  tmp$V3 <- chr1[inds_maxima_2,2] +  chr1[inds_maxima_2,7]
  tmp$V5 <- Li_123_r[inds_maxima_2[which(p_lm <= 1)]]
  tmp$V20 <- p_lm[which(p_lm <= 1)]
  #tmp <- tmp[which(p.adjust(tmp$V20,'fdr') <= 0.1),]
  write.table(tmp,paste(fname,'_maxima_step_10_',cond2,'_',cond1,'.lm.new.last.bed',sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
  write.table(tmp,paste('rG4/',fname,'_maxima_step_10_',cond2,'_',cond1,'.lm.new.last.bed',sep=''), col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
  
  }
  }
