# TODO: change functions to accept either one IRinfo format or the other 
#       (gbfiles/dogma) and not be repeated.

#### This section contains the functions related to the gbfiles input ####

#' Find the genes near the junction coordinate based on the output of 
#' the \code{\link{gene.cordinate}} function.
#'
#' @param gene.cordinates cordinates given by the \code{\link{gene.cordinate}} 
#' function.
#' @param Irinfo vector with ira start, irb start, ir length and genome length.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively
#' @param radius number with the radius wanted
#' @param silence boolean
#' @return gene list
JunctRadiusGeneFinder<- function(gene.cordinates, IRinfo, J.pos, radius, silence=TRUE){
  if(J.pos==1){
    J<- IRinfo[1]
  }
  else if(J.pos==2){
    J<- IRinfo[1]+IRinfo[3]
  }
  else if(J.pos==3){
    J<- IRinfo[2]
  }
  else if(J.pos==4){
    J<- IRinfo[2]+IRinfo[3]
  }
  else {
    stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  g<-gene.cordinates
  r<- radius
  gs<- IRinfo[4]
  t1<- subset(g, J-r<as.numeric(g[,2]) & as.numeric(g[,2])<J+r)
  t2<- subset(g, J-r<as.numeric(g[,3]) & as.numeric(g[,3])<J+r)
  t3<- subset(g, J-r<as.numeric(g[,4]) & as.numeric(g[,4])<J+r)
  t4<- subset(g, J-r<as.numeric(g[,5]) & as.numeric(g[,5])<J+r)
  for (i in 1:nrow(g)){
    if(sum(g[i, ]=="0")==2){
      g[i,4:5]<-c(NA, NA)
    }
  }
  t5<- subset(g, J-r<(as.numeric(g[,2])+gs) & (as.numeric(g[,2])+gs)<J+r)
  t6<- subset(g, J-r<(as.numeric(g[,3])+gs) & (as.numeric(g[,3])+gs)<J+r)
  t7<- subset(g, J-r<(as.numeric(g[,4])+gs) & (as.numeric(g[,4])+gs)<J+r)
  t8<- subset(g, J-r<(as.numeric(g[,5])+gs) & (as.numeric(g[,5])+gs)<J+r)
  t<-unique(do.call("rbind", list(t1, t2, t3, t4, t5, t6, t7, t8)))
  if((length(t)+silence)==0){
    warning("There is no gene on the specified junction with the given radius, increasing the radius might help")
  }
  
  # If there are genes with the same name, we join them if they are in line
  if(nrow(t) > 0){
    n_occur <- data.frame(table(t[,1]))
    if(nrow(n_occur[n_occur$Freq > 1,]) > 0){
      df <- data.frame(t)
      # We takes those genes which are repeated
      df <- df[df$gene %in% n_occur[n_occur$Freq > 1, 'Var1'], ]
      df$start <- as.numeric(df$start)
      df$end <- as.numeric(df$end)
      df <- df[order(df$gene, df$end, df$start), ]
      
      # We calculate the differences between their end and the next start
      rows <- nrow(df)
      diffs <- rep(NA, rows)
      for(i in 1:rows){
        if(df$gene[(i)%%rows+1] == df$gene[i]){
          diffs[[(i)%%rows+1]] <- (df$start[(i)%%rows+1] - df$end[i])%%IRinfo[4]
        }
      }
      # We change the end value for the genes that are together
      df$interval <- diffs
      for(i in 1:nrow(df)){
        if(!is.na(df$interval[i]) && df$interval[i]<=1){
          if(i == 1 & df$gene[1] == df$gene[rows]){
            df$end[rows] <- df$end[1]
          } else if(df$gene[i] == df$gene[i-1]){
            df$end[i-1] <- df$end[i]
          }
        }
      }
      # We filter the genes that we don't need anymore
      df %<>% filter(!interval<=1) %>% select(-interval)
      t <- as.matrix(df)
    }
  }
  return(t)
}

#' Junction site gene plotter
#'
#' Plotting the genes in the vicinity of the junction site of the chloroplast.
#' 
#' @param Radius number with the radius wanted.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param track The track on which the genes are to be plotted, starting from the 
#' bottom to up as integers 1,2,...
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @param theme dataframe with the colors chosen on the web.
#' @return None.
JG.plotter<- function(Radius, J.pos, track, jlens, theme){
  if(J.pos==1){
    pc<-jlens$jlb.len
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-jlens$jsb.len
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-jlens$jsa.len
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-jlens$jla.len
    J<- IRList[[track]][2]+IRList[[track]][3]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, 
             or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  
  t<- JunctRadiusGeneFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup <- matrix(0, n, 4)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5,6)]
    }
    else {
      tup[i,]<-t[i, c(1,2,3,6)]
    }
  }
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  gcol<- function(tup){
    l<-length(tup[,1])
    col<-numeric(l)
    for (i in 1:l){
      if(str_detect(tup[i,1], "trn")){col[i]<- theme$trn.color}
      else if(str_detect(tup[i,1], "psb")){col[i]<- theme$psb.color}
      else if(str_detect(tup[i,1], "psa")){col[i]<- theme$psa.color}
      else if(str_detect(tup[i,1], "rps")){col[i]<- theme$rps.color}
      else if(str_detect(tup[i,1], "rrn")){col[i]<- theme$rrn.color}
      else if(str_detect(tup[i,1], "rpl")){col[i]<- theme$rpl.color}
      else if(str_detect(tup[i,1], "rbc")){col[i]<- theme$rbc.color}
      else if(str_detect(tup[i,1], "rpo")){col[i]<- theme$rpo.color}
      else if(str_detect(tup[i,1], "ycf")){col[i]<- theme$ycf.color}
      else if(str_detect(tup[i,1], "ndh")){col[i]<- theme$ndh.color}
      else if(str_detect(tup[i,1], "chl")){col[i]<- theme$chl.color}
      else if(str_detect(tup[i,1], "ccs")){col[i]<- theme$ccs.color}
      else if(str_detect(tup[i,1], "pet")){col[i]<- theme$pet.color}
      else if(str_detect(tup[i,1], "atp")){col[i]<- theme$atp.color}
      else {col[i]<- theme$other_gene.color}
    }
    return(col)
  }
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
  }
  for (i in 1:n){
    x1 <- max(Pcord[i,1], pc-10)
    x2 <- min(Pcord[i,2], pc+10)
    if(tup[i,4] == '-'){
      for (j in seq(0.10, 0.70, 0.05)){
        segments(x1, track*5+j+5, x2, track*5+j+5, lwd=1, col=paste(gcol(tup)[i]))
      }
    }
    else {
      for (j in seq(1.1, 1.7, 0.05)){
        segments(x1, track*5-j+5, x2, track*5-j+5, lwd=1, col=paste(gcol(tup)[i]))
      }
    }
  }
}


#' Find the mismatches near the junction cordinate based on the indel_table 
#' given by \code{\link{get_ir}} function.
#'
#' @param indel_table dataframe with mismatches given by \code{\link{get_ir}} 
#' function.
#' @param IRinfo vector with ira start, irb start, ir length and genome length.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively
#' @param radius number with the radius wanted
#' @return list with indel inside radius
JunctRadiusIndelFinder<- function(indel_table, IRinfo, J.pos, radius){
  if(J.pos==1){
    J<- IRinfo[1]
  } else if(J.pos==2){
    J<- IRinfo[1]+IRinfo[3]
  } else if(J.pos==3){
    J<- IRinfo[2]
  } else if(J.pos==4){
    J<- IRinfo[2]+IRinfo[3]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  
  inside <- vector()
  n <- nrow(indel_table)
  for(i in 1:n){
    indel_pos <- indel_table$position[i]
    if((J-radius < indel_pos) & (J+radius > indel_pos)
       | (J-radius < indel_pos + IRinfo[4]) & (J+radius > indel_pos + IRinfo[4])){
      inside <- c(inside, list(indel_table[indel_table$position == indel_pos, ]))
    }
  }
  if(length(inside)==0){
    return(NULL)
  }
  return(inside)
}

#' Junction site indel plotter
#' 
#' Plotting the indel mismatches.
#' 
#' @param Radius number with the radius wanted. 
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param track The track on which the genes are to be plotted, starting from the 
#' bottom to up as integers 1,2,...
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @param theme dataframe with the colors chosen on the web.
#' @return None.
JI.plotter<- function(Radius, J.pos, track, jlens, theme){
  indel_table <- IndelList[[track]]
  ir_info <- IRList[[track]]
  mid_inside <- NULL
  half_ir_plotdist <- (jlens$jsb.len - jlens$jlb.len)/2
  
  if(J.pos==1){
    pc<-jlens$jlb.len
    J<- ir_info[1]
    mid_inside <- subset(indel_table, position > (J + Radius) & position < (J+ir_info[3] - Radius))
  } else if(J.pos==2){
    pc<-jlens$jsb.len
    J<- ir_info[1]+ir_info[3]
  } else if(J.pos==3){
    pc<-jlens$jsa.len
    J<- ir_info[2]
    mid_inside <- subset(indel_table, position > (J + Radius) & position < (J+ir_info[3] - Radius))
  } else if(J.pos==4){
    pc<-jlens$jla.len
    J<- ir_info[2]+ir_info[3]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, 
             or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  
  if(NROW(mid_inside) > 0){
    mismatches <- unique(mid_inside$mismatch_type)
    txt <- NULL
    if(is.element('replace', mismatches)){
      nmis <- nrow(subset(mid_inside, mismatch_type == 'replace'))
      txt <- paste(txt, paste0(nmis, '/'))
    }
    if(is.element('delete', mismatches)){
      nmis <- nrow(subset(mid_inside, mismatch_type == 'delete'))
      if (is.null(txt)){
        txt <- paste(txt, paste0(nmis, '-'))
      } else {
        txt <- paste(txt, paste0(nmis, '-'), sep = ', ')
      }
    }
    if(is.element('insert', mismatches)){
      nmis <- nrow(subset(mid_inside, mismatch_type == 'insert'))
      if (is.null(txt)){
        txt <- paste(txt, paste0(nmis, '+'))
      } else {
        txt <- paste(txt, paste0(nmis, '+'), sep = ', ')
      }
    }
    points(half_ir_plotdist+pc, track*5+2.3+5, cex=0.7, 
           pch=6, col=theme$midmis.color)
    text(half_ir_plotdist+pc, track*5+4.5, txt, cex=0.25, font=4)
    segments(half_ir_plotdist+pc, track*5+7, half_ir_plotdist+pc,
             track*5+5.2, lty='dashed', lwd=1, col = theme$misline.color)
  } else if(J.pos == 1 | J.pos == 3) {
    text(half_ir_plotdist+pc, track*5+4.5, "//", cex=0.95, font=1)
  }
  
  inside <- JunctRadiusIndelFinder(indel_table, ir_info, J.pos, Radius)
  n <- length(inside)
  
  if(n > 0){
    for (i in 1:n){
      inside[[i]]$dista <- abs(inside[[i]]$position - J)
      if(J.pos == 1 | J.pos == 3){
        inside[[i]]$Pcord <- pc + (inside[[i]]$dista*half_ir_plotdist/Radius)
      } else {
        inside[[i]]$Pcord <- pc - (inside[[i]]$dista*half_ir_plotdist/Radius)
      }
      if(inside[[i]]$mismatch_type == 'replace'){
        inside[[i]]$pch <- 16
        inside[[i]]$col <- theme$replace.color
      } else if(inside[[i]]$mismatch_type == 'delete') {
        inside[[i]]$pch <- 17
        inside[[i]]$col <- theme$delete.color
      } else if(inside[[i]]$mismatch_type == 'insert') {
        inside[[i]]$pch <- 18
        inside[[i]]$col <- theme$insert.color
      } else {
        inside[[i]]$pch <- 19
      }
    }
    
    for (i in 1:n){
      points(max(inside[[i]]$Pcord, pc-10), track*5+2.3+5, cex=0.7, 
             pch=inside[[i]]$pch, col=inside[[i]]$col)
      # text(max(inside[[i]]$Pcord, pc-10), track*5+1.6+5, 
      #      paste(inside[[i]]$position, inside[[i]]$string), 
      #      cex=0.4, col='black', pos=3)
      segments(max(inside[[i]]$Pcord, pc-10), track*5+7, max(inside[[i]]$Pcord, pc-10), 
               track*5+5.2, lty='dashed', lwd=1, col = theme$misline.color)
    }
  }
}


#' Gene Name plotter
#'
#' Plotting the gene names on a given gene which is already plotted on the 
#' tracks of the IR plot.
#' 
#' @param Radius number with the radius wanted.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param track The track on which the genes are to be plotted, starting from the 
#' bottom to up as integers 1,2,...
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @param theme dataframe with the colors chosen on the web.
#' @return None.
GN.plotter<- function(Radius, J.pos, track, jlens, theme){  
  txtin.color <- theme$txtin.color
  txtout.color <- theme$txtout.color
  
  if(J.pos==1){
    pc<-jlens$jlb.len
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-jlens$jsb.len
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-jlens$jsa.len
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-jlens$jla.len
    J<- IRList[[track]][2]+IRList[[track]][3]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, 
             or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusGeneFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup <- matrix(0, n, 4)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5,6)]
    }
    else {
      tup[i,]<-t[i, c(1,2,3,6)]
    }
  }
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
  }
  for (i in 1:n){
    x1 <- max(Pcord[i,1], pc-10)
    x2 <- min(Pcord[i,2], pc+10)
    if(tup[i,4] == '-'){
      if (min(x1, x2) >= 104){
        text(103.5, track*5+1.05+5, tup[i,1], cex=txtcex, col=txtout.color, font=txtfont)
      }
      else if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5+0.35+5, tup[i,1], cex=txtcex, col=txtin.color, 
             font=txtfont)
        text(max(x1, x2)+1, track*5+0.3+5, paste(Rcord[i,2]-Rcord[i,1], "bp"), 
             cex=numcex, col=txtin.color, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+0.35+5, tup[i,1], cex=txtcex, 
             col=txtin.color, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+1.15+5, tup[i,1], cex=txtcex, 
             col=txtout.color, font=txtfont)
      }
    }
    else {
      if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5-1.40+5, tup[i,1], cex=txtcex, 
             col=txtin.color, font=txtfont)
        text(max(x1, x2)+1, track*5-1.45+5, paste(Rcord[i,2]-Rcord[i,1], "bp"), 
             cex=numcex, col=txtin.color, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-1.40+5, tup[i,1], cex=txtcex, 
             col=txtin.color, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-2.25+5, tup[i,1], cex=txtcex, 
             col=txtout.color, font=txtfont)
      }
    }
  }
}

#' On Junction plotter
#'
#' Plotting the fine tuned narrow lines showing the limits of the genes which 
#' are passing through the junction sites with their bp.
#' 
#' @param Radius number with the radius wanted.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param track The track on which the genes are to be plotted, starting from the 
#' bottom to up as integers 1,2,...
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @return None.
OJ.plotter <- function(Radius, J.pos, track, jlens){
  if(J.pos==1){
    pc<-jlens$jlb.len
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-jlens$jsb.len
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-jlens$jsa.len
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-jlens$jla.len
    J<- IRList[[track]][2]+IRList[[track]][3]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, 
             or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  
  t <- JunctRadiusGeneFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
  t[is.na(t)] <- "0"
  n <- length(t[,1])
  tup <- matrix(0, n, 4)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5,6)]
    }
    else {
      tup[i,]<-t[i, c(1,2,3,6)]
    }
  }
  bw <- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1] <- as.numeric(tup[,2])
  Rcord[,2] <- as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    Pcord[Pcord < 0] <- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
  }
  for (i in 1:n){
    x1 <- max(Pcord[i,1], pc-10)
    x2 <- min(Pcord[i,2], pc+10)
    if (Rcord[i,2] > J & Rcord[i,1] <J){
      if(tup[i,4] == '-'){
        Arrows(min(x1, x2), track*5+1.1+5, pc-0.15, track*5+1.1+5, 
               arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5+1.1+5, min(x1, x2), track*5+1.1+5, 
               arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5+1.1+5, max(x1, x2), track*5+1.1+5, 
               arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(max(x1, x2), track*5+1.1+5, pc+0.15, track*5+1.1+5, 
               arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5+1.4+5, paste(Rcord[i, 2]-J, "bp"), cex=0.4, pos=4)#up-right
        text(pc+1.4, track*5+1.4+5, paste(J-Rcord[i, 1], "bp"), cex=0.4, pos=2)#up-left
      } else {
        Arrows(max(x1, x2), track*5-2.1+5, pc+0.15, track*5-2.1+5, 
               arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5-2.1+5, max(x1, x2), track*5-2.1+5, 
               arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5-2.1+5, min(x1, x2), track*5-2.1+5, 
               arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(min(x1, x2), track*5-2.1+5, pc-0.15, track*5-2.1+5, 
               arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5-2.6+5, paste(Rcord[i, 2]-J, "bp"), cex=0.4, pos=4)#low-right
        text(pc+1.4, track*5-2.6+5, paste(J-Rcord[i, 1], "bp"), cex=0.4, pos=2)#low-left
      }
    }
  }
}

#' Junction Distance plotter
#'
#' plotting the narrow lines of the distance of the genes for the junction sites 
#' which are not passing through any gene and their bp.
#' 
#' @param Radius number with the radius wanted.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param track The track on which the genes are to be plotted, starting from the 
#' bottom to up as integers 1,2,...
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @return None.
JD.plotter <- function(Radius, J.pos, track, jlens){
  if(J.pos==1){
    pc<-jlens$jlb.len
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-jlens$jsb.len
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-jlens$jsa.len
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-jlens$jla.len
    J<- IRList[[track]][2]+IRList[[track]][3]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, 
             or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  t<- JunctRadiusGeneFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup <- matrix(0, n, 4)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5,6)]
    }
    else {
      tup[i,]<-t[i, c(1,2,3,6)]
    }
  }
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    ind<-Pcord < 0
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
    Rcord[ind]<- Rcord[ind] + IRList[[track]][4]
  }
  counter<- 0
  for (i in 1:n){
    if (Rcord[i,2] > J & Rcord[i,1] <J){
      counter<- counter+1
    }
  }
  if (counter==0){#find the closest gene to the junction site
    nearest<- which(abs(Pcord-pc)==min(abs(Pcord-pc)))
    col.cor<- floor((nearest-0.01)/n)+1
    row.cor<- nearest-n*floor((nearest-0.01)/n)###Now we have the row of the nearest gene
    #with this setting the position of the zero (genes tangent to the junction site) will not be plotted.
    #If interested either put "=" for the middle condition of the left or right binary operator or better develop zero only handling if function
    x1 <- max(Pcord[row.cor,1], pc-10)
    x2 <- min(Pcord[row.cor,2], pc+10)
    
    if (tup[row.cor, 4] == '-' & pc-Pcord[row.cor, col.cor] < 0 & x2 - x1 > 3){#top,right, big
      curvedarrow(from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), 
                  to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5+1.3+5), 
                  curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", 
                  arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5+1.3+0.2+5, 
           paste(Rcord[row.cor, 1]-J, "bp"), cex=0.4)
    }
    else if(tup[row.cor, 4] == '-' & pc-Pcord[row.cor, col.cor] < 0 & x2-x1 <= 3 ) {#top, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , 
             pc+(Pcord[row.cor, col.cor]-pc)/2 -2 , track*5+1.3+5, angle = 15, 
             length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, 
           paste(Rcord[row.cor, 1]-J, "bp"), cex=0.4)
    }
    else if(tup[row.cor, 4] == '+' & pc-Pcord[row.cor, col.cor] < 0 & x2-x1 > 3) {#low, right, big
      curvedarrow(from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), 
                   to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5-2.3+5), 
                   curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", 
                   arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5-2.3-0.2+5, 
           paste(Rcord[row.cor, 1]-J, "bp"), cex=0.4)
    }
    else if(tup[row.cor, 4] == '+' & pc-Pcord[row.cor, col.cor] < 0 & x2-x1 <= 3) {#low, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5, 
             pc+(Pcord[row.cor, col.cor]-pc)/2 -3, track*5-2.3+5, angle = 15, 
             length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 -4.5 , track*5-2.3-0.2+5, 
           paste(Rcord[row.cor, 1]-J, "bp"), cex=0.4)
    }
    else if(tup[row.cor, 4] == '+' & pc-Pcord[row.cor, col.cor] > 0 & x2-x1 > 3) {#low, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), 
                   to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5-2.3+5), 
                   curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", 
                   arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4, track*5-2.3-0.2+5, 
           paste(J-Rcord[row.cor, 2], "bp"), cex=0.4)
    }
    else if(tup[row.cor, 4] == '+' & pc-Pcord[row.cor, col.cor] > 0 & x2-x1 <= 3) {#low, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5, 
             pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5-2.3+5, angle = 15, 
             length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4.5 , track*5-2.3-0.2+5, 
           paste(J-Rcord[row.cor, 2], "bp"), cex=0.4)
    }
    else if(tup[row.cor, 4] == '-' & pc-Pcord[row.cor, col.cor] > 0 & x2-x1 > 3) {#top, left, big
      curvedarrow(from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), 
                  to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5+1.3+5), 
                  curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", 
                  arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, 
           paste(J-Rcord[row.cor, 2], "bp"), cex=0.4)
    }
    else if(tup[row.cor, 4] == '-' & pc-Pcord[row.cor, col.cor] > 0 & x2-x1 <= 3) {#top, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5, 
             pc+(Pcord[row.cor, col.cor]-pc)/2 + 2 , track*5+1.3+5, angle = 15, 
             length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5+1.3+0.2+5, 
           paste(J-Rcord[row.cor, 2], "bp"), cex=0.4)
    }
  }
}

#' Max radius finder
#'
#' Finding the best radius for the vicinity of a junction site. This radius 
#' determines which genes appear on the plot. The radius is decided according to
#' when a gene can be plotted. If the species are similar, the radius chosen 
#' will be the same for all of them. Else, each species will have its own radius 
#' and the plot won't be in scale.
#' 
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param l number of species to plot.
#' @param genelist list with the genes.
#' @param IRlist vector with ira start, irb start, ir length and genome length.
#' @return list with the radius for each species (the same if they are similar,
#' different if they are not similar species).
#' @export
Max.Radius <- function(J.pos, l, genelist, irlist){
  if(J.pos==1){ # TODO: Radius0 which one?
    Radius0<-550#680
  }
  else if(J.pos==2){
    Radius0<-100
  }
  else if(J.pos==3){
    Radius0<-700#1800
  }
  else if(J.pos==4){
    Radius0<-700#1000
  }
  R <- numeric(l)
  for (track in 1:l){
    Radius <- Radius0
    t<- JunctRadiusGeneFinder(genelist[[track]], irlist[[track]], J.pos , Radius)
    # if it doesn't find any gene, it tries with a bigger radius
    while(nrow(t)==0){
      Radius<- 1.2*Radius
      t<- JunctRadiusGeneFinder(genelist[[track]], irlist[[track]], J.pos , Radius)
    }
    R[track]<-Radius
  }
  # If the difference is not too big, choose the same radius for all
  if(max(R)-min(R) < 300){ # what difference is good?
    Radius <- round(max(R)+1)
    for(track in 1:l){
      R[track] <- Radius
    }
  }
  return(R)
}

#' Character count
#'
#' Counts how many characters a word has
#' 
#' @param word character.
#' @return number. TODO: ?
chr.count<- function(word){
  if (length(word)==1){
    return(length(strsplit(word, "")[[1]]))
  }
  else {
    t<- numeric(length(word))
    for (i in 1:length(word)){
      t[i]<- length(strsplit(word[i], "")[[1]])
    }
    return(t)
  }
}

#' LSC place
#'
#' Gives where LSC starts
#' 
#' @param IRinfo vector with ira start, irb start, ir length and genome length.
#' @return number: place where LSC starts.
LSC<- function(IRinfo){
  c<-as.vector(IRinfo)
  IR<- c[3]
  SSC<- c[2]-(c[1]+IR)
  return(c[4]-(SSC+2*IR))
}


#' SSC place
#'
#' Gives where SSC starts
#' 
#' @param IRinfo vector with ira start, irb start, ir length and genome length.
#' @return number: place where SSC starts.
SSC<- function(IRinfo){
  c<-as.vector(IRinfo)
  return(c[2]-(c[1]+c[3]))
}


#' Plotting data auxiliar function
#'
#' Calls all the needed functions to plot the information. It dies it for all
#' the species.
#' 
#' @param radius list with the calculated radius for the species and the junction.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param l number of species.
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @param theme dataframe with the colors chosen on the web.
#' @return None.
plot.data.aux <- function(radius, J.pos, l, jlens, theme){
  for (i in 1:l){JG.plotter(radius[i], J.pos, i, jlens, theme)}
  for (i in 1:l){GN.plotter(radius[i], J.pos, i, jlens, theme)}
  for (i in 1:l){OJ.plotter(radius[i], J.pos, i, jlens)}
  for (i in 1:l){JD.plotter(radius[i], J.pos, i, jlens)}
  for (i in 1:l){
    if(!is.na(IndelList[[i]])){
      JI.plotter(radius[i], J.pos, i, jlens, theme)
    }
  }
}

#' Plotting data function
#'
#' Creates the entire plot.
#' 
#' @param file filename where the plot is done.
#' @param theme dataframe with the colors chosen on the web.
#' @param sample boolean. TRUE if the plot is just a sample (2 species), FALSE
#' if all the species are to be plotted. 
#' @return None.
#' @export
IRs2<- function(file="IR_out", theme, sample = FALSE){
  if(sample == TRUE){
    l_aux <- 2
    spnames_aux <- spnames[1:2]
  } else {
    l_aux <- l
    spnames_aux <- spnames
  }
  
  lsc1.len <- 15
  ir.len <- 25
  ssc.len <- 30
  lsc2.len <- 11
  
  jlens <- list("jlb.len" = lsc1.len, "jsb.len" = lsc1.len + ir.len,
                "jsa.len" = lsc1.len + ir.len + ssc.len,
                "jla.len" = lsc1.len + 2*ir.len + ssc.len)
  
  names(FasList)<-  unlist(spnames)
  names(IRList)<-   unlist(spnames)
  names(GeneList)<- unlist(spnames)
  
  if(anyNA(IRList, recursive= TRUE) == FALSE) {
    # plot the lines to be filled with data
    m<-matrix(rep(c(lsc1.len, ir.len, ssc.len, ir.len, lsc2.len), l_aux+2), 5, l_aux+2)
    m[,l_aux+2]<-rep(NA,5)
    m[,1]<-rep(0,5)
    
    par(mai=c(0.5, 0.6+max(chr.count(unlist(spnames_aux)))/20, 0.7, 0.5))# c(bottom, left, top, right)
    barplot(m, horiz = T, lwd = 2, cex.names = 0.46, space = 4, border = T, axes = F, 
            col=c(theme$lsc.color, theme$ir.color, theme$ssc.color, 
                  theme$ir.color, theme$lsc.color))
    title(main="Inverted Repeats")
    
    # Separate each section with a line up down and name the junctions
    segments(jlens$jlb.len, 5*l_aux+8, jlens$jlb.len, 7, lty=1, lwd=1, col = theme$divline.color)
    segments(jlens$jsb.len, 5*l_aux+8, jlens$jsb.len, 7, lty=1, lwd=1, col = theme$divline.color)
    segments(jlens$jsa.len, 5*l_aux+8, jlens$jsa.len, 7, lty=1, lwd=1, col = theme$divline.color)
    segments(jlens$jla.len, 5*l_aux+8, jlens$jla.len, 7, lty=1, lwd=1, col = theme$divline.color)
    
    text(jlens$jlb.len, 5*l_aux+9, "JLB", font=4, cex=0.7)
    text(jlens$jsb.len, 5*l_aux+9, "JSB", font=4, cex=0.7)
    text(jlens$jsa.len, 5*l_aux+9, "JSA", font=4, cex=0.7)
    text(jlens$jla.len, 5*l_aux+9, "JLA", font=4, cex=0.7)
    
    for (i in seq(9.5, 5*l_aux+9, 5)){
      count<-which(seq(9.5, 5*l_aux+9, 5)==i)
      
      # sections names
      text(101, i, "LSC", cex=0.57, font=2, col=theme$txtin.color)
      text(2.5, i, "LSC", cex=0.57, font=2, col=theme$txtin.color)
      text(18, i, "IRb", cex=0.57, font=2, col=theme$txtin.color)
      text(43, i, "SSC", cex=0.57, font=2, col=theme$txtin.color)
      text(73, i, "IRa", cex=0.57, font=2, col=theme$txtin.color)
      
      # separation between sections
      points(27.5, i, cex=2.3, pch=18, col='white')
      points(55, i, cex=2.3, pch=18, col='white')
      points(82.5, i, cex=2.3, pch=18, col='white')
      
      text(55, i, "//", cex=0.95, font=1)
      if(is.na(IndelList[[count]])){
        text(27.5, i, "//", cex=0.95, font=1)
        text(82.5, i, "//", cex=0.95, font=1)
      }
      
      # bp count in each section
      LSCp<-paste(prettyNum(LSC(IRList[[count]]), big.mark = ","), "bp", "")
      SSCp<-paste(prettyNum(SSC(IRList[[count]]), big.mark = ","), "bp", "")
      IRp<- paste(prettyNum(IRList[[count]][3], big.mark = ","), "bp", "")
      
      text(10.5, i-0.09, LSCp, font=4, cex=0.52, col=theme$inbp.color)
      text(64.5, i-0.09, SSCp, font=4, cex=0.52, col=theme$inbp.color)
      text(35.5, i-0.09, IRp, font=4, cex=0.52, col=theme$inbp.color)
      text(90, i-0.09, IRp, font=4, cex=0.52, col=theme$inbp.color)
      
      # total bp count
      axis(2, i-1.5, labels = paste0(prettyNum(IRList[[count]][4], big.mark=","), "bp"), 
           font = 4, cex.axis=0.7, las=2, tick=F)
    }
    # species names
    axis(2, seq(9.5, 5*l_aux+9, 5), labels = spnames_aux, font = 4, cex.axis=0.7, las=2, tick=F)
    points(0, 4.5, cex=2.3, pch=18, col="white")
    
    # radius for each junction
    I<-   Max.Radius(1, l_aux, genelist = GeneList, irlist = IRList)
    II<-  Max.Radius(2, l_aux, genelist = GeneList, irlist = IRList)
    III<- Max.Radius(3, l_aux, genelist = GeneList, irlist = IRList)
    IV<-  Max.Radius(4, l_aux, genelist = GeneList, irlist = IRList)
    
    # plot genes and information around the junction
    plot.data.aux(I, 1, l_aux, jlens, theme)
    plot.data.aux(II, 2, l_aux, jlens, theme)
    plot.data.aux(III, 3, l_aux, jlens, theme)
    plot.data.aux(IV, 4, l_aux, jlens, theme)
    
  } else {
    # See for which samples IR was not found and print that on the plot.
    noIR <- c()
    for (i in 1:l){
      if(anyNA(IRList[i], recursive= TRUE) == TRUE){
        noIR <- append(noIR, paste('\n', paste0(i, '.'), names(IRList)[i]))
      }
    }
    
    text <- "No IR was found for the following samples:"
    for (i in 1:length(noIR)){
      text <- paste(text, noIR[i])
    }
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, text, cex = 1.6, col = "black")
    par(mar = c(5, 4, 4, 2) + 0.1)
  }
}

#### This section contains the functions related to the manual files input ####
# TODO: join this functions to the ones above.


#' #' Junction site gene plotter for dogma and fasta files
#'
#' Plotting the genes in the vicinity of the junction site of the chloroplast.
#' 
#' @param Radius number with the radius wanted.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param track The track on which the genes are to be plotted, starting from the 
#' bottom to up as integers 1,2,...
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @param theme dataframe with the colors chosen on the web.
#' @return None.
JG.plotterDinp<- function(Radius, J.pos, track, jlens, theme){
  if(J.pos==1){
    pc<-jlens$jlb.len
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-jlens$jsb.len
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-jlens$jsa.len
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-jlens$jla.len
    J<- IRListDinp[[track]][4]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, 
             or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  
  t<- JunctRadiusGeneFinder(GeneList[[track]], toIRinfo(IRListDinp[[track]]), J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)]
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  gcol<- function(tup){
    l<-length(tup[,1])
    col<-numeric(l)
    for (i in 1:l){
      if(str_detect(tup[i,1], "trn")){col[i]<- theme$trn.color}
      else if(str_detect(tup[i,1], "psb")){col[i]<- theme$psb.color}
      else if(str_detect(tup[i,1], "psa")){col[i]<- theme$psa.color}
      else if(str_detect(tup[i,1], "rps")){col[i]<- theme$rps.color}
      else if(str_detect(tup[i,1], "rrn")){col[i]<- theme$rrn.color}
      else if(str_detect(tup[i,1], "rpl")){col[i]<- theme$rpl.color}
      else if(str_detect(tup[i,1], "rbc")){col[i]<- theme$rbc.color}
      else if(str_detect(tup[i,1], "rpo")){col[i]<- theme$rpo.color}
      else if(str_detect(tup[i,1], "ycf")){col[i]<- theme$ycf.color}
      else if(str_detect(tup[i,1], "ndh")){col[i]<- theme$ndh.color}
      else if(str_detect(tup[i,1], "chl")){col[i]<- theme$chl.color}
      else if(str_detect(tup[i,1], "ccs")){col[i]<- theme$ccs.color}
      else if(str_detect(tup[i,1], "pet")){col[i]<- theme$pet.color}
      else if(str_detect(tup[i,1], "atp")){col[i]<- theme$atp.color}
      else {col[i]<- theme$other_gene.color}
    }
    return(col)
  }
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRListDinp[[track]][5])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      for (j in seq(0.10, 0.70, 0.05)){
        segments(x1, track*5+j+5, x2, track*5+j+5, lwd=1, col=paste(gcol(tup)[i]))
      }
    }
    else {
      x1 <- max(Pcord[i,1], pc-10)
      x2 <- min(Pcord[i,2], pc+10)
      for (j in seq(1.1, 1.7, 0.05)){
        segments(x1, track*5-j+5, x2, track*5-j+5, lwd=1, col=paste(gcol(tup)[i]))
      }
    }
  }
}

#' Junction site indel plotter for dogma and fasta files
#' 
#' Plotting the indel mismatches.
#' 
#' @param Radius number with the radius wanted. 
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param track The track on which the genes are to be plotted, starting from the 
#' bottom to up as integers 1,2,...
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @param theme dataframe with the colors chosen on the web.
#' @return None.
JI.plotterDinp<- function(Radius, J.pos, track, jlens, theme){
  indel_table <- IndelList[[track]]
  ir_info <- toIRinfo(IRListDinp[[track]])
  mid_inside <- NULL
  half_ir_plotdist <- (jlens$jsb.len - jlens$jlb.len)/2
  
  if(J.pos==1){
    pc<-jlens$jlb.len
    J<- ir_info[1]
    mid_inside <- subset(indel_table, position > (J + Radius) & position < (J+ir_info[3] - Radius))
  } else if(J.pos==2){
    pc<-jlens$jsb.len
    J<- ir_info[1]+ir_info[3]
  } else if(J.pos==3){
    pc<-jlens$jsa.len
    J<- ir_info[2]
    mid_inside <- subset(indel_table, position > (J + Radius) & position < (J+ir_info[3] - Radius))
  } else if(J.pos==4){
    pc<-jlens$jla.len
    J<- ir_info[2]+ir_info[3]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, 
             or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  
  if(!is.null(mid_inside)){
    mismatches <- unique(mid_inside$mismatch_type)
    txt <- NULL
    if(is.element('replace', mismatches)){
      nmis <- nrow(subset(mid_inside, mismatch_type == 'replace'))
      txt <- paste(txt, paste0(nmis, '/'))
    }
    if(is.element('delete', mismatches)){
      nmis <- nrow(subset(mid_inside, mismatch_type == 'delete'))
      if (is.null(txt)){
        txt <- paste(txt, paste0(nmis, '-'))
      } else {
        txt <- paste(txt, paste0(nmis, '-'), sep = ', ')
      }
    }
    if(is.element('insert', mismatches)){
      nmis <- nrow(subset(mid_inside, mismatch_type == 'insert'))
      if (is.null(txt)){
        txt <- paste(txt, paste0(nmis, '+'))
      } else {
        txt <- paste(txt, paste0(nmis, '+'), sep = ', ')
      }
    }
    points(half_ir_plotdist+pc, track*5+2.3+5, cex=0.7, 
           pch=6, col=theme$midmis.color)
    text(half_ir_plotdist+pc, track*5+4.5, txt, cex=0.25, font=4)
    segments(half_ir_plotdist+pc, track*5+7, half_ir_plotdist+pc,
             track*5+5.2, lty='dashed', lwd=1, col = theme$misline.color)
  } else if(J.pos == 1 | J.pos == 3) {
    text(half_ir_plotdist+pc, track*5+4.5, "//", cex=0.95, font=1)
  }
  
  inside <- JunctRadiusIndelFinder(indel_table, ir_info, J.pos, Radius)
  n <- length(inside)
  
  if(n > 0){
    for (i in 1:n){
      inside[[i]]$dista <- abs(inside[[i]]$position - J)
      if(J.pos == 1 | J.pos == 3){
        inside[[i]]$Pcord <- pc + (inside[[i]]$dista*half_ir_plotdist/Radius)
      } else {
        inside[[i]]$Pcord <- pc - (inside[[i]]$dista*half_ir_plotdist/Radius)
      }
      if(inside[[i]]$mismatch_type == 'replace'){
        inside[[i]]$pch <- 16
        inside[[i]]$col <- theme$replace.color
      } else if(inside[[i]]$mismatch_type == 'delete') {
        inside[[i]]$pch <- 17
        inside[[i]]$col <- theme$delete.color
      } else if(inside[[i]]$mismatch_type == 'insert') {
        inside[[i]]$pch <- 18
        inside[[i]]$col <- theme$insert.color
      } else {
        inside[[i]]$pch <- 19
      }
    }
    
    for (i in 1:n){
      points(max(inside[[i]]$Pcord, pc-10), track*5+2.3+5, cex=0.7, 
             pch=inside[[i]]$pch, col=inside[[i]]$col)
      # text(max(inside[[i]]$Pcord, pc-10), track*5+1.6+5, 
      #      paste(inside[[i]]$position, inside[[i]]$string), 
      #      cex=0.4, col='black', pos=3)
      segments(max(inside[[i]]$Pcord, pc-10), track*5+7, max(inside[[i]]$Pcord, pc-10), 
               track*5+5.2, lty='dashed', lwd=1, col = theme$misline.color)
    }
  }
}

#' Gene Name plotter
#'
#' Plotting the gene names on a given gene which is already plotted on the 
#' tracks of the IR plot.
#' 
#' @param Radius number with the radius wanted.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param track The track on which the genes are to be plotted, starting from the 
#' bottom to up as integers 1,2,...
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @param theme dataframe with the colors chosen on the web.
#' @return None.
GN.plotterDinp<- function(Radius, J.pos, track, jlens, theme){
  txtin.color <- theme$txtin.color
  txtout.color <- theme$txtout.color
  
  if(J.pos==1){
    pc<-jlens$jlb.len
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-jlens$jsb.len
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-jlens$jsa.len
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-jlens$jla.len
    J<- IRListDinp[[track]][4]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, 
             or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusGeneFinder(GeneList[[track]], toIRinfo(IRListDinp[[track]]), J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)]
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRListDinp[[track]][5])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (min(x1, x2) >= 104){
        text(103.5, track*5+1.05+5, tup[i,1], cex=txtcex, 
             col=txtout.color, font=txtfont)
      }
      else if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5+0.35+5, tup[i,1], cex=txtcex, 
             col=txtin.color, font=txtfont)
        text(max(x1, x2)+1, track*5+0.3+5, paste(Rcord[i,1]-Rcord[i,2], "bp"),
             cex=numcex, col=txtin.color, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+0.35+5, tup[i,1], 
             cex=txtcex, col=txtin.color, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+1.15+5, tup[i,1], 
             cex=txtcex, col=txtout.color, font=txtfont)
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5-1.40+5, tup[i,1], 
             cex=txtcex, col=txtin.color, font=txtfont)
        text(max(x1, x2)+1, track*5-1.45+5, paste(Rcord[i,2]-Rcord[i,1], "bp"), 
             cex=numcex, col=txtin.color, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-1.40+5, tup[i,1], 
             cex=txtcex, col=txtin.color, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-2.25+5, tup[i,1], 
             cex=txtcex, col=txtout.color, font=txtfont)
      }
    }
  }
}

#' On Junction plotter
#'
#' Plotting the fine tuned narrow lines showing the limits of the genes which 
#' are passing through the junction sites with their bp.
#' 
#' @param Radius number with the radius wanted.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param track The track on which the genes are to be plotted, starting from the 
#' bottom to up as integers 1,2,...
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @return None.
OJ.plotterDinp<- function(Radius, J.pos, track, jlens){#GeneName plotter
  if(J.pos==1){
    pc<-jlens$jlb.len
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-jlens$jsb.len
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-jlens$jsa.len
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-jlens$jla.len
    J<- IRListDinp[[track]][4]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, 
             or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  
  t<- JunctRadiusGeneFinder(GeneList[[track]], toIRinfo(IRListDinp[[track]]), J.pos , Radius)
  t[is.na(t)]<- "0"
  n<- length(t[,1])
  tup<- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)]
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRListDinp[[track]][5])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (Rcord[i,1] > J & Rcord[i,2] <J ){
        Arrows(min(x1, x2), track*5+1.1+5, pc-0.15, track*5+1.1+5, arr.type = "T", 
               cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5+1.1+5, min(x1, x2), track*5+1.1+5, arr.type = "T", 
               cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5+1.1+5, max(x1, x2), track*5+1.1+5, arr.type = "T", 
               cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(max(x1, x2), track*5+1.1+5, pc+0.15, track*5+1.1+5, arr.type = "T", 
               cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5+1.4+5, paste(Rcord[i, 1]-J+1, "bp"), cex=0.4, pos=4)#up-right
        text(pc+1.4, track*5+1.4+5, paste(J-Rcord[i, 2], "bp"), cex=0.4, pos=2)#up-left
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (Rcord[i,2] > J & Rcord[i,1] <J ){
        Arrows(max(x1, x2), track*5-2.1+5, pc+0.15, track*5-2.1+5, arr.type = "T", 
               cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5-2.1+5, max(x1, x2), track*5-2.1+5, arr.type = "T", 
               cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5-2.1+5, min(x1, x2), track*5-2.1+5, arr.type = "T", 
               cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(min(x1, x2), track*5-2.1+5, pc-0.15, track*5-2.1+5, arr.type = "T", 
               cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5-2.6+5, paste(Rcord[i, 2]-J, "bp"), cex=0.4, pos=4)#low-right
        text(pc+1.4, track*5-2.6+5, paste(J-Rcord[i, 1], "bp"), cex=0.4, pos=2)#low-left
      }
    }
  }
}

#' Junction Distance plotter
#'
#' plotting the narrow lines of the distance of the genes for the junction sites 
#' which are not passing through any gene and their bp.
#' 
#' @param Radius number with the radius wanted.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param track The track on which the genes are to be plotted, starting from the 
#' bottom to up as integers 1,2,...
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @return None.
JD.plotterDinp<- function(Radius, J.pos, track, jlens){
  if(J.pos==1){
    pc<-jlens$jlb.len
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-jlens$jsb.len
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-jlens$jsa.len
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-jlens$jla.len
    J<- IRListDinp[[track]][4]
  } else {
    stop("J.pos missing or out bound. It should be either 1,2,3, 
             or 4 for JLB, JSB, JSA, and JLA, recpectively ")
  }
  
  t<- JunctRadiusGeneFinder(GeneList[[track]], toIRinfo(IRListDinp[[track]]), J.pos , Radius)
  t[is.na(t)]<- "0"
  n <- length(t[,1])
  tup <- matrix(0, n, 3)
  for (i in 1:n){
    dist<-abs(as.numeric(t[i,][2:5])-J)
    if (which(dist==min(dist))>3){
      tup[i,]<-t[i, c(1,4,5)]
    }
    else {
      tup[i,]<-t[i, c(1,2,3)]
    }
  }
  bw<- 10/Radius
  Rcord<- matrix(0, n, 2)
  Rcord[,1]<-as.numeric(tup[,2])
  Rcord[,2]<-as.numeric(tup[,3])
  Pcord<- (Rcord-J)*bw+pc
  if (J.pos==4){
    ind<-Pcord < 0
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRListDinp[[track]][5])-J)*bw+pc
    Rcord[ind]<- Rcord[ind] + IRListDinp[[track]][5]
  }
  counter<- 0
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (Rcord[i,1] > J & Rcord[i,2] <J ){
        counter<- counter+1
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (Rcord[i,2] > J & Rcord[i,1] <J ){
        counter<- counter+1
      }
    }
  }
  if (counter==0){#find the closest gene to the junction site
    nearest<- which(abs(Pcord-pc)==min(abs(Pcord-pc)))
    col.cor<- floor((nearest-0.01)/n)+1
    row.cor<- nearest-n*floor((nearest-0.01)/n)###Now we have the row of the nearest gene
    #with this setting the position of the zero (genes tangent to the junction site) 
    #will not be plotted. If interested either put "=" for the middle condition of 
    #the left or right binary operator or better develop zero only handling if function
    if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & 
       min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top,right, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), 
                   to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5+1.3+5), 
                   curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", 
                   arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5+1.3+0.2+5, 
           paste(Rcord[row.cor, 2]-J, "bp"), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 
            & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5, 
             pc+(Pcord[row.cor, col.cor]-pc)/2 -2, track*5+1.3+5, angle = 15, 
             length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5, 
           track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp"), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 
            & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, right, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), 
                   to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5-2.3+5), 
                   curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", 
                   arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5-2.3-0.2+5, 
           paste(Rcord[row.cor, 1]-J, "bp"), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 
            & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , 
             pc+(Pcord[row.cor, col.cor]-pc)/2 -3 , track*5-2.3+5, angle = 15, 
             length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 -4.5 , track*5-2.3-0.2+5, 
           paste(Rcord[row.cor, 1]-J, "bp"), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 
            & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), 
                   to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5-2.3+5), 
                   curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", 
                   arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4, track*5-2.3-0.2+5, 
           paste(J-Rcord[row.cor, 2], "bp"), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & 
            abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , 
             pc+(Pcord[row.cor, col.cor]-pc)/2 + 3 , track*5-2.3+5, angle = 15, 
             length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4.5 , track*5-2.3-0.2+5, 
           paste(J-Rcord[row.cor, 2], "bp"), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & 
            min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), 
                   to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5+1.3+5), 
                   curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", 
                   arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, 
           paste(J-Rcord[row.cor, 1], "bp"), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & 
            min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , 
             pc+(Pcord[row.cor, col.cor]-pc)/2 + 2 , track*5+1.3+5, 
             angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4 , track*5+1.3+0.2+5, 
           paste(J-Rcord[row.cor, 1], "bp"), cex=0.4)
    }
  }
}

#' Max radius finder
#'
#' Finding the best radius for the vicinity of a junction site. This radius 
#' determines which genes appear on the plot. The radius is decided according to
#' when a gene can be plotted. If the species are similar, the radius chosen 
#' will be the same for all of them. Else, each species will have its own radius 
#' and the plot won't be in scale.
#' 
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param l number of species to plot.
#' @param genelist list with the genes.
#' @param IRlist vector with ira start, irb start, ir length and genome length.
#' @return list with the radius for each species (the same if they are similar,
#' different if they are not similar species).
#' @export
Max.RadiusDinp<-function(J.pos, l, genelist, irlistDinp){
  if(J.pos==1){
    Radius0<-680#680#100
  }
  if(J.pos==2){
    Radius0<-100#100#100
  }
  if(J.pos==3){
    Radius0<-1000#1800#100
  }
  if(J.pos==4){
    Radius0<-1000#1000#100
  }
  R<- numeric(l)
  for (track in 1:l){
    Radius <- Radius0
    t<- JunctRadiusGeneFinder(genelist[[track]], toIRinfo(irlistDinp[[track]]), J.pos , Radius)
    while(nrow(t)==0){
      Radius<- Radius+(0.2)*Radius
      t<- JunctRadiusGeneFinder(genelist[[track]], toIRinfo(irlistDinp[[track]]), J.pos , Radius)
    }
    R[track]<-Radius
  }
  # If the difference is not too big, choose the same radius for all
  if(max(R)-min(R) < 300){ # TODO: decide what difference is better
    Radius <- round(max(R)+1)
    for(track in 1:l){
      R[track] <- Radius
    }
  }
  return(R)
}


#' LSC place
#'
#' Gives where LSC starts
#' 
#' @param IRinfo vector with ira start, irb start, ir length and genome length.
#' @return number: place where LSC starts.
LSCDinp<- function(IRinp){#a vecotr as a list
  c<-as.vector(IRinp)
  return(c[5]-((c[2]-c[1])+(c[4]-c[3])+(c[3]-c[2])))
}


#' SSC place
#'
#' Gives where SSC starts
#' 
#' @param IRinfo vector with ira start, irb start, ir length and genome length.
#' @return number: place where SSC starts.
SSCDinp<- function(IRinp){
  c<-as.vector(IRinp)
  return(c[3]-c[2])
}

#' Plotting data auxiliar function
#'
#' Calls all the needed functions to plot the information. It dies it for all
#' the species.
#' 
#' @param radius list with the calculated radius for the species and the junction.
#' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, 
#' and 4 respectively.
#' @param l number of species.
#' @param jlens named list with the distance in the plot for the junction sites, 
#' in the same order as J.pos: jlb.len, jsb.len, jsa.len, jla.len.
#' @param theme dataframe with the colors chosen on the web.
#' @return None.
plot.data.aux.D <- function(radius, J.pos, l, jlens, theme){
  for (i in 1:l){JG.plotterDinp(radius[i], J.pos, i, jlens, theme)}
  for (i in 1:l){GN.plotterDinp(radius[i], J.pos, i, jlens, theme)}
  for (i in 1:l){OJ.plotterDinp(radius[i], J.pos, i, jlens)}
  for (i in 1:l){JD.plotterDinp(radius[i], J.pos, i, jlens)}
  for (i in 1:l){
    if(!is.na(IndelList[[i]])){
      JI.plotterDinp(radius[i], J.pos, i, jlens, theme)
    }
  }
}

#' Plotting data function for dogma and fasta files
#'
#' Creates the entire plot.
#' 
#' @param file filename where the plot is done.
#' @param theme dataframe with the colors chosen on the web.
#' @param sample boolean. TRUE if the plot is just a sample (2 species), FALSE
#' if all the species are to be plotted. 
#' @return None.
#' @export
IRsD2<- function(file="IR_out", theme, sample = FALSE){
  if(sample == TRUE){
    l_aux <- 2
    spnames_aux <- spnames[1:2]
  } else {
    l_aux <- l
    spnames_aux <- spnames
  }
  
  lsc1.len <- 15
  ir.len <- 25
  ssc.len <- 30
  lsc2.len <- 11
  
  jlens <- list("jlb.len" = lsc1.len, "jsb.len" = lsc1.len + ir.len,
                "jsa.len" = lsc1.len + ir.len + ssc.len,
                "jla.len" = lsc1.len + 2*ir.len + ssc.len)
  
  names(FasList)<- unlist(spnames)
  names(IRListDinp)<- unlist(spnames)
  names(GeneList)<- unlist(spnames)
  
  if(anyNA(IRListDinp, recursive= TRUE) == FALSE) {
    # plot the lines to be filled with data
    m<-matrix(rep(c(lsc1.len, ir.len, ssc.len, ir.len, lsc2.len), l_aux+2), 5, l_aux+2)
    m[,l+2]<-rep(NA,5)
    m[,1]<-rep(0,5)
    
    par(mai=c(0.5, 0.6+max(chr.count(unlist(spnames_aux)))/20, 0.7, 0.5))# c(bottom, left, top, right)
    barplot(m, horiz = T, lwd = 2, cex.names = 0.46, space = 4, border = T, axes = F, 
            col=c(theme$lsc.color, theme$ir.color, theme$ssc.color, 
                  theme$ir.color, theme$lsc.color))
    title(main="Inverted Repeats")
    
    # Separate each section with a line up down and name the junctions
    segments(jlens$jlb.len, 5*l_aux+8, jlens$jlb.len, 7, lty=1, lwd=1, col = theme$divline.color)
    segments(jlens$jsb.len, 5*l_aux+8, jlens$jsb.len, 7, lty=1, lwd=1, col = theme$divline.color)
    segments(jlens$jsa.len, 5*l_aux+8, jlens$jsa.len, 7, lty=1, lwd=1, col = theme$divline.color)
    segments(jlens$jla.len, 5*l_aux+8, jlens$jla.len, 7, lty=1, lwd=1, col = theme$divline.color)
    
    text(jlens$jlb.len, 5*l_aux+9, "JLB", font=4, cex=0.7)
    text(jlens$jsb.len, 5*l_aux+9, "JSB", font=4, cex=0.7)
    text(jlens$jsa.len, 5*l_aux+9, "JSA", font=4, cex=0.7)
    text(jlens$jla.len, 5*l_aux+9, "JLA", font=4, cex=0.7)
    
    for (i in seq(9.5, 5*l_aux+9, 5)){
      count<-which(seq(9.5, 5*l_aux+9, 5)==i)
      
      # sections names
      text(101, i, "LSC", cex=0.57, font=2, col=theme$txtin.color)
      text(2.5, i, "LSC", cex=0.57, font=2, col=theme$txtin.color)
      text(18, i, "IRb", cex=0.57, font=2, col=theme$txtin.color)
      text(43, i, "SSC", cex=0.57, font=2, col=theme$txtin.color)
      text(73, i, "IRa", cex=0.57, font=2, col=theme$txtin.color)
      
      # separation between sections
      points(27.5, i, cex=2.3, pch=18, col='white')
      points(55, i, cex=2.3, pch=18, col='white')
      points(82.5, i, cex=2.3, pch=18, col='white')
      
      text(55, i, "//", cex=0.95, font=1)
      if(is.na(IndelList[[count]])){
        text(27.5, i, "//", cex=0.95, font=1)
        text(82.5, i, "//", cex=0.95, font=1)
      }
      
      # bp count in each section
      LSCp<-paste(prettyNum(LSCDinp(IRListDinp[[count]]), big.mark = ","), "bp", "")
      SSCp<-paste(prettyNum(SSCDinp(IRListDinp[[count]]), big.mark = ","), "bp", "")
      IRpb<- paste(prettyNum((IRListDinp[[count]][2]-IRListDinp[[count]][1]), big.mark = ","), "bp", "")
      IRpa<- paste(prettyNum((IRListDinp[[count]][4]-IRListDinp[[count]][3]), big.mark = ","), "bp", "")
      
      text(10.5, i-0.09, LSCp, font=4, cex=0.52, col=theme$inbp.color)
      text(64.5, i-0.09, SSCp, font=4, cex=0.52, col=theme$inbp.color)
      text(35.5, i-0.09, IRpb, font=4, cex=0.52, col=theme$inbp.color)
      text(90, i-0.09, IRpa, font=4, cex=0.52, col=theme$inbp.color)
      
      # total bp count
      axis(2, i-1.5, labels = paste0(prettyNum(IRListDinp[[count]][4], big.mark=","), "bp"), 
           font = 4, cex.axis=0.7, las=2, tick=F)
    }
    
    # species names
    axis(2, seq(9.5, 5*l_aux+9, 5), labels = spnames_aux, font = 4, cex.axis=0.7, las=2, tick=F)
    points(0, 4.5, cex=2.3, pch=18, col="white")
    
    # radius for each junction
    I<-   Max.RadiusDinp(1, l, genelist = GeneList, irlistDinp = IRListDinp)
    II<-  Max.RadiusDinp(2, l, genelist = GeneList, irlistDinp = IRListDinp)
    III<- Max.RadiusDinp(3, l, genelist = GeneList, irlistDinp = IRListDinp)
    IV<-  Max.RadiusDinp(4, l, genelist = GeneList, irlistDinp = IRListDinp)
    
    # plot genes and information around the junction
    plot.data.aux.D(I, 1, l, jlens, theme)
    plot.data.aux.D(II, 2, l, jlens, theme)
    plot.data.aux.D(III, 3, l, jlens, theme)
    plot.data.aux.D(IV, 4, l, jlens, theme)
    
  } else {
    # See for which samples IR was not found and print that on the plot.
    noIR <- c()
    for (i in 1:l){
      if(anyNA(IRListDinp[i], recursive= TRUE) == TRUE){
        noIR <- append(noIR, paste('\n', paste0(i, '.'), names(IRListDinp)[i]))
      }
    }
    
    text <- "No IR was found for the following samples:"
    for (i in 1:length(noIR)){
      text <- paste(text, noIR[i])
    }
    
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, text, cex = 1.6, col = "black")
    par(mar = c(5, 4, 4, 2) + 0.1)
  }
}##For the Fasta and Dogma Files
