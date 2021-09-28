###Under the license of the University of Helsinki
#Author: Ali Amiryousefi
#email: ali.amiryousefi@helisnki.fi


#loading required packages

# library(seqinr)
# library(ape)
# library(reutils)
# library(snow)
# library(snowfall)
# library(knitr)
# library(magrittr)
# library(dplyr)

library(shape)
library(jpeg)
library(diagram)
library(grDevices)
library(gplots)


gene.name<-function(gb, type){
  if(type=="gene"){
    t <- gb[grep("  gene  ", gb)+1]
    if(length(t) > 0){
      for (i in 1:length(t)){
        crude<-strsplit(t[i], " ")[[1]]
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  gene  ", gb)+2][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  gene  ", gb)+3][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  gene  ", gb)+4][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        crude.name<-crude[which(crude!=rep("", length(crude)))]
        t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
      }
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else if(type=="tRNA"){
    t <- gb[grep(" tRNA  ", gb)+1]
    if(length(t) > 0){
      for (i in 1:length(t)){
        crude<-strsplit(t[i], " ")[[1]]
        if(length(grep("=", crude))==0){
          t[i] <- gb[grep("tRNA  ", gb)+2][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  tRNA  ", gb)+3][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  tRNA  ", gb)+4][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        crude.name<-crude[which(crude!=rep("", length(crude)))]
        t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
      }
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("tRNA ", gb)+1]), ""), ""), "\n")
      }
    }
  }
  else if(type=="rRNA"){
    t <- gb[grep("rRNA  ", gb)+1]
    if(length(t) > 0){
      for (i in 1:length(t)){
        crude<-strsplit(t[i], " ")[[1]]
        if(length(grep("=", crude))==0){
          t[i] <- gb[grep("rRNA  ", gb)+2][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  rRNA  ", gb)+3][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  rRNA  ", gb)+4][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        crude.name<-crude[which(crude!=rep("", length(crude)))]
        t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
      }
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else if(type=="mRNA"){
    t <- gb[grep("mRNA  ", gb)+1]
    if(length(t) > 0){
      for (i in 1:length(t)){
        crude<-strsplit(t[i], " ")[[1]]
        if(length(grep("=", crude))==0){
          t[i] <- gb[grep("mRNA  ", gb)+2][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  mRNA  ", gb)+3][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
          t[i] <- gb[grep("  mRNA  ", gb)+4][i]
          crude<- strsplit(t[i], " ")[[1]]
        }
        crude.name<-crude[which(crude!=rep("", length(crude)))]
        t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
      }
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else {
    stop("The type should be defined as either gene or tRNA")
  }
  t<-t[!is.na(t)]
  return(t)
  #intermediate gene name function to substract the gene names of either gene or tRNA, mRNA or rRNA from their second line information(or third), the input is the gb file.
}


gb.gene.cor<- function(gb){
  gene<- gb[grep("  gene  ", gb)]
  trna<- gb[grep("  tRNA ", gb)]
  rrna<- gb[grep("  rRNA ", gb)]
  t<-c(gene, trna, rrna)
  m<- matrix(0, length(t), 2)
  for (i in 1:length(t)){
    if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]==","){#check the last characted if its ',' and add the next one from gb file to that
      #previous code--> if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]=="," && strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])-1]==")")
      s<-t[i]
      t[i]<-paste(t[i], gsub(" ", "", gb[which(t[i]==gb)+1]), "")
    }
    t[i]<- gsub(", ", ",", t[i])
    t[i]<- gsub(") ", ")", t[i])
    if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]==","){#check again if the last characted if its ',' and add the next one from gb file to that
      #previous code--> if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]=="," && strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])-1]==")")
      t[i]<-paste(t[i], gsub(" ", "", gb[which(s==gb)+2]), "")
    }
    t[i]<- gsub(", ", ",", t[i])
    t[i]<- gsub(") ", ")", t[i])
    m[i,2]<-strsplit(t[i], " ")[[1]][length(strsplit(t[i], " ")[[1]])]
  }
  names<-  c(gene.name(gb, "gene"), gene.name(gb, "tRNA"), gene.name(gb, "rRNA") )
  if(length(t)!=length(names)){
    stop("Error: Check the gene names of the gb file")
  }
  m[,1]<-names
  m #gives the crude table extracted from gb file for genes
}


gene.cordinates<- function(gb){
  gb.gene.cor.out<-gb.gene.cor(gb)
  l<- length(gb.gene.cor.out[,1])
  m<- matrix(0, l, 5)
  m[,1]<- gb.gene.cor.out[,1]
  for (i in 1:l){
    if(length(strsplit(gb.gene.cor.out[i,2], ",")[[1]])==2){
      m[i,2]<-sub("join\\(", "", strsplit(gb.gene.cor.out[i,2], ",")[[1]][1])
      m[i,2]<-sub("order\\(", "", m[i,2])
      m[i,4]<-sub(")", "", strsplit(gb.gene.cor.out[i,2], ",")[[1]][2])
      cord<-m[i,2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
      cord<-m[i,4]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(4,5)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(4,5)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
    else {
      m[i,2]<-gb.gene.cor.out[i,2]
      cord<-m[i,2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        print(cord)
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
  }
  m<- gsub(")", "", m)
  m<- gsub("\\(", "", m)
  m<- gsub("<", "", m)
  m<- gsub(">", "", m)
  return(m)#gives a polished table of genes and their cordinates, as the gene names their start and end of first part and the second parts as the first to fifth columns. The orders are reflected with the swapping of values
}


JunctRadiusFinder<- function(gene.cordinates, IRinfo, J.pos, radius, silence=TRUE){#function to find the genes near the junction cordinate based on the output of the gene.cordinate function
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
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  g<-gene.cordinates
  l<- ncol(g)
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
  return(t)
}


JunctNoFinder<- function(gene.cordinates, IRinfo, J.pos ,Number){#the same function like it precusor but returns the number of genes found near the junction site
  g<-gene.cordinates
  l<- ncol(g)
  n<-Number
  counter=0
  count=0
  while (nrow(JunctRadiusFinder(g, IRinfo, J.pos ,count)) <= n-1){
    count<-count+1
  }
  t<-JunctRadiusFinder(g, IRinfo, J.pos ,count)
  t[is.na(t)]<- "0"
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
  return(tup)
}


JG.plotter<- function(Radius, J.pos, track){#JunctionGene plotter
  #' Junction site gene plotter
  #'
  #' Plotting the genes in the vicinity of the junction site of the chloroplast.
  #' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, and 4 respectively
  #' @param margin The margin value to be added to the interval as the proportion of the minimum lenght of the nth gene in vicinity
  #' @param track The track on which the genes are to be plotted, strating from the bottom to up as integers 1,2,...
  #' @export
  if(J.pos==1){
    pc<-15
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRList[[track]][2]+IRList[[track]][3]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
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
      if(tup[i,1]=="trnH"){col[i]<- "burlywood4"}
      else if(tup[i,1]=="trnN"){col[i]<- "violet"}
      else if(tup[i,1]=="psbA"){col[i]<- "purple2"}
      else if(tup[i,1]=="rps19"){col[i]<- "firebrick1"}
      else if(tup[i,1]=="ycf1"){col[i]<- "dodgerblue3"}
      else if(tup[i,1]=="ndhF"){col[i]<- "darkred"}
      else if(tup[i,1]=="rpl22"){col[i]<- "navy"}
      else if(tup[i,1]=="rpl2"){col[i]<- "forestgreen"}
      else if(tup[i,1]=="rpl23"){col[i]<- "palevioletred4"}
      else if(tup[i,1]=="ycf2"){col[i]<- "blue3"}
      else if(tup[i,1]=="rrn23"){col[i]<- "deeppink"}
      else if(tup[i,1]=="rrn4"){col[i]<- "deeppink3"}
      else if(tup[i,1]=="rrn5"){col[i]<- "deeppink4"}
      else if(tup[i,1]=="chlL"){col[i]<- "darkorange1"}
      else if(tup[i,1]=="chlN"){col[i]<- "darkorange3"}
      else if(tup[i,1]=="rps12"){col[i]<- "slategray1"}
      else if(tup[i,1]=="rps7"){col[i]<- "slategray3"}
      else if(tup[i,1]=="rps3"){col[i]<- "slategrey"}
      else if(tup[i,1]=="16rrn"){col[i]<- "darkkhaki"}
      else if(tup[i,1]=="trnV"){col[i]<- "slateblue4"}
      else if(tup[i,1]=="trnM"){col[i]<- "steelblue4"}
      else if(tup[i,1]=="trnL"){col[i]<- "steelblue"}
      else if(tup[i,1]=="ccsA"){col[i]<- "goldenrod1"}
      else if(tup[i,1]=="rpl32"){col[i]<- "darkolivegreen4"}
      else if(tup[i,1]=="ndhB"){col[i]<- "lightseagreen"}
      else if(tup[i,1]=="rps15"){col[i]<- "lightsalmon1"}
      else if(tup[i,1]=="ndhH"){col[i]<- "navajowhite"}
      else if(tup[i,1]=="ndhA"){col[i]<- "palevioletred2"}
      else {col[i]<- "black"}
    }
    return(col)
  }
  if (J.pos==4){
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
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
      for (j in seq(1.1, 1.7, 0.05)){
        segments(max(Pcord[i,1], pc-10), track*5-j+5, min(Pcord[i,2], pc+10), track*5-j+5, lwd=1, col=paste(gcol(tup)[i]))
      }
    }
  }
}


GN.plotter<- function(Radius, J.pos, track){#GeneName plotter
  #' Gene Name plotter
  #'
  #' Plotting the gene names on a given gene which is already plotted on the tracks of the IR plot
  if(J.pos==1){
    pc<-15
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRList[[track]][2]+IRList[[track]][3]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
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
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (min(x1, x2) >= 104){
        text(103.5, track*5+1.05+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
      else if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5+0.35+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
        text(max(x1, x2)+1, track*5+0.3+5, paste(Rcord[i,1]-Rcord[i,2], "bp", " "), cex=numcex, col=txtincol, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+0.35+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+1.15+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5-1.40+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
        text(max(x1, x2)+1, track*5-1.45+5, paste(Rcord[i,2]-Rcord[i,1], "bp", " "), cex=numcex, col=txtincol, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-1.40+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-2.25+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
    }
  }
}


OJ.plotter<- function(Radius, J.pos, track){#GeneName plotter
  #' On Junction plotter
  #'
  #' Plotting the fine tuned narrow lines showing the limits of the genes which are passing through the junction sites with their bp
  if(J.pos==1){
    pc<-15
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRList[[track]][2]+IRList[[track]][3]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
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
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
  }
  for (i in 1:n){
    if(Rcord[i,1]>Rcord[i,2]){
      x1<-min(Pcord[i,1], pc+10)
      x2<-max(Pcord[i,2], pc-10)
      if (Rcord[i,1] > J & Rcord[i,2] <J ){
        Arrows(min(x1, x2), track*5+1.1+5, pc-0.15, track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5+1.1+5, min(x1, x2), track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5+1.1+5, max(x1, x2), track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(max(x1, x2), track*5+1.1+5, pc+0.15, track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5+1.4+5, paste(Rcord[i, 1]-J+1, "bp", sep=" "), cex=0.4, pos=4)#up-right
        text(pc+1.4, track*5+1.4+5, paste(J-Rcord[i, 2], "bp", sep=" "), cex=0.4, pos=2)#up-left
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (Rcord[i,2] > J & Rcord[i,1] <J ){
        Arrows(max(x1, x2), track*5-2.1+5, pc+0.15, track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5-2.1+5, max(x1, x2), track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5-2.1+5, min(x1, x2), track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(min(x1, x2), track*5-2.1+5, pc-0.15, track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5-2.6+5, paste(Rcord[i, 2]-J, "bp", sep=" "), cex=0.4, pos=4)#low-right
        text(pc+1.4, track*5-2.6+5, paste(J-Rcord[i, 1], "bp", sep=" "), cex=0.4, pos=2)#low-left
      }
    }
  }
}


JD.plotter<- function(Radius, J.pos, track){#GeneName plotter
  #' Junction Distance plotter
  #'
  #' plotting the narrow lines of the distance of the genes for the junction sites which are not passing through any gene and their bp
  if(J.pos==1){
    pc<-15
    J<- IRList[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRList[[track]][1]+IRList[[track]][3]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRList[[track]][2]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRList[[track]][2]+IRList[[track]][3]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
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
    ind<-Pcord < 0
    Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
    Rcord[ind]<- Rcord[ind] + IRList[[track]][4]
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
    #with this setting the position of the zero (genes tangant to the junction site) will not be plotted. If interested either put "=" for the middle condition of the left or right binary operator or better develop zero only handling if function
    if     (Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top,right, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5+1.3+5), curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 -2 , track*5+1.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, right, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5-2.3+5), curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5-2.3-0.2+5, paste(Rcord[row.cor, 1]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 -3 , track*5-2.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 -4.5 , track*5-2.3-0.2+5, paste(Rcord[row.cor, 1]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5-2.3+5), curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4, track*5-2.3-0.2+5, paste(J-Rcord[row.cor, 2], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 + 3 , track*5-2.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4.5 , track*5-2.3-0.2+5, paste(J-Rcord[row.cor, 2], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5+1.3+5), curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, paste(J-Rcord[row.cor, 1], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 + 2 , track*5+1.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4 , track*5+1.3+0.2+5, paste(J-Rcord[row.cor, 1], "bp", sep=" "), cex=0.4)
    }
  }
}


Max.Radius<-function(J.pos, l, genelist, IRlist){
  if(J.pos==1){
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
    t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
    while(nrow(t)==0){
      Radius<- Radius+(0.2)*Radius#0.2
      t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
    }
    R[track]<-Radius
  }
  # If the difference is not too big, choose the same radius for all
  if(max(R)-min(R) < 100){ # what difference TODO
    Radius <- round(max(R)+1)
    for(track in 1:l){
      R[track] <- Radius
    }
  }
  return(R)
}


trnfixer<- function(Genelist){
  for (i in 1:length(Genelist[,1])){

    if (Genelist[i,1] %in% c("trnK-UUU", "tRNA-Lys")){
      Genelist[i,1]<- "trnK"
    }

    if (Genelist[i,1] %in% c("trnN-UUA", "tRNA-Asn")){
      Genelist[i,1]<- "trnN"
    }

    if (Genelist[i,1] %in% c("trnT-UGA", "tRNA-Thr")){
      Genelist[i,1]<- "trnT"
    }

    if (Genelist[i,1] %in% c("trnR-UCU", "tRNA-Arg")){
      Genelist[i,1]<- "trnR"
    }

    if (Genelist[i,1] %in% c("trnM-UAC", "tRNA-Met")){
      Genelist[i,1]<- "trnM"
    }

    if (Genelist[i,1] %in% c("trnI-UAA", "tRNA-Ile")){
      Genelist[i,1]<- "trnI"
    }

    if (Genelist[i,1] %in% c("trnQ-GUU", "tRNA-Gln")){
      Genelist[i,1]<- "trnQ"
    }

    if (Genelist[i,1] %in% c("trnH-GUA", "tRNA-His", "trnH-GUG")){
      Genelist[i,1]<- "trnH"
    }

    if (Genelist[i,1] %in% c("trnP-GGG", "tRNA-Pro")){
      Genelist[i,1]<- "trnP"
    }

    if (Genelist[i,1] %in% c("trnE-CUC", "tRNA-Glu")){
      Genelist[i,1]<- "trnE"
    }

    if (Genelist[i,1] %in% c("trnD-CUA", "tRNA-Asp")){
      Genelist[i,1]<- "trnD"
    }

    if (Genelist[i,1] %in% c("trnA-CGC", "tRNA-Ala")){
      Genelist[i,1]<- "trnA"
    }

    if (Genelist[i,1] %in% c("trnG-CCC", "tRNA-Gly")){
      Genelist[i,1]<- "trnG"
    }

    if (Genelist[i,1] %in% c("trnV-CAC", "tRNA-Val")){
      Genelist[i,1]<- "trnV"
    }

    if (Genelist[i,1] %in% c("trnY-AUA", "tRNA-Tyr")){
      Genelist[i,1]<- "trnY"
    }

    if (Genelist[i,1] %in% c("trnS-AGA", "tRNA-Ser")){
      Genelist[i,1]<- "trnS"
    }

    if (Genelist[i,1] %in% c("trnW-ACC", "tRNA-Trp")){
      Genelist[i,1]<- "trnW"
    }

    if (Genelist[i,1] %in% c("trnC-ACA", "tRNA-Cys")){
      Genelist[i,1]<- "trnC"
    }

    if (Genelist[i,1] %in% c("trnL-AAU", "tRNA-Leu")){
      Genelist[i,1]<- "trnL"
    }

    if (Genelist[i,1] %in% c("trnF-AAA", "tRNA-Phe")){
      Genelist[i,1]<- "trnF"
    }
  }
  return(Genelist)
}


trnCut<- function(GL){
  name<- GL[,1]
  for (i in 1:length(name)){
    if (length(strsplit(as.character(name[i]), split="")[[1]])>5){
      name[i]<-paste(strsplit(as.character(name[i]), split="")[[1]][1:4], collapse = "")
    }
  }
  GL[,1]<- as.character(name)
  return(GL)
}#second version for dogma input and fix and output


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


LSC<- function(IRinfo){#a vecotr as a list
  c<-as.vector(IRinfo)
  IR<- c[3]
  SSC<- c[2]-(c[1]+IR)
  return(c[4]-(SSC+2*IR))
}


SSC<- function(IRinfo){
  c<-as.vector(IRinfo)
  return(c[2]-(c[1]+c[3]))
}


## TODO comment gbfiles named list (gbfile, local.file)
IRs<- function(gbfiles, Sfiles, progress){
  l<<- length(gbfiles)#STOP with the proper length condition to be added
  
  FasList<<-  list()
  IRList<<-   list()
  GeneList<<- list()
  spnames<<-  list()
  nuclw<<- numeric(l)
  for (i in 1:l){
    gb.info <- gbfiles[[i]]
    
    if(gb.info$local.file) {
      gb.data <- read.gb(gb.info$gbfile)
    } else {
      gb.data <- fetch.gb(gb.info$gbfile)
    }
    GeneList[[i]] <<-gene.cordinates(gb.data)
    
    IRList[[i]] <<-  get_ir(gb.info$gbfile, gb.info$local.file)
    
    rev <- Sfiles[[i]]
    #form one IRListDimp
    
    if (rev){
      GeneList[[i]]<<- SSCrev(GeneList[[i]], (IRList[[i]][1]+IRList[[i]][3]), IRList[[i]][2])
    }
    spnames[[i]]<<- sp.name(gb.data)
    
    FasList[[i]]<<- rdnFixer(gb.data)
    nuclw[i]<<- 100/IRList[[i]][4]
    
    GeneList[[i]]<<- trnfixer(GeneList[[i]])
    GeneList[[i]]<<- trnCut(GeneList[[i]])
    
    progress$set(value = i, 
                 detail = paste0('This may take a while... ', i, '/', l))
  }
}

calculate.max.radius <- function(l, genelist, IRlist){
  radius <- list()
  for(i in 1:l){
    r <- Max.Radius(i, l, genelist, IRlist)
    print(r)
    print(radius)
    radius <-append(radius, r)
  }
  return(radius)
}

plot.data.aux <- function(radius, no.junction, l){
  for (i in 1:l){JG.plotter(radius[i], no.junction, i)}
  for (i in 1:l){GN.plotter(radius[i], no.junction, i)}
  for (i in 1:l){OJ.plotter(radius[i], no.junction, i)}
  for (i in 1:l){JD.plotter(radius[i], no.junction, i)}
}

# TODO comment decir que es para graficos
IRs2<- function(file="IR_out"){
  names(FasList)<-  unlist(spnames)
  names(IRList)<-   unlist(spnames)
  names(GeneList)<- unlist(spnames)
  
  if(anyNA(IRList, recursive= TRUE) == FALSE) {
    width<-8.3
    height<- l
    # dev.new(width, height)
    m<-matrix(rep(c(15, 25, 30, 25, 10), l+2), 5, l+2)
    m[,l+2]<-rep(NA,5)
    m[,1]<-rep(0,5)
    # if (file_type  == "pdf") {
    #   grDevices::pdf(file, width=8.3, height=(l+2)*8.3/12)#I changed the track to 12 from 11#alternative togglier for height can be l*.83 /// or l ///(l*7.47/10+0.83) ///(l+1)*8.3/11
    # } else {
    #   do.call(file_type, args=list(filename=paste('IR', file_type, sep = '.'),
    #                                width = 8.3, height = (l+2)*8.3/12, 
    #                                units = "in", res = 300))
    # }
    par(mai=c(0.5, 0.6+max(chr.count(unlist(spnames)))/20, 0.7, 0.5))# c(bottom, left, top, right)
    barplot(m, horiz=T, lwd=2, cex.names=0.46, space = 4, border = T, axes = F, 
            col=c("lightblue", "orange1", "lightgreen", "orange1", "lightblue"))
    title(main="Inverted Repeats")
    segments(15, 5*l+8, 15, 7, lty=1, lwd=1, col = "darkgrey")
    segments(40, 5*l+8, 40, 7, lty=1, lwd=1, col = "darkgrey")
    segments(70, 5*l+8, 70, 7, lty=1, lwd=1, col = "darkgrey")
    segments(95, 5*l+8, 95, 7, lty=1, lwd=1, col = "darkgrey")
    text(15, 5*l+9, "JLB", font=4, cex=0.7)
    text(40, 5*l+9, "JSB", font=4, cex=0.7)
    text(70, 5*l+9, "JSA", font=4, cex=0.7)
    text(95, 5*l+9, "JLA", font=4, cex=0.7)
    
    intextcol<- "white"
    innocol<- "blue"
    for (i in seq(9.5, 5*l+9, 5)){
      text(101, i, "LSC", cex=0.57, font=2, col=intextcol)
      text(2.5, i, "LSC", cex=0.57, font=2, col=intextcol)
      text(18, i, "IRb", cex=0.57, font=2, col=intextcol)
      text(43, i, "SSC", cex=0.57, font=2, col=intextcol)
      text(73, i, "IRa", cex=0.57, font=2, col=intextcol)
      points(27.5, i, cex=2.3, pch=18, col="white")
      points(55, i, cex=2.3, pch=18, col="white")
      points(82.5, i, cex=2.3, pch=18, col="white")
      text(27.5, i, "//", cex=0.95, font=1)
      text(55, i, "//", cex=0.95, font=1)
      text(82.5, i, "//", cex=0.95, font=1)
      count<-which(seq(9.5, 5*l+9, 5)==i)
      LSCp<-paste(prettyNum(LSC(IRList[[count]]), big.mark = ","), "bp", "")
      SSCp<-paste(prettyNum(SSC(IRList[[count]]), big.mark = ","), "bp", "")
      IRp<- paste(prettyNum(IRList[[count]][3], big.mark = ","), "bp", "")
      text(10.5, i-0.09, LSCp, font=4, cex=0.52, col=innocol)
      text(64.5, i-0.09, SSCp, font=4, cex=0.52, col=innocol)
      text(35.5, i-0.09, IRp, font=4, cex=0.52, col=innocol)
      text(90, i-0.09, IRp, font=4, cex=0.52, col=innocol)
      axis(2, i-1.5, labels = paste(prettyNum(IRList[[count]][4], big.mark=","), "bp", ""), font = 4, cex.axis=0.7, las=2, tick=F)
    }
    axis(2, seq(9.5, 5*l+9, 5), labels = spnames, font = 4, cex.axis=0.7, las=2, tick=F)
    points(0, 4.5, cex=2.3, pch=18, col="white")
    

    I<-   Max.Radius(1, l, genelist = GeneList, IRlist = IRList)
    II<-  Max.Radius(2, l, genelist = GeneList, IRlist = IRList)
    III<- Max.Radius(3, l, genelist = GeneList, IRlist = IRList)
    IV<-  Max.Radius(4, l, genelist = GeneList, IRlist = IRList)

    plot.data.aux(I, 1, l)
    plot.data.aux(II, 2, l)
    plot.data.aux(III, 3, l)
    plot.data.aux(IV, 4, l)

    # radius <- calculate.max.radius(l, GeneList, IRList)
    # print(radius)
    # 
    # for(i in 1:l){
    #   plot.data.aux(radius[[i]], 1, l)
    # }
    
    if (length(unique(I)) != 1 || length(unique(II)) != 1
        || length(unique(III)) != 1 || length(unique(IV)) != 1){
      # text <- "\nWARNING: The plots are not in scale!\n"
      # textplot(text, cex = 0.9)
      # return(TRUE)
    } else {
      # return(FALSE)
    }

    
    # while (!is.null(dev.list()))  dev.off()    
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
    
    width<-8.3
    height<- l
    dev.new(width, height)
    m<-matrix(rep(c(15, 25, 30, 25, 10), l+2), 5, l+2)
    m[,l+2]<-rep(NA,5)
    m[,1]<-rep(0,5)
    grDevices::pdf(file, width=8.3, height=(l+2)*8.3/12)
    
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, text, cex = 1.6, col = "black")
    par(mar = c(5, 4, 4, 2) + 0.1)
    dev.off()
    # return(FALSE)
  }
}##For the GB files


###This section contains the functions related to the dogma input###
rdnFixerD<- function(FasGenome){
  seq<- FasGenome
  seq[which(seq=="u")]<-sample(c("t"), length(which(seq=="u")), TRUE)
  seq[which(seq=="r")]<-sample(c("a", "g"), length(which(seq=="r")), TRUE)
  seq[which(seq=="y")]<-sample(c("c", "t"), length(which(seq=="y")), TRUE)
  seq[which(seq=="s")]<-sample(c("c", "g"), length(which(seq=="s")), TRUE)
  seq[which(seq=="w")]<-sample(c("a", "t"), length(which(seq=="w")), TRUE)
  seq[which(seq=="k")]<-sample(c("g", "t"), length(which(seq=="k")), TRUE)
  seq[which(seq=="m")]<-sample(c("c", "a"), length(which(seq=="m")), TRUE)
  seq[which(seq=="b")]<-sample(c("c", "g", "t"), length(which(seq=="b")), TRUE)
  seq[which(seq=="d")]<-sample(c("a", "g", "t"), length(which(seq=="d")), TRUE)
  seq[which(seq=="h")]<-sample(c("c", "a", "t"), length(which(seq=="h")), TRUE)
  seq[which(seq=="v")]<-sample(c("c", "a", "g"), length(which(seq=="v")), TRUE)
  seq[which(seq=="n")]<-sample(c("c", "g", "t", "a"), length(which(seq=="n")), TRUE)
  seq[which(seq=="-")]<-sample(c("c", "g", "t", "a"), length(which(seq=="-")), TRUE)
  return(seq)
}##for the Dogma and fasta

# TODO borrar
# sp.nameD<- function(gb){
#   paste(strsplit(gb[2], " ")[[1]][3], strsplit(gb[2], " ")[[1]][4], sep=" ")
# }##for the dogma input, read the name from the header of fasta files

# TODO borrar
# trnfixerD<- function(Genelist){
#   Genelist<- Genelist[,c(3,2,1,4)]
#   for (i in 1:length(Genelist[,1])){
# 
#     if (Genelist[i,1] %in% c("trnK-UUU", "tRNA-Lys")){
#       Genelist[i,1]<- "trnK"
#     }
# 
#     if (Genelist[i,1] %in% c("trnN-UUA", "tRNA-Asn")){
#       Genelist[i,1]<- "trnN"
#     }
# 
#     if (Genelist[i,1] %in% c("trnT-UGA", "tRNA-Thr")){
#       Genelist[i,1]<- "trnT"
#     }
# 
#     if (Genelist[i,1] %in% c("trnR-UCU", "tRNA-Arg")){
#       Genelist[i,1]<- "trnR"
#     }
# 
#     if (Genelist[i,1] %in% c("trnM-UAC", "tRNA-Met")){
#       Genelist[i,1]<- "trnM"
#     }
# 
#     if (Genelist[i,1] %in% c("trnI-UAA", "tRNA-Ile")){
#       Genelist[i,1]<- "trnI"
#     }
# 
#     if (Genelist[i,1] %in% c("trnQ-GUU", "tRNA-Gln")){
#       Genelist[i,1]<- "trnQ"
#     }
# 
#     if (Genelist[i,1] %in% c("trnH-GUA", "tRNA-His")){
#       Genelist[i,1]<- "trnH"
#     }
# 
#     if (Genelist[i,1] %in% c("trnP-GGG", "tRNA-Pro")){
#       Genelist[i,1]<- "trnP"
#     }
# 
#     if (Genelist[i,1] %in% c("trnE-CUC", "tRNA-Glu")){
#       Genelist[i,1]<- "trnE"
#     }
# 
#     if (Genelist[i,1] %in% c("trnD-CUA", "tRNA-Asp")){
#       Genelist[i,1]<- "trnD"
#     }
# 
#     if (Genelist[i,1] %in% c("trnA-CGC", "tRNA-Ala")){
#       Genelist[i,1]<- "trnA"
#     }
# 
#     if (Genelist[i,1] %in% c("trnG-CCC", "tRNA-Gly")){
#       Genelist[i,1]<- "trnG"
#     }
# 
#     if (Genelist[i,1] %in% c("trnV-CAC", "tRNA-Val")){
#       Genelist[i,1]<- "trnV"
#     }
# 
#     if (Genelist[i,1] %in% c("trnY-AUA", "tRNA-Tyr")){
#       Genelist[i,1]<- "trnY"
#     }
# 
#     if (Genelist[i,1] %in% c("trnS-AGA", "tRNA-Ser")){
#       Genelist[i,1]<- "trnS"
#     }
# 
#     if (Genelist[i,1] %in% c("trnW-ACC", "tRNA-Trp")){
#       Genelist[i,1]<- "trnW"
#     }
# 
#     if (Genelist[i,1] %in% c("trnC-ACA", "tRNA-Cys")){
#       Genelist[i,1]<- "trnC"
#     }
# 
#     if (Genelist[i,1] %in% c("trnL-AAU", "tRNA-Leu")){
#       Genelist[i,1]<- "trnL"
#     }
# 
#     if (Genelist[i,1] %in% c("trnF-AAA", "tRNA-Phe")){
#       Genelist[i,1]<- "trnF"
#     }
#   }
#   return(Genelist)
# }##For the dogma and fasta


trnDogma<- function(dogma){##takes in the dogma input and return the same output but with the gene names fixed
  name<- as.character(dogma[,3])
  for (i in 1:length(name)){
    if (length(strsplit(as.character(name[i]), split="")[[1]])>5){
      name[i]<-paste(strsplit(as.character(name[i]), split="")[[1]][1:4], collapse = "")
    }
  }
  dogma[,3]<- as.character(name)
  return(dogma)
}#second version for dogma input and fix and output


GnlBuilder<- function(dogma){##takes the dogma input and return the type of the list suited for the IR and other functions
  l<- length(dogma[,1])
  list<- cbind(as.character(rep(0, l)),as.character(rep(0, l)),as.character(rep(0, l)),as.character(rep(0, l)),as.character(rep(0, l)))
  for (i in 1:l){
    if (as.character(dogma[i,4])=="-"){
      p<-dogma[i, 2]
      dogma[i,2]<- dogma[i,1]
      dogma[i,1]<- p
    }
    list[i,1]<-as.character(dogma[i,3])
    list[i,2]<-as.character(dogma[i,1])
    list[i,3]<-as.character(dogma[i,2])
  }
  return(list)
}
###Suppose that you are given two vectors of IRinp and one value of the genome size.
###IRinp going to be the vector with input cordinates of the JLB, JSB, JSA, JLA, and genome size
#for its first, second, third, fourth, and fifth values respectively
JunctRadiusFinderDinp<- function(gene.cordinates, IRinp, J.pos, radius, silence=TRUE){#function to find the genes near the junction cordinate based on the output of the gene.cordinate function
  if(J.pos==1){
    J<- IRinp[1]
  }
  else if(J.pos==2){
    J<- IRinp[2]
  }
  else if(J.pos==3){
    J<- IRinp[3]
  }
  else if(J.pos==4){
    J<- IRinp[4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  g<-gene.cordinates
  l<- ncol(g)
  r<- radius
  gs<- IRinp[5]
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
  return(t)
}

# TODO borrar?
# JunctNoFinderDinp<- function(gene.cordinates, IRinp, J.pos ,Number){#the same function like it precusor but returns the number of genes found near the junction site
#   if(J.pos==1){
#     J<- IRinp[1]
#   }
#   else if(J.pos==2){
#     J<- IRinp[2]
#   }
#   else if(J.pos==3){
#     J<- IRinp[3]
#   }
#   else if(J.pos==4){
#     J<- IRinp[4]
#   }
#   else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
#   g<-gene.cordinates
#   l<- ncol(g)
#   n<-Number
#   counter=0
#   count=0
#   while (nrow(JunctRadiusFinderDinp(g, IRinp, J.pos ,count)) <= n-1){
#     count<-count+1
#   }
#   t<-JunctRadiusFinderDinp(g, IRinp, J.pos ,count)
#   t[is.na(t)]<- "0"
#   tup<- matrix(0, n, 3)
#   for (i in 1:n){
#     dist<-abs(as.numeric(t[i,][2:5])-J)
#     if (which(dist==min(dist))>3){
#       tup[i,]<-t[i, c(1,4,5)]
#     }
#     else {
#       tup[i,]<-t[i, c(1,2,3)]
#     }
#   }
#   return(tup)
# }

# TODO borrar?
# thFinderDinp<- function(gene.cordinates, IRinp, J.pos, n){#finding the nth closest gene to the junctions site
#   g<-gene.cordinates
#   l<- ncol(g)
#   if(n==1){
#     return(JunctNoFinderDinp(g, IRinp, J.pos, n))
#   }
#   else{
#     previous<-JunctNoFinderDinp(g, IRinp, J.pos, n-1)
#     current<- JunctNoFinderDinp(g, IRinp, J.pos, n  )
#     return(current[!current[,2] %in% previous[,2],])
#   }
# }

# TODO borrar
# RadiusNoFinderDinp<- function(gene.cordinates, IRinp, J.pos, Number){
#   g<-gene.cordinates
#   l<- ncol(g)
#   n<-Number
#   counter=0
#   count=0
#   while (nrow(JunctRadiusFinderDinp(g, IRinp, J.pos, count)) <= n-1){
#     count<-count+1
#   }
#   return(count)
# }


JG.plotterDinp<- function(Radius, J.pos, track){#JunctionGene plotter given an global IRListDinp
  #' Junction site gene plotter
  #'
  #' Plotting the genes in the vicinity of the junction site of the chloroplast.
  #' @param gene.cordinate The output of the gene.cordinate function with the name of the genes and their cordinates on genome
  #' @param junction.cordinate The cordinate of the junction site for which the genes are to be plotted
  #' @param n The number of the genes that is needed to be plotted
  #' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, and 4 respectively
  #' @param margin The margin value to be added to the interval as the proportion of the minimum lenght of the nth gene in vicinity
  #' @param track The track on which the genes are to be plotted, strating from the bottom to up as integers 1,2,...
  #' @export
  if(J.pos==1){
    pc<-15
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRListDinp[[track]][4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
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
      if(tup[i,1]=="trnH"){col[i]<- "burlywood4"}
      else if(tup[i,1]=="trnN"){col[i]<- "violet"}
      else if(tup[i,1]=="psbA"){col[i]<- "purple2"}
      else if(tup[i,1]=="rps19"){col[i]<- "firebrick1"}
      else if(tup[i,1]=="ycf1"){col[i]<- "dodgerblue3"}
      else if(tup[i,1]=="ndhF"){col[i]<- "darkred"}
      else if(tup[i,1]=="rpl22"){col[i]<- "navy"}
      else if(tup[i,1]=="rpl2"){col[i]<- "forestgreen"}
      else if(tup[i,1]=="rpl23"){col[i]<- "blue3"}
      else {col[i]<- "black"}
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
      for (j in seq(1.1, 1.7, 0.05)){
        segments(max(Pcord[i,1], pc-10), track*5-j+5, min(Pcord[i,2], pc+10), track*5-j+5, lwd=1, col=paste(gcol(tup)[i]))
      }
    }
  }
}


GN.plotterDinp<- function(Radius, J.pos, track){#GeneName plotter
  #' Gene Name plotter
  #'
  #' Plotting the gene names on a given gene which is already plotted on the tracks of the IR plot
  if(J.pos==1){
    pc<-15
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRListDinp[[track]][4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
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
        text(103.5, track*5+1.05+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
      else if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5+0.35+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
        text(max(x1, x2)+1, track*5+0.3+5, paste(Rcord[i,1]-Rcord[i,2], "bp", " "), cex=numcex, col=txtincol, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+0.35+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5+1.15+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (abs(x1-x2) >= 9.7){
        text(min(x1, x2)+1.8, track*5-1.40+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
        text(max(x1, x2)+1, track*5-1.45+5, paste(Rcord[i,2]-Rcord[i,1], "bp", " "), cex=numcex, col=txtincol, font=numfont, pos=2)
      }
      else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-1.40+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
      }
      else if (abs(x1-x2) < 3){
        text(((min(x1, x2)+max(x1, x2))/2), track*5-2.25+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
      }
    }
  }
}


OJ.plotterDinp<- function(Radius, J.pos, track){#GeneName plotter
  #' On Junction plotter
  #'
  #' Plotting the fine tuned narrow lines showing the limits of the genes which are passing through the junction sites with their bp
  if(J.pos==1){
    pc<-15
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRListDinp[[track]][4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
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
        Arrows(min(x1, x2), track*5+1.1+5, pc-0.15, track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5+1.1+5, min(x1, x2), track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5+1.1+5, max(x1, x2), track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(max(x1, x2), track*5+1.1+5, pc+0.15, track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5+1.4+5, paste(Rcord[i, 1]-J+1, "bp", sep=" "), cex=0.4, pos=4)#up-right
        text(pc+1.4, track*5+1.4+5, paste(J-Rcord[i, 2], "bp", sep=" "), cex=0.4, pos=2)#up-left
      }
    }
    else {
      x1<-max(Pcord[i,1], pc-10)
      x2<-min(Pcord[i,2], pc+10)
      if (Rcord[i,2] > J & Rcord[i,1] <J ){
        Arrows(max(x1, x2), track*5-2.1+5, pc+0.15, track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc+0.15, track*5-2.1+5, max(x1, x2), track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(pc-0.15, track*5-2.1+5, min(x1, x2), track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        Arrows(min(x1, x2), track*5-2.1+5, pc-0.15, track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
        text(pc-1.4, track*5-2.6+5, paste(Rcord[i, 2]-J, "bp", sep=" "), cex=0.4, pos=4)#low-right
        text(pc+1.4, track*5-2.6+5, paste(J-Rcord[i, 1], "bp", sep=" "), cex=0.4, pos=2)#low-left
      }
    }
  }
}


JD.plotterDinp<- function(Radius, J.pos, track){#GeneName plotter
  #' Junction Distance plotter
  #'
  #' plotting the narrow lines of the distance of the genes for the junction sites which are not passing through any gene and their bp
  if(J.pos==1){
    pc<-15
    J<- IRListDinp[[track]][1]
  }
  else if(J.pos==2){
    pc<-40
    J<- IRListDinp[[track]][2]
  }
  else if(J.pos==3){
    pc<-70
    J<- IRListDinp[[track]][3]
  }
  else if(J.pos==4){
    pc<-95
    J<- IRListDinp[[track]][4]
  }
  else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
  txtincol<- "white"
  txtoutcol="black"
  txtfont<- 4
  numfont<- 1
  numcex<- 0.44
  txtcex<- 0.46
  t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
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
    #with this setting the position of the zero (genes tangant to the junction site) will not be plotted. If interested either put "=" for the middle condition of the left or right binary operator or better develop zero only handling if function
    if     (Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top,right, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5+1.3+5), curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 -2 , track*5+1.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, right, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5-2.3+5), curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5-2.3-0.2+5, paste(Rcord[row.cor, 1]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, right, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 -3 , track*5-2.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 -4.5 , track*5-2.3-0.2+5, paste(Rcord[row.cor, 1]-J, "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5-2.3+5), curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4, track*5-2.3-0.2+5, paste(J-Rcord[row.cor, 2], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 + 3 , track*5-2.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4.5 , track*5-2.3-0.2+5, paste(J-Rcord[row.cor, 2], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top, left, big
      curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5+1.3+5), curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, paste(J-Rcord[row.cor, 1], "bp", sep=" "), cex=0.4)
    }
    else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, left, small
      arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 + 2 , track*5+1.3+5, angle = 15, length = 0.05, lwd=0.6)
      text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4 , track*5+1.3+0.2+5, paste(J-Rcord[row.cor, 1], "bp", sep=" "), cex=0.4)
    }
  }
}


Max.RadiusDinp<-function(J.pos, l, genelist, IRlistDinp){
  if(J.pos==1){
    Radius<-680#680#100
  }
  if(J.pos==2){
    Radius<-100#100#100
  }
  if(J.pos==3){
    Radius<-1000#1800#100
  }
  if(J.pos==4){
    Radius<-1000#1000#100
  }
  R<- numeric(l)
  for (track in 1:l){
    t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
    while(nrow(t)==0){
      Radius<- Radius+(0.2)*Radius
      t<- JunctRadiusFinderDinp(GeneList[[track]], IRListDinp[[track]], J.pos , Radius)
    }
    R[track]<-Radius
  }
  return(round(max(R)+1))
}


LSCDinp<- function(IRinp){#a vecotr as a list
  c<-as.vector(IRinp)
  return(c[5]-((c[2]-c[1])+(c[4]-c[3])+(c[3]-c[2])))
}


SSCDinp<- function(IRinp){
  c<-as.vector(IRinp)
  return(c[3]-c[2])
}

## CHANGE THIS ONE TODO IRinfo < get_ir
###everything is cleared for Dogma
IRsD<- function(dgfiles, fastafiles, irfiles, nfiles, file="IR_out"){
  l<<- length(dgfiles)#STOP with the proper length condition to be added
  FasList<<-  list()
  IRListDinp<<- list()
  GeneList<<- list()
  spnames<<-  list()
  nuclw<<- numeric(l)
  for (i in 1:l){
    dg<-dgfiles[[i]]
    FasList[[i]]<<- fastafiles[[i]]
    GeneList[[i]]<<-dgfiles[[i]]
    spnames[[i]]<<- nfiles[[i]]
    if (irfiles[[i]][1]==999){
      sss<- IRinfo(as.vector(FasList[[i]]))
      IRListDinp[[i]]<<-  c(sss[1],sss[1]+sss[3],sss[2],sss[2]+sss[3],sss[4])
    }
    else {IRListDinp[[i]]<<- irfiles[[i]]}
    nuclw[i]<<- 100/length(FasList[[1]])
  }
  for (i in 1:l){
    GeneList[[i]]<<- GeneList[[i]]
  }
}

IRsD2<- function(dgfiles, file="IR_out"){
  names(FasList)<- unlist(spnames)
  names(IRListDinp)<- unlist(spnames)
  names(GeneList)<- unlist(spnames)
  #copied and paste the following line to make sure produce the graph
  width<-8.3
  height<- l
  dev.new(width, height)
  m<-matrix(rep(c(15, 25, 30, 25, 10), l+2), 5, l+2)
  m[,l+2]<-rep(NA,5)
  m[,1]<-rep(0,5)
  grDevices::pdf(file, width=8.3, height=(l+2)*8.3/12)#I changed the track to 12 from 11#alternative togglier for height can be l*.83 /// or l ///(l*7.47/10+0.83) ///(l+1)*8.3/11
  par(mai=c(0.5, 0.6+max(chr.count(unlist(spnames)))/20, 0.7, 0.5))# c(bottom, left, top, right)
  barplot(m, horiz=T, lwd=1,cex.names=0.46, space = 4, border = T, axes = F, col=c("lightblue", "orange1", "lightgreen", "orange1", "lightblue"))
  title(main="Inverted Repeats")
  segments(15, 5*l+8, 15, 7, lty=1, lwd=1, col = "darkgrey")
  segments(40, 5*l+8, 40, 7, lty=1, lwd=1, col = "darkgrey")
  segments(70, 5*l+8, 70, 7, lty=1, lwd=1, col = "darkgrey")
  segments(95, 5*l+8, 95, 7, lty=1, lwd=1, col = "darkgrey")
  text(15, 5*l+9, "JLB", font=4, cex=0.7)
  text(40, 5*l+9, "JSB", font=4, cex=0.7)
  text(70, 5*l+9, "JSA", font=4, cex=0.7)
  text(95, 5*l+9, "JLA", font=4, cex=0.7)

  intextcol<- "white"
  innocol<- "blue"
  for (i in seq(9.5, 5*l+9, 5)){
    text(101, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(2.5, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(18, i, "IRb", cex=0.57, font=2, col=intextcol)
    text(43, i, "SSC", cex=0.57, font=2, col=intextcol)
    text(73, i, "IRa", cex=0.57, font=2, col=intextcol)
    points(27.5, i, cex=2.3, pch=18, col="white")
    points(55, i, cex=2.3, pch=18, col="white")
    points(82.5, i, cex=2.3, pch=18, col="white")
    text(27.5, i, "//", cex=0.95, font=1)
    text(55, i, "//", cex=0.95, font=1)
    text(82.5, i, "//", cex=0.95, font=1)
    count<-which(seq(9.5, 5*l+9, 5)==i)
    LSCp<-paste(prettyNum(LSCDinp(IRListDinp[[count]]), big.mark = ","), "bp", "")
    SSCp<-paste(prettyNum(SSCDinp(IRListDinp[[count]]), big.mark = ","), "bp", "")
    IRpb<- paste(prettyNum((IRListDinp[[count]][2]-IRListDinp[[count]][1]), big.mark = ","), "bp", "")
    IRpa<- paste(prettyNum((IRListDinp[[count]][4]-IRListDinp[[count]][3]), big.mark = ","), "bp", "")
    text(10.5, i-0.09, LSCp, font=4, cex=0.52, col=innocol)
    text(64.5, i-0.09, SSCp, font=4, cex=0.52, col=innocol)
    text(35.5, i-0.09, IRpb, font=4, cex=0.52, col=innocol)
    text(90, i-0.09, IRpa, font=4, cex=0.52, col=innocol)
    axis(2, i-1.5, labels = paste(prettyNum(IRListDinp[[count]][5], big.mark=","), "bp", ""), font = 4, cex.axis=0.7, las=2, tick=F)
  }
  axis(2, seq(9.5, 5*l+9, 5), labels = spnames, font = 4, cex.axis=0.7, las=2, tick=F)
  points(0, 4.5, cex=2.3, pch=18, col="white")

  I<-   Max.RadiusDinp(1, l, genelist = GeneList, IRlistDinp = IRListDinp)
  II<-  Max.RadiusDinp(2, l, genelist = GeneList, IRlistDinp = IRListDinp)
  III<- Max.RadiusDinp(3, l, genelist = GeneList, IRlistDinp = IRListDinp)
  IV<-  Max.RadiusDinp(4, l, genelist = GeneList, IRlistDinp = IRListDinp)


  for (i in 1:l){JG.plotterDinp(I, 1, i)}
  for (i in 1:l){GN.plotterDinp(I, 1, i)}
  for (i in 1:l){JG.plotterDinp(II, 2, i)}#default has been 100
  for (i in 1:l){GN.plotterDinp(II, 2, i)}#
  for (i in 1:l){JG.plotterDinp(III, 3, i)}
  for (i in 1:l){GN.plotterDinp(III, 3, i)}
  for (i in 1:l){JG.plotterDinp(IV, 4, i)}
  for (i in 1:l){GN.plotterDinp(IV, 4, i)}
  for (i in 1:l){OJ.plotterDinp(I, 1, i)}
  for (i in 1:l){OJ.plotterDinp(III, 3, i)}
  for (i in 1:l){OJ.plotterDinp(II, 2, i)}
  for (i in 1:l){OJ.plotterDinp(IV, 4, i)}
  for (i in 1:l){JD.plotterDinp(I, 1, i)}
  for (i in 1:l){JD.plotterDinp(III, 3, i)}
  for (i in 1:l){JD.plotterDinp(II, 2, i)}
  for (i in 1:l){JD.plotterDinp(IV, 4, i)}
  #copied here
  width<-8.3
  height<- l
  dev.new(width, height)
  m<-matrix(rep(c(15, 25, 30, 25, 10), l+2), 5, l+2)
  m[,l+2]<-rep(NA,5)
  m[,1]<-rep(0,5)
  grDevices::pdf(file, width=8.3, height=(l+2)*8.3/12)#I changed the track to 12 from 11#alternative togglier for height can be l*.83 /// or l ///(l*7.47/10+0.83) ///(l+1)*8.3/11
  par(mai=c(0.5, 0.6+max(chr.count(unlist(spnames)))/20, 0.7, 0.5))# c(bottom, left, top, right)
  barplot(m, horiz=T, lwd=1,cex.names=0.46, space = 4, border = T, axes = F, col=c("lightblue", "orange1", "lightgreen", "orange1", "lightblue"))
  title(main="Inverted Repeats")
  segments(15, 5*l+8, 15, 7, lty=1, lwd=1, col = "darkgrey")
  segments(40, 5*l+8, 40, 7, lty=1, lwd=1, col = "darkgrey")
  segments(70, 5*l+8, 70, 7, lty=1, lwd=1, col = "darkgrey")
  segments(95, 5*l+8, 95, 7, lty=1, lwd=1, col = "darkgrey")
  text(15, 5*l+9, "JLB", font=4, cex=0.7)
  text(40, 5*l+9, "JSB", font=4, cex=0.7)
  text(70, 5*l+9, "JSA", font=4, cex=0.7)
  text(95, 5*l+9, "JLA", font=4, cex=0.7)

  intextcol<- "white"
  innocol<- "blue"
  for (i in seq(9.5, 5*l+9, 5)){
    text(101, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(2.5, i, "LSC", cex=0.57, font=2, col=intextcol)
    text(18, i, "IRb", cex=0.57, font=2, col=intextcol)
    text(43, i, "SSC", cex=0.57, font=2, col=intextcol)
    text(73, i, "IRa", cex=0.57, font=2, col=intextcol)
    points(27.5, i, cex=2.3, pch=18, col="white")
    points(55, i, cex=2.3, pch=18, col="white")
    points(82.5, i, cex=2.3, pch=18, col="white")
    text(27.5, i, "//", cex=0.95, font=1)
    text(55, i, "//", cex=0.95, font=1)
    text(82.5, i, "//", cex=0.95, font=1)
    count<-which(seq(9.5, 5*l+9, 5)==i)
    LSCp<-paste(prettyNum(LSCDinp(IRListDinp[[count]]), big.mark = ","), "bp", "")
    SSCp<-paste(prettyNum(SSCDinp(IRListDinp[[count]]), big.mark = ","), "bp", "")
    IRpb<- paste(prettyNum((IRListDinp[[count]][2]-IRListDinp[[count]][1]), big.mark = ","), "bp", "")
    IRpa<- paste(prettyNum((IRListDinp[[count]][4]-IRListDinp[[count]][3]), big.mark = ","), "bp", "")
    text(10.5, i-0.09, LSCp, font=4, cex=0.52, col=innocol)
    text(64.5, i-0.09, SSCp, font=4, cex=0.52, col=innocol)
    text(35.5, i-0.09, IRpb, font=4, cex=0.52, col=innocol)
    text(90, i-0.09, IRpa, font=4, cex=0.52, col=innocol)
    axis(2, i-1.5, labels = paste(prettyNum(IRListDinp[[count]][5], big.mark=","), "bp", ""), font = 4, cex.axis=0.7, las=2, tick=F)
  }
  axis(2, seq(9.5, 5*l+9, 5), labels = spnames, font = 4, cex.axis=0.7, las=2, tick=F)
  points(0, 4.5, cex=2.3, pch=18, col="white")

  I<-   Max.RadiusDinp(1, l, genelist = GeneList, IRlistDinp = IRListDinp)
  II<-  Max.RadiusDinp(2, l, genelist = GeneList, IRlistDinp = IRListDinp)
  III<- Max.RadiusDinp(3, l, genelist = GeneList, IRlistDinp = IRListDinp)
  IV<-  Max.RadiusDinp(4, l, genelist = GeneList, IRlistDinp = IRListDinp)


  for (i in 1:l){JG.plotterDinp(I, 1, i)}
  for (i in 1:l){GN.plotterDinp(I, 1, i)}
  for (i in 1:l){JG.plotterDinp(II, 2, i)}#default has been 100
  for (i in 1:l){GN.plotterDinp(II, 2, i)}#
  for (i in 1:l){JG.plotterDinp(III, 3, i)}
  for (i in 1:l){GN.plotterDinp(III, 3, i)}
  for (i in 1:l){JG.plotterDinp(IV, 4, i)}
  for (i in 1:l){GN.plotterDinp(IV, 4, i)}
  for (i in 1:l){OJ.plotterDinp(I, 1, i)}
  for (i in 1:l){OJ.plotterDinp(III, 3, i)}
  for (i in 1:l){OJ.plotterDinp(II, 2, i)}
  for (i in 1:l){OJ.plotterDinp(IV, 4, i)}
  for (i in 1:l){JD.plotterDinp(I, 1, i)}
  for (i in 1:l){JD.plotterDinp(III, 3, i)}
  for (i in 1:l){JD.plotterDinp(II, 2, i)}
  for (i in 1:l){JD.plotterDinp(IV, 4, i)}
  dev.off()
  dev.off()
  dev.off()
  dev.off()
}##For the Fasta and Dogma Files
