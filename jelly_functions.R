library(vegan)
library(MASS)
library(ape)
library(reshape2)
library(ggplot2)
library(ggvenn)
library(RColorBrewer)

### Set up the color pallette
col.a <- colorRampPalette(c("peachpuff", "sandybrown", "tomato","tomato3","lightblue","skyblue","skyblue4", "darkslategray"))

##To group samples into the respective locations (Hisdalen/Onarøy, Farm/Control)
fact.indx <- function(x, data.eco=eco){
    list(HF=x[which(data.eco$Fish.farm=="Hisdalen" & data.eco$Sampling.Site=="Near")],
         HC=x[which(data.eco$Fish.farm=="Hisdalen" & data.eco$Sampling.Site!="Near")],
         OF=x[which(data.eco$Fish.farm != "Hisdalen" & data.eco$Sampling.Site=="Near")],
         OC=x[which(data.eco$Fish.farm!="Hisdalen" & data.eco$Sampling.Site!="Near")])
}

### Group distribitions using piecharts, input: a matrix with rows of taxa and columns of samples
### Parameters: if read.count is false (default), it will count the number of unique phyla within each phylum, if TRUE, it will sum up all the aligned reads within each phylum. 

molto.pies <- function(tax.group=mtz, fam.count=1, other.thresh=3, readcount=FALSE){
    col.a <- colorRampPalette(c("sandybrown", "tomato","tomato4","skyblue","skyblue4", "darkslategray"))
    par(mar=c(5,5,5,5))  
    if(readcount==FALSE){
        fam.no <- (rev(sort(table(ph.fam[ph.fam$phylum%in%tax.group,"phylum"])))/nrow(ph.fam[ph.fam$phylum%in%tax.group,]))*100
        fam.other <- sum(fam.no[fam.no<other.thresh])
        fam.h <- c(fam.no[fam.no>other.thresh], fam.other)
        names(fam.h)[length(fam.h)] <- paste0("Other(n=", sum(fam.no<other.thresh),")")
        pie(fam.h, labels=paste0(names(fam.h), "\n", round(fam.h,1),"%"),
            col=col.a(length(fam.h)), las=2, cex=1.5, radius=.7)
        text(-0.7,-1,paste("Total genera:", sum(ph.fam$phylum%in%tax.group)), font=2, cex=2)
    }
    else if(readcount==TRUE){ 
        ph.sum.p <- rev(sort((rowSums(ph.matrix[rownames(ph.matrix)%in%tax.group,])/sum(rowSums(ph.matrix[rownames(ph.matrix)%in%tax.group,])))*100))
        ph.other <- sum(ph.sum.p[ph.sum.p<other.thresh]) 
        ph.h <- c(ph.sum.p[ph.sum.p>other.thresh], "Other"=ph.other) 
        pie(ph.h, labels=paste0(names(ph.h), "\n", round(ph.h,1),"%"),
            col=col.a(length(ph.h)), las=2, cex=1.7, radius=0.6)
        text(-0.7,-1,paste("Total reads:", round(sum(rowSums(ph.matrix[rownames(ph.matrix)%in%tax.group,])),2)), font=2, cex=1.7)
    }
}

### Getting the taxonomic distributions in relative abundance barplot form, emphasizing the domingant groups. Input is a matrix with taxa for rows and samples for columns.
### If plot.brplot is false, it returns the distributions in numerical form. 

barplot.gr <-function(tax.type="fam", tax.matrix=pol.matrix, tax.group=mtz, phyl="", other.thresh=3, plot.brplt=TRUE){
    if(tax.type=="phylum"){
        ph.m.tax <- ph.m.f[rownames(ph.m.f)%in%tax.group,]
        gr.m.tax <- ph.m.tax}
    else  if(tax.type=="fam"){
            fam.m.tax <- fam.m.f[rownames(fam.matrix)%in%sort(ph.fam[ph.fam$phylum==phyl,"family"]),]
            gr.m.tax <- fam.m.tax
        }
    else  if(tax.type=="genus"){
            gen.m.tax <- gen.m.f[rownames(gen.matrix)%in%sort(ph.fam[ph.fam$phylum==phyl,"genus"]),]
            gr.m.tax <- gen.m.tax
        }
    else  if(tax.type=="morpho"){
            gr.m.tax <-tax.matrix
        }
    ##Make proportions
    gr.m.tax.p <- round(sapply(1:ncol(gr.m.tax), function(x){(gr.m.tax[,x]/sum(gr.m.tax[,x]))*100}),3)
    ##Find the most abundant taxa among all groups
    gr.sums <- rowSums(gr.m.tax.p)
    ##Setting up the "other" category (groups with low relative abundance)
    gr.other <- names(gr.sums[gr.sums<=other.thresh])
    gr.h <- gr.m.tax.p[rownames(gr.m.tax.p)%in%names(gr.sums[gr.sums>other.thresh]),]
    gr.h <- rbind(gr.h, "Other"=100-colSums(gr.h))
    if(plot.brplt==FALSE){
        gr.h}
    else{
            ##Plot the barplots
            par(mar=c(6,8,1,1))
            layout(matrix(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2), nrow = 2))
            plot(1:10, xlim=c(1,5), ylim=c(0,107), type="n", ylab="",
                 axes=F, xlab="", xaxs="i",yaxs="i", cex.main=1.5, cex.axis=2.4, cex.lab=2.4)
            box(lty=1, lwd=3)
            mtext("Relative abundance (%)", side = 2, line = 5, cex=2.2)
            axis(1, seq(1.5, length(factor.names)+0.5,1), labels=c("HF", "HC", "OF", "OC"),
                 cex.axis=3, tick=F, font.axis=1, las=1, line=2)
            axis(2, seq(0, 100, 20), labels=(seq(0,100,20)), las=2, cex.axis=2.4)
            abline(h=seq(20, 100, 20), col="grey", lwd=3)
            group.names<-sort(unique(unlist(lapply(factor.sum,function(x){names(x)}))))
            group.names <- c(group.names[group.names!="Other"],"Other")
            rowSums(gr.m.tax.p)
            sapply(1:ncol(gr.m.tax), function(x){
                col.br=col.a(nrow(gr.h))
                rect(x+0.1, 0, x+0.9, gr.h[,x][1], col=col.br[1], border="black")
                for(i in 2:nrow(gr.h)){
                    rect(x+0.1, sum(gr.h[,x][1:(i-1)]), x+0.9, sum(gr.h[,x][1:i]), col=col.br[i], border="black"
                         )}    
                text(x+.5, 104, sum(gr.m.tax[,x]>0), cex=2.4)
            })
            par(mar = c(6, 0, 0, 0))
            plot.new()
            legend("bottomleft", legend=rev(c(names(gr.sums[gr.sums>other.thresh]), "Other")),
                   col=rev(col.a(length(gr.sums[gr.sums>other.thresh])+1)),
                   lty=1, box.lwd=3, lwd=8, bg="white", cex=3, border="white")
            
        }   
}

### Alpha diversity measurements function
### This command is quite complicated, but it basically can do all the diversity measurements together. The input is a matrix of a certain taxonmic level (specified in tax.level), and it assumes that there is an index dataframe that includes all the taxonomies (ph.fam).
### For taxonomic level (tax.level), if phylum is used, the taxonomic group has to be specified ("mtz", "non.mtz")
### Using the action function, you can specify whether to plot boxplots of the indices ("boxplot", default), make different measurements ("anova", "wilcoxon", "kruskall", "tukey") or just return the raw values from each group ("raw") based on the given diversity index.
### Diversity index in this case is either taxonomic richness ("tax.rich") or shannon diversity index ("shannon").
### If the morphological data is used, the mol.data should be switched to false. 

tax.group <- mtz
diversity.analysis <- function(tax.level="family", tax.group=mtz, phyl="Annelida", div.index="shannon", action="boxplot", p.adj="none", mol.data=TRUE, morph.matrix=pol.matrix){
    if(mol.data==FALSE){
        tax.m <- morph.matrix
    }
    else if (mol.data==TRUE){
        if(tax.level=="family"){
            tax.m <- fam.matrix[rownames(fam.matrix)%in%sort(ph.fam[ph.fam$phylum==phyl,"family"]),]
        }else if(tax.level=="genus"){
            tax.m <- gen.matrix[rownames(gen.matrix)%in%sort(ph.fam[ph.fam$phylum==phyl,"genus"]),]
        } else if(tax.level=="phylum"){
            tax.m <- ph.matrix[rownames(ph.matrix)%in%tax.group,]}
    }
    ##Choosing between diversity indices 
    if(div.index=="tax.rich"){
        ##Family richness
        tax.sum <- sapply(1:ncol(tax.m), function(x){ length(tax.m[,x][tax.m[,x]>0])})
        div.df <- as.data.frame(cbind(Fish.farm=eco.s$Fish.farm,
                                      Sampling.Site=eco.s$Sampling.Site,
                                      indx=as.numeric(tax.sum)))
        div.d <- lapply(factor.grouping, function(x){tax.sum[x]})
        names(div.d) <- colnames(fam.m.f)
    } else if(div.index=="shannon"){
        ##Family diversity
        shan.d <- sapply(1:ncol(tax.m), function(x){ diversity(tax.m[,x][tax.m[,x]>0])})
        div.d <- lapply(factor.grouping, function(x){shan.d[x]})
        div.df <- as.data.frame(cbind(Fish.farm=eco.s$Fish.farm,Sampling.Site=eco.s$Sampling.Site,indx=round(as.numeric(shan.d), 2)))
    }
    ##With raw, you get the dataframe that is used for plotting and statistical analysis
    if(action=="raw"){
        return(div.d)
    }
    ##Boxplotting
    if(action=="boxplot"){
        boxplot(div.d, col="skyblue4", lwd=2, jitter=T, las=1,
                                        #                ylab=ifelse(div.index=="shannon","Shannon diverity (H')","Taxonomic richness"),
                ylab="", cex.axis=1.8, cex.lab=1.8,
                cex.main=.85)
        ##Statistical analysis
    }else if(action=="wilcoxon"){
        groups <- rep(names(div.d), times=c(as.numeric(unlist(lapply(div.d, length)))))
        pairwise.wilcox.test(c(as.numeric(unlist(div.d))), groups, p.adjust.method =p.adj)
        ##Kruskal is for comparing multiple samples 
    }else if(action=="kruskal"){
        groups <- rep(names(div.d), times=c(as.numeric(unlist(lapply(div.d, length)))))
        kruskal.test(c(as.numeric(unlist(div.d))), groups)
    }else if(action=="anova"){
        summary(aov(indx~Fish.farm*Sampling.Site, data=div.df))            
    }else if(action=="tukey"){
        TukeyHSD(aov(indx~Fish.farm*Sampling.Site, data=div.df))
    }
}

###########################################
### Beta diversity measuerements      
### NMDS plots
### For the color and scatter type factor (in the NMDS plot)
col.factor <-  function(ds, eco.df=eco.df){
    ss <- eco.df[eco.df$Sample %in% rownames(ds),]
    sapply(1:nrow(ds), function(i){
        if(ss[i,]$Fish.farm=="Hisdalen"&ss[i,]$Sampling.Site=="Near"){
            return("tomato")}
        else if(ss[i,]$Fish.farm=="Hisdalen"&ss[i,]$Sampling.Site=="Far"){
                return("tomato4")}
        else if(ss[i,]$Fish.farm=="Onarøy"&ss[i,]$Sampling.Site=="Near"){
                    return("skyblue")}
        else if(ss[i,]$Fish.farm=="Onarøy"&ss[i,]$Sampling.Site=="Far"){
                        return("skyblue4")}
    })
}           
pch.factor <-  function(ds, eco.df=eco.df){
    ss <- eco.df[eco.df$Sample %in% rownames(ds),]
    sapply(1:nrow(ds), function(i){
        if(ss[i,]$Sampling.Site=="Near"){
            return(19)}
           else if(ss[i,]$Sampling.Site=="Far"){
                    return(18)}
       })
}

### The plotting command for NMDS, input: a matrix with rows of samples and columns of taxa 
### If pres.abs is true, it reverts all non-zero values to 1 (presence/absence values)
           
plot.nmds <- function(gr.matrix.t, pres.abs=FALSE, pres.thresh=0, phyl.name="",  pos.l.1="topleft", sqrt.tr=TRUE, pos.l.2="bottomright", st.mth="bray", iso=TRUE, main.title="", eco.df=eco.s){
    set.seed(2)
    gr.matrix.t <- gr.matrix.t[rowSums(gr.matrix.t)>0,]
    if(sqrt.tr==TRUE){
        gr.matrix.t <- sqrt(gr.matrix.t)
    }
    if (pres.abs==TRUE){
        gr.matrix.t <- ifelse((gr.matrix.t>pres.thresh)==TRUE, 1, 0)
        gr.matrix.t <- gr.matrix.t[rowSums(gr.matrix.t)>0,]
    }
    swiss.dist <- vegdist(gr.matrix.t, method=st.mth)
        swiss.mds <- isoMDS(swiss.dist)
    }
    else{swiss.mds <- metaMDS(swiss.dist, distance=st.mth)}
    par(mar=c(4,4,3,1))
    if(main.title==""){
        plot(swiss.mds$points, type="n", main=paste("NMDS,", st.mth,",",phyl.name), xlab="Coordinate 1", ylab="Coordinate 2", cex=1.5, cex.lab=1.5, cex.main=1.5)
    }
    else{
        plot(swiss.mds$points, type="n", main=main.title, xlab="Coordinate 1", ylab="Coordinate 2", cex=1.5,
             cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
    }
    stress <- paste0("stress: ", round(metaMDS(swiss.dist)$stress,2))
    lines(swiss.mds$points, type="p",  pch=pch.factor(gr.matrix.t, eco.df), col=col.factor(gr.matrix.t, eco.df), cex=3)
    text(par()$usr[1]*(ifelse(par()$usr[1]<0,0.7,1.2)), par()$usr[3]*(ifelse(par()$usr[3]<0,0.9,1.1)), stress,bg="white", font=2, cex=1.3)   
    legend(pos.l.2, legend=c("HF", "HC", "OF", "OC"), col=c("tomato","tomato4","skyblue","skyblue4"),pch=c(19,18,19,18), bg="white", cex=1.3, pt.cex=2)
}

### Get PERMANOVA values from adonis input: a matrix with rows of taxa and columns of samples 
### If pres.abs is true, it reverts all non-zero values to 1 (presence/absence values)
adonise.ph.mtrx <- function(the.matrix, sqr.tr=TRUE, pres.abs=FALSE, pres.thresh=0, eco.df=eco.s){
    ##    the.matrix <- matrix.phyl(x)
    if(pres.abs==TRUE){
        the.matrix <-ifelse((the.matrix>pres.thresh)==TRUE, 1, 0)}
    if(sqr.tr==TRUE){
        the.matrix <- sqrt(the.matrix)}
    adonis(the.matrix~Fish.farm*Sampling.Site,data=eco.df[eco.df$Sample%in%rownames(the.matrix),] )[[1]]}

########################################
##These are generally useful functions for any analyses

### For a faster reset of the display issues on ESS
set.display <- function(number){ 
    local.h <- paste0("localhost:", number,".0")
    Sys.setenv('DISPLAY'=local.h)
}
      
###Get the size of all the variables
format_bytes <- function(bytes) {
  if (bytes >= 2^30) {
    return(sprintf("%.2f GB", bytes / 2^30))
  } else if (bytes >= 2^20) {
    return(sprintf("%.2f MB", bytes / 2^20))
  } else if (bytes >= 2^10) {
    return(sprintf("%.2f KB", bytes / 2^10))
  } else {
    return(sprintf("%d bytes", bytes))
  }
}

##############Get sizes of all objects
# Get the names of all objects in the current environment
object_names <- ls()

# Loop through each object name and print its size along with the name
for (obj_name in object_names) {
  obj <- get(obj_name)  # Get the object associated with the name
  obj_size <- object.size(obj)  # Get the size of the object
  cat(obj_name, "- Size:", format_bytes(obj_size),"\n" )
}
