library(vegan)
library(scales)
library(MASS)
library(ape)
library(reshape2)
library(parallel)
library(ggplot2)
library(ggvenn)
library(RColorBrewer)

#################################################################
###1. Startup functions and variable sets
#These are the colors I most typically go for
##This is the one I use for the pies
#col.a <- colorRampPalette(c("sandybrown", "tomato","tomato4","skyblue","skyblue4", "darkslategray"))
##And that is for the other plots
col.a <- colorRampPalette(c("peachpuff", "sandybrown", "tomato","tomato3","lightblue","skyblue","skyblue4", "darkslategray"))

###List of metazoans that appear as phyla in NCBI
mtz <- sort(c("Mollusca","Echinodermata","Arthropoda", "Cnidaria", "Annelida", "Nematoda", "Nemertea", "Platyhelminthes", "Priapulida", "Brachiopoda", "Xenacoelomorpha","Hemichordata","Bryozoa", "Porifera","Sipuncula","Chaetognatha", "Tardigradia","Chordata", "Rotifera", "Gastrotricha", "Kinorhyncha", "Ctenophora", "Entoprocta"))

##Since I am working on emacs ess, this is a shortcut for when display changes. 
set.display <- function(number){ 
    local.h <- paste0("localhost:", number,".0")
    Sys.setenv('DISPLAY'=local.h)
}

####Here is a reset function for mfrow, not sure if it is faster to write it but still..
restore.mfrow <- function(x){
    par(mfrow=c(1,1))}

########################################################################
##For the sample
oversamples <- sub("([0-9]+).*", "\\1", sample.names)
#factor.names <- c("Hisdalen-farm", "Hisdalen-control", "Onarøy-farm", "Onarøy-control")
factor.names <- c("HF", "HC", "OF", "OC")
overeco <- eco[!duplicated(eco$Location),]

##########################################################################
###2. For an automatic generation of the dataframes used in the data visualization steps:
## Input file 1: BLAST tabulated file that has at least the following columns:
## query, subject, title, qlen, slen, qstart, qend, sstart, ssend, evalue, score, length, pident, family, phylum
## Taxonomic information was taken from a custom script that uses the ncbi taxonomy database information on each species.
## Input file 2: A tabulated file that contains quantitative data for each sample, this particular one was made using Salmon. Should contain the following columns:
## Name(the query), Length, EffectiveLength, TPM, NumReads

nextseq.vars <- function(pid.thresh=93, alilength=100, pid.thresh.2=95, count.mode="read.nmbr", filter=TRUE, read.nmbr.thr=1){
    ##Step 1: Use parameters to filter out low alignment quality sequences, remove queries that matched to multiple phyla, report on their proportion against the total number of matched queries.
    tax.info.hp <- tax.info.na[tax.info.na$pid>=pid.thresh & tax.info.na$length>=alilength,]
    m.ph <- sapply(unique(tax.info.hp$query), function(x){length(unique(tax.info.hp[tax.info.hp$query==x,"phylum"]))})
    print(paste0("Proportion of single phyla per contig; ", round((sum(m.ph[m.ph==1])/length(m.ph))*100, 2), "%"))
    s.ph <- m.ph[m.ph==1]
    tax.info.sph <- tax.info.hp[tax.info.hp$query %in% names(s.ph), ]
    
###Step 2: Order query matches based on score and only keep the best scoring one (this needs a bit more complexity in the future as I think it makes for a very absolute arbitration)
    tax.sc.sph <- tapply(1:nrow(tax.info.sph),tax.info.sph$query, function(i){
        i[ order(tax.info.sph[i, 'score'], decreasing=T) ]
    })
    
    tax.t.sc.sph <- sapply(tax.sc.sph, function(x){ x[1] })
tax.ts.sph <- tax.info.sph[tax.t.sc.sph,]

###Step 3: Filter our some queries based on keyword search in the title column. For example, genome and chromosome based contigs seem to create issues. Also remove query entries with both the family and phylum having "NULL" values 
    tax.ts.sph <- tax.ts.sph[grep("chromosome|genome",tax.ts.sph$title, invert=TRUE),]
    tax.ts.sph <- tax.ts.sph[tax.ts.sph$family!="NULL"|tax.ts.sph$phylum!="NULL",]

###Step 4: Merge the alignment info dataset with the quantiative one. Output will be a list of each sample. 
    l.tx.sph <- lapply(names(sal.l),function(x){
        colnames(sal.l[[x]]) <- c("query", "length", "eff.length","tpm","read.nmbr")
        merge(sal.l[[x]], tax.ts.sph[,c("query","subject", "title", "family","pident","score","phylum")],by="query")
    })

###Step 5: If required, perform a second filter step, one using alignment quality parameters and another using minimum quantity parameters. 
    list.len <- length(sample.names)
    l.tx.p <- lapply(1:list.len, function(x){l.tx.sph[[x]][l.tx.sph[[x]]$pident >= pid.thresh.2, ]})
    l.tx.p <- lapply(l.tx.p, function(x){x[x$read.nmbr>=read.nmbr.thr,]})
    
###Step 6: Remove family groups that are irrelevant to the dataset (with caution)
    
    l.tx.p <- lapply(l.tx.p, function(x){x[x$phylum!="Basidiomycota"&
                                           x$family!="Schistosomatidae"&
                                           x$family!="Rhabditidae"&
                                           x$family!="Aplysiidae"&
                                           x$family!="Malasseziaceae",]})

###Step 7: Create two variables, one that summarizes the number of alignments per phyla per sample and another for allignments per families per sample
    
    l.tax.agg.ph <- lapply(1:list.len,function(x){
        agg.ph <- aggregate(l.tx.p[[x]][,count.mode], by=list(Category=l.tx.p[[x]]$phylum), FUN=sum)
        colnames(agg.ph) <- c("phylum", "ph.sum")
        l.t <- merge(l.tx.p[[x]], agg.ph, by=c("phylum"))
        l.t[l.t$ph.sum>10,]
    })
    ##Unduplicate phylum sums
    l.tax.ph.u<- lapply(l.tax.agg.ph,function(x){x[!duplicated(x$phylum),]})
    
    l.tax.agg <- lapply(1:list.len,function(x){
        agg.f <- aggregate(l.tx.p[[x]][,count.mode], by=list(Category=l.tx.p[[x]]$family), FUN=sum)
        colnames(agg.f) <- c("family", "fam.sum")
        l.t <- merge(l.tx.p[[x]], agg.f, by=c("family"))
        l.t
        })
    
    ##Only keep one family per entry (since they all should have the sams read).
    l.tax.u<- lapply(l.tax.agg,function(x){x[!duplicated(x$family),]})
    
    
    ##Make an array with all the phyla 
    all.ph <- sort(unique(unlist(lapply(l.tax.ph.u, function(x){x$phylum}))))
   
    ##And an array with all the families
    all.fams <- sort(unique(unlist(lapply(l.tax.u, function(x){unique(x$family)}))))
    ##Step 8: Export the data by creating a list of all the lists and variables
    output.l <- list(l.tx.p, l.tax.agg.ph, l.tax.agg, l.tax.u, l.tax.ph.u, all.ph, all.fams)
    names(output.l) <- c("l.tx.p", "l.tax.agg.ph", "l.tax.agg", "l.tax.u", "l.tax.ph.u", "all.ph", "all.fams")
    output.l
    ##To execute this function, you should write it like this :
    ##outputs <- nextseq.vars(95, 100, 95)
    ##And then for each output you output the data like this:
                                        #l.tx.p <- outputs[[1]]
                                        #l.tax.agg.ph<- outputs[[2]]
                                        #l.tax.agg <- outputs[[3]]
                                        #l.tax.u <- outputs[[4]]
                                        #l.tax.ph.u <- outputs[[5]]
                                        #all.ph <- outputs[[6]]
                                        #all.fams <- outputs[[7]]
}

###Function below is used to sort each sample into the corect sampling location using a reference dataset with location/environmental info

fact.indx <- function(x, data.eco=eco){
    list(Hisd_farm=x[which(data.eco$Fish.farm=="Hisdalen" & data.eco$Sampling.Site=="Near")],
         Hisd_control=x[which(data.eco$Fish.farm=="Hisdalen" & data.eco$Sampling.Site!="Near")],
         Onar_farm=x[which(data.eco$Fish.farm != "Hisdalen" & data.eco$Sampling.Site=="Near")],
         Onar_control=x[which(data.eco$Fish.farm!="Hisdalen" & data.eco$Sampling.Site!="Near")])
}

## 3. Function for creating pie charts
#Input: A dataset (l.tax.u) that contains the sum of all alignments per family per sample (see nextseq.vars function)
#"famsum" variable sets whether we want to count the number of families detected per phylum (fam.sum= FALSE) or count the number of alignments per phyla
#"famcount" vaiable indicates a minimum number of counts per family alignment

molto.pies <- function(tax.group=mtz, fam.count=1, other.thresh=3, famsum=FALSE){
    if(famsum==FALSE){
        ##Get all the unique families 
        all.groups <- cbind(unlist(lapply(l.tax.u, function(x){x[x$fam.sum>fam.count,c("phylum")]})),unlist(lapply(l.tax.u, function(x){x[x$fam.sum>fam.count,c("family")]})))
        all.groups <- all.groups[!duplicated(all.groups[,2]),]
        ##Enumerate the families and change them into proportions
        ph.count.abs <- rev(sort(table(all.groups[all.groups[,1]%in%tax.group,1])))
        ph.count <- ph.count.abs/sum(ph.count.abs)*100
        ##Create an "other" group for all the groups whose proportion is smaller than "other.thresh", to avoid name overal on the pie
        other.name <- paste0("Other (n=", length(ph.count[ph.count<=other.thresh]),")")
        ph.count.h <- c(ph.count[ph.count>other.thresh], Other=100-sum(ph.count[ph.count>other.thresh]))
        names(ph.count.h)[length(ph.count.h)] <- other.name
        ##Pie ploting
        par(mar=c(5,5,5,5))
        pie(ph.count.h, labels=paste0(names(ph.count.h), "\n", round(ph.count.h,1),"%"),
            col=col.a(length(ph.count.h)), las=2, cex=1.5)
        text(-0.7,-1,paste("Total families:", sum(ph.count.abs)), font=2, cex=2)
    }
    else if(famsum==TRUE){
        ##Get the sum of alignments per phylum and change them into proportions

        all.ph.sum <- unlist(lapply(l.tax.u, function(x){
            sapply(tax.group, function(y){
                sum(x[x$phylum==y&x$fam.sum>fam.count,"fam.sum"])})
        }))
        al.ph.count <- sapply(tax.group, function(x){
            sum(all.ph.sum[names(all.ph.sum)==x])
        })
        al.ph.count <- sort(al.ph.count[al.ph.count>0])/sum(al.ph.count)*100
        ##Create an "other" group for all the groups whose proportion is smaller than "other.thresh", to avoid name overal on the pie
        other.name <- paste0("Other (n=", length(al.ph.count[al.ph.count<=other.thresh]),")")
        al.ph.count.h <- c(al.ph.count[al.ph.count>other.thresh], Other=sum(al.ph.count[al.ph.count<=other.thresh]))
        names(al.ph.count.h)[length(al.ph.count.h)] <- other.name
        ##Pie ploting
        par(mar=c(5,5,5,5))
        pie(al.ph.count.h, labels=paste0(names(al.ph.count.h), "\n", round(al.ph.count.h,1),"%"),col=rev(col.a(length(al.ph.count.h))),cex=1.5)
        text(-0.7,-1,paste("Total reads:", round(sum(all.ph.sum),2)), font=2, cex=2)
    }
}

## 4.Creating barplots of taxonomic distributions per sampling location

###4.1 : A function that makes family distribution barplots per sampling location for the molecular datasets. Variables, "selected.ph" phylum of interest, e.g.: "Annelida", "other.thresh" the proportion for groups to be inserted in the other "other" category, "eco.df", the dataset that is used as a reference to index the samples to the correct sampling location
bar.ph <- function(selected.ph, other.thresh=1, eco.df=eco){
    factor.names <- c("HF", "HC", "OF", "OC")
    factor.list <- fact.indx(1:nrow(eco.df), data.eco=eco.df)
###Keep samples that had entries for the selected phylum
    m.factor.list <- lapply(1:length(factor.list), function(y){
        factor.list[[y]][which(unlist(lapply(factor.list[[y]], function(x){
            nrow(l.tx.p[[x]][l.tx.p[[x]]$read.nmbr> 1 & l.tx.p[[x]]$phylum%in%selected.ph,])
        }))!=0)]
    })
###Create a list of family abundances per sampling location
    factor.distr <- lapply(m.factor.list, function(f.l){
        fact.l <- lapply(f.l, function(x){
            a <- l.tx.p[[x]][l.tx.p[[x]]$read.nmbr>=1 & l.tx.p[[x]]$phylum%in%selected.ph,]
            cc <- aggregate(a[,"read.nmbr"], by=list(a[,"family"]),sum)
            colnames(cc) <- c("family","count")
            cc$perc <- (cc$count/sum(cc[,"count"]))*100
            cc})
        fact.names <- sort(unique(unlist(lapply(fact.l, function(x){x$family}))))
        distr <- sapply(fact.names, function(y){
            sum(unlist(lapply(1:length(fact.l), function(x){
                fact.l[[x]][fact.l[[x]]$family==y,"count"]}
                )))
        })
        distr
    })
    ##Proportionalize the distribution as well as created an "Other" category for alignments that fall below the theshold number    
    factor.sum <- lapply(factor.distr, function(distr){
        distr.h <- distr[distr/sum(distr)*100>=other.thresh]
        distr.h <- c(distr.h[order(names(distr.h))], Other=sum(distr[distr/sum(distr)*100<other.thresh]))
        distr.h/sum(distr.h)*100
    })
    ##Plot the distributions into a barplot
    par(mar=c(10,5,5,1))
    plot(1:10, xlim=c(1,7), ylim=c(0,107), type="n", ylab="Relative abundance (%)", axes=F, xlab="", xaxs="i",yaxs="i", cex.main=1.5, cex.axis=1.4, cex.lab=1.4)
    axis(1, seq(1.5, length(factor.names)+0.5,1), labels=factor.names, cex.axis=1.4, tick=F, font.axis=1.8, las=1)
    axis(2, seq(0, 100, 20), labels=(seq(0,100,20)), las=2, cex.axis=1.5)
    abline(h=seq(20, 100, 20), col="grey", lwd=3)
    group.names<-sort(unique(names(unlist(lapply(factor.sum,function(x){x})))))
    ##This is a bit convoluted but it is here for when there are no other variables created in the family distribution
    others=length(unique(unlist(lapply(l.tx.p,function(x){x[x$phylum==selected.ph&x$read.nmbr>1, "family"]}))))-length(group.names)-1
    if(length(others)>0){
        group.names <- c(group.names[group.names!="Other"],"Other")}
    legend("bottomright", legend=rev(group.names), col=rev(col.a(length(group.names))), lty=1, lwd=8, bg="white",cex=1.1, border="white")
    sapply(1:length(factor.sum), function(x){
        ##This color variable makes it so that the coloring of each family is consistent among the sample locations
        col.ph=col.a(length(group.names))[group.names %in% names(factor.sum[[x]])]
        rect(x+0.1, 0, x+0.9, factor.sum[[x]][1], col=col.ph[1], border="black")
        for(i in 2:length(factor.sum[[x]])){
            rect(x+0.1, sum(factor.sum[[x]][1:(i-1)]), x+0.9, sum(factor.sum[[x]][1:i]), col=col.ph[i], border="black"
                 )}    
        text(x+.5, 104, length(factor.distr[[x]]))
    })    
    box(lty=1, lwd=3)
}

### 4.2: A function that makes family distribution barplots per sampling location for the morphological datasets. Variables, "selected.ph" phylum of interest, e.g.: "Annelida", "aggtype": aggregation factor, usually count or biomass, "other.thresh" the proportion for groups to be inserted in the other "other" category.

bar.morpho <- function(selected.ph, aggtype="count", other.thresh=1, eco.df=eco.s){
    factor.list <- fact.indx(1:nrow(eco.df), data.eco=eco.df)
    m.factor.list <- lapply(1:length(factor.list), function(y){
        factor.list[[y]][which(unlist(lapply(factor.list[[y]], function(x){
            nrow(morpho.l[[x]][morpho.l[[x]][,aggtype]> 0 & morpho.l[[x]]$phylum==selected.ph,])
        }))!=0)]
    })
    ###Create a list of family abundances per sampling location
    factor.distr <- lapply(m.factor.list, function(f.l){
        fact.l <- lapply(f.l, function(x){
            a <- morpho.l[[x]][morpho.l[[x]][,aggtype]>0 & morpho.l[[x]]$phylum==selected.ph,]
            cc <- aggregate(a[,aggtype], by=list(a[,"family"]),sum)
            colnames(cc) <- c("family","count")
            cc$perc <- (cc$count/sum(cc[,"count"]))*100
            cc})
        fact.names <- sort(unique(unlist(lapply(fact.l, function(x){x$family}))))
        distr <- sapply(fact.names, function(y){
            sum(unlist(lapply(1:length(fact.l), function(x){
                fact.l[[x]][fact.l[[x]]$family==y,"count"]}
                )))
        })
        distr
    })
    ##Proportionalize the distribution as well as created an "Other" category for alignments that fall below the theshold number    
    factor.sum <- lapply(factor.distr, function(distr){
        distr.h <- distr[distr/sum(distr)*100>=other.thresh]
        distr.h <- c(distr.h[order(names(distr.h))], Other=sum(distr[distr/sum(distr)*100<other.thresh]))
        distr.h/sum(distr.h)*100
    })
    ##Plot the distributions into a barplot
    par(mar=c(10,5,5,1))
    plot(1:10, xlim=c(1,7), ylim=c(0,107), type="n", ylab="Relative abundance (%)", axes=F, xlab="", xaxs="i",yaxs="i", cex.main=1.5, cex.axis=1.4, cex.lab=1.4, main=paste(selected.ph, "Morpho"))
    axis(1, seq(1.5, length(factor.names)+0.5,1), labels=factor.names, cex.axis=1.4, tick=F, font.axis=1.8, las=1)
    axis(2, seq(0, 100, 20), labels=(seq(0,100,20)), las=2, cex.axis=1.5)
    abline(h=seq(20, 100, 20), col="grey", lwd=3)
    group.names<-sort(unique(names(unlist(lapply(factor.sum,function(x){x})))))
    others=length(unique(unlist(lapply(morpho.l,function(x){x[x$phylum==selected.ph&x[,aggtype]>1, "family"]}))))-length(group.names)-1
    if(length(others)>0){
        group.names <- c(group.names[group.names!="Other"],"Other")}
    legend("bottomright", legend=rev(group.names), col=rev(col.a(length(group.names))), lty=1, lwd=8, bg="white",cex=1.1, border="white")
    sapply(1:length(factor.sum), function(x){
        col.ph=col.a(length(group.names))[group.names %in% names(factor.sum[[x]])]
        rect(x+0.1, 0, x+0.9, factor.sum[[x]][1], col=col.ph[1], border="black")
        for(i in 2:length(factor.sum[[x]])){
            rect(x+0.1, sum(factor.sum[[x]][1:(i-1)]), x+0.9, sum(factor.sum[[x]][1:i]), col=col.ph[i], border="black"
                 )}    
        text(x+.5, 104, length(factor.distr[[x]]))
    })    
    box(lty=1, lwd=3)
}

##4.3 Very similar barplot distribution to bar.ph() but this one is for metazoan and non-metazoans
bar.gr <- function(group, eco.df=eco, o.t=2){
    not.empty <- unlist(lapply(1:length(l.tax.u), function(x){nrow(l.tax.u[[x]][l.tax.u[[x]]$phylum%in%group,])>0}))
    factor.list <- fact.indx(1:nrow(eco.df), data.eco=eco.df)
    m.factor.list <- lapply(factor.list, function(x){x[x%in%which(not.empty==T)]})
    factor.distr <- lapply(m.factor.list, function(f.l){
        fact.l <- lapply(f.l, function(x){
            a <- l.tax.u[[x]][l.tax.u[[x]]$fam.sum>0,]
            group.phyl <- a[a$phylum %in% group,]
            cc <- aggregate(group.phyl[,"fam.sum"], by=list(group.phyl[,"phylum"]),sum)
            colnames(cc) <- c("phylum","ph_sum")
            cc$perc <- (cc$ph_sum/sum(cc[,"ph_sum"]))*100
            cc
        })
        fact.names <- sort(unique(unlist(lapply(fact.l, function(x){x$phylum}))))
        distr <- sapply(fact.names, function(y){
            sum(unlist(lapply(1:length(fact.l), function(x){
            fact.l[[x]][fact.l[[x]]$phylum==y,"ph_sum"]}
            )))
        })
        distr
    })
    factor.sum <- lapply(factor.distr,function(distr, other.thresh=o.t){
        distr.h <- distr[distr/sum(distr)*100>=other.thresh]
        distr.h <- c(distr.h[order(names(distr.h))], Other=sum(distr[distr/sum(distr)*100<other.thresh]))
        distr.h/sum(distr.h)*100
    })
    par(mar=c(10,5,5,1))
    plot(1:10, xlim=c(1,7), ylim=c(0,107), type="n", ylab="Relative abundance (%)", axes=F, xlab="", xaxs="i",yaxs="i", cex.main=1.5, cex.axis=1.4, cex.lab=1.4)
    axis(1, seq(1.5, length(factor.names)+0.5,1), labels=c("HF", "HC", "OF", "OC"), cex.axis=1.4, tick=F, font.axis=1.8, las=1)
    axis(2, seq(0, 100, 20), labels=(seq(0,100,20)), las=2, cex.axis=1.5)
    abline(h=seq(20, 100, 20), col="grey", lwd=3)
    group.names<-sort(unique(unlist(lapply(factor.sum,function(x){names(x)}))))
    group.names <- c(group.names[group.names!="Other"],"Other")
    legend("bottomright", legend=rev(group.names), col=rev(col.a(length(group.names))), lty=1, lwd=8, bg="white", cex=1.2, border="white")
    sapply(1:length(factor.sum), function(x){
        col.ph=col.a(length(group.names))[group.names %in% names(factor.sum[[x]])]
        rect(x+0.1, 0, x+0.9, factor.sum[[x]][1], col=col.ph[1], border="black")
        for(i in 2:length(factor.sum[[x]])){
            rect(x+0.1, sum(factor.sum[[x]][1:(i-1)]), x+0.9, sum(factor.sum[[x]][1:i]), col=col.ph[i], border="black"
                 )}    
        text(x+.5, 104, length(factor.distr[[x]]), cex=1.3)
    })
    box(lty=1, lwd=3)
}

###############################################################
### 5. Alpha diversity metrics
##5.1 Getting shapiro restults based on abundance of a specific taxonomic phylum (probably wrong way of doing it tough..)

shapiro.abundance <- function(tx.gr){
    read.sum.gr <- lapply(l.tx.p, function(x){
        a <- x[x$read.nmbr>0 & x$phylum%in%tx.gr,"read.nmbr"]
        a})
    shapiro.test(unlist(lapply(read.sum.gr, sum)))
}

##5.2 Getting diversity measurements and plots, This is a complex function that can provide a few things: i)shannon H diversity value, simpson diversity values, family richness values, and boxplots for each such function. Uses l.tax.u famsum as input

##Options:
##tax.group: phylogenetic group,
##div.index: diversity index used, options are famsum (sum of total families per location) readsum(mean/sum of reads in location, Shannon (for shannon H' diversity), set to shannon,
##tax.group.name: group specified in the plot title/axis,
##select.fun: function selection for the read sum parameter (mean, sum ect)
##action: a few options here again, i)raw: provides the direct results from the metrics by the diversity index, ii) boxplot: returns a boxplot of the four locations, iii) wilcoxon/kruskall for non-parametric results, anova/tukey for parametric statistical results
## main.in: include a main title, set to true,
##p.adj: this is for adjusting the p values for the wilcoxon test (eg if bonferroni should be used), set to none

##Some examples of syntax
##For getting tukey of family sums
##    diversity.core.l(mtz, div.index="famsum", action="tukey")$`Fish.farm:Sampling.Site`,
##For getting anova stats for shannon
##   as.data.frame(diversity.core.l(not.mtz,  action="anova")[[1]])
##For getting boxplots of family richness of Annelids
##   diversity.core.l("Annelida", div.index="famsum", tax.group.name="Annelida", action="boxplot")

diversity.core.l <- function(tax.group, div.index="shannon", tax.group.name="", select.fun=sum, action="boxplot", main.in=T, p.adj="none"){
    if(div.index=="famsum"){
        sh.indx.gr <- lapply(l.tax.u, function(x){
            a <- x[x$fam.sum>0 & x$phylum%in%tax.group,"family"]
            a
        })
        sh.indx.l <- sapply(1:nrow(eco.s),function(y){length(unique(unlist(lapply(y, function(x){sh.indx.gr[[x]]}))))})
        sh.mol.l <- fact.indx(sh.indx.l,eco.s)
        aov.df <- as.data.frame(cbind(Fish.farm=eco.s$Fish.farm,Sampling.Site=eco.s$Sampling.Site, indx=as.numeric(sh.indx.l)))
    }
    else if(div.index=="readsum"){
        sh.indx.gr <- lapply(l.tax.u, function(x){
            a <- x[x$read.nmbr>0 & x$phylum%in%tax.group,"fam.sum"]
            a
        })
        sh.indx.l <- sapply(sh.indx.gr, select.fun)
        sh.mol.l <- fact.indx(sh.indx.l,eco.s)
        aov.df <- as.data.frame(cbind(Fish.farm=eco.s$Fish.farm,Sampling.Site=eco.s$Sampling.Site, indx=as.numeric(sh.indx.l)))
    }
    else if(div.index=="shannon"){
        sh.indx.gr <- unlist(lapply(l.tax.u, function(x){round(diversity(x[x$phylum%in%tax.group,"fam.sum"], index=div.index),2)}))
        sh.indx.loc <- data.frame(sh.indx=sh.indx.gr, loc=eco.s$Sample)
        sh.indx.loc <- sh.indx.loc[sh.indx.loc$sh.indx>-1,]
        sh.mol.l <- fact.indx(sh.indx.loc$sh.indx,eco.s)
        aov.df <- as.data.frame(cbind(Fish.farm=eco.s$Fish.farm,Sampling.Site=eco.s$Sampling.Site, indx=sh.indx.gr))
    }
    if(action=="raw"){
        return(sh.mol.l)
    }
    else if(action=="boxplot"){
        if(main.in==T){
            boxplot(sh.mol.l, names=paste(names(sh.morph.pol.l)), col="skyblue4", lwd=2,
                    main=ifelse(div.index=="shannon",paste("Shannon index", tax.group.name, pid.thresh, ",", alilength),paste("Family richness", tax.group.name, pid.thresh, ",", alilength)), jitter=T, las=2, cex.axis=2, cex.main=.85)
        }else{
            boxplot(sh.mol.l, names=paste(names(sh.morph.pol.l)), col="skyblue4", lwd=2, jitter=T, las=1, ylab="Shannon index", cex.main=.85)}
    }else if(action=="wilcoxon"){
        groups <- rep(names(sh.mol.l), times=c(as.numeric(unlist(lapply(sh.mol.l, length)))))
        pairwise.wilcox.test(c(as.numeric(unlist(sh.mol.l))), groups, p.adjust.method =p.adj)
###Kruskal is for comparing multiple samples 
    }else if(action=="kruskal"){
        groups <- rep(names(sh.mol.l), times=c(as.numeric(unlist(lapply(sh.mol.l, length)))))
        kruskal.test(c(as.numeric(unlist(sh.mol.l))), groups)
    }else if(action=="anova"){
        summary(aov(indx~Fish.farm*Sampling.Site, data=aov.df))            
    }else if(action=="tukey"){
        TukeyHSD(aov(indx~Fish.farm*Sampling.Site, data=aov.df))
    }
}

##################################################
### 6. NMDS plots:

## 6.1 Function that creates an abundance matrix based on family read sums of a particular phylum for each sample
###Input dataframe is l.tax.u because it has the fam sum.column.
matrix.phyl <- function(phyl, sample.list=l.tax.u, lt.len=length(l.tax.u), sample.n=sample.names){
    phyl.u <- lapply(1:lt.len, function(x){
        sample.list[[x]][sample.list[[x]]$phylum%in%phyl,]})
    all.fams <- sort(unique(unlist(lapply(1:length(phyl.u), function(x){
        unique(phyl.u[[x]]$family)}
        ))))
    gr.matrix <- matrix(0, nrow=length(all.fams), ncol=lt.len)
    rownames(gr.matrix) <- all.fams
    colnames(gr.matrix) <- sample.n
    ##Start with the phylum sums..
    for(i in 1:lt.len){
        for(j in phyl.u[[i]]$family){
            if(which(rownames(gr.matrix)==j)>0){
                gr.matrix[which(rownames(gr.matrix)==j),i] <- phyl.u[[i]][phyl.u[[i]]$family==j, "fam.sum"]
            }}
    }
    gr.matrix.t <- t(gr.matrix)
    gr.matrix.t <- gr.matrix.t[rowSums(gr.matrix.t)>0,]
    gr.matrix.t
}

### 6.2
###Plot nmds, it provides processing of data and plot options for nmds analysis. Input, gr.matrix.t (tranposed), i.e.: matrix containing abundance  information per location per sampling site.
###options:
## pres.abs: use only perence absence information
## pres.thresh: minimum abundance for a species to be counted as present
## phyl.name: phylum in question for the plot title/axis
## pos.l.1/pos.l.2: where to position the legends on the plot so as to avoid colliding with the points
## sqrt.tr: square root transform data
## st.mth: nmds statistical method used, set to bray, can be jaccard
## iso: whether to use isoMDS or metaMDS (very slight differences)
## main.title: title of the plto
## eco.df: ecology reference tab option

plot.nmds.phyl.mf <- function(gr.matrix.t, pres.abs=FALSE, pres.thresh=0, phyl.name="", pos.l.1="topleft", sqrt.tr=TRUE, pos.l.2="bottomright", st.mth="bray", iso=TRUE, main.title="", eco.df=eco){
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
    if(iso==TRUE){
        swiss.mds <- isoMDS(swiss.dist)
    }
    else{swiss.mds <- metaMDS(swiss.dist, distance=st.mth)}
    par(mar=c(4,4,3,1))
    if(main.title==""){
        plot(swiss.mds$points, type="n", main=paste("NMDS,", st.mth,",",phyl.name), xlab="Coordinate 1", ylab="Coordinate 2", cex=1.5, cex.main=1.5)
    }
    else{
        plot(swiss.mds$points, type="n", main=main.title, xlab="Coordinate 1", ylab="Coordinate 2", cex=1.5, cex.main=1.5, cex.axis=1.5)
    }
    ##Add adonis pvalues
    a <- paste0("Fish farm: ", adonis(gr.matrix.t~Fish.farm*Sampling.Site,data=eco.df[eco.df$Sample%in%rownames(gr.matrix.t),] )[[1]][1,][6])
    b <- paste0("Distance: ", adonis(gr.matrix.t~Fish.farm*Sampling.Site,data=eco.df[eco.df$Sample%in%rownames(gr.matrix.t),] )[[1]][2,][6])
    cc <- paste0("F*D: ", adonis(gr.matrix.t~Fish.farm*Sampling.Site,data=eco.df[eco.df$Sample%in%rownames(gr.matrix.t),] )[[1]][3,][6])
                                        #    stress <- paste0("stress: ", round(swiss.mds$stress,2))
    stress <- paste0("stress: ", round(metaMDS(swiss.dist)$stress,2))
    lines(swiss.mds$points, type="p",  pch=col.factor.pch(gr.matrix.t, eco.df), col=col.factor.mf(gr.matrix.t, eco.df), cex=3)
    text(par()$usr[1]*(ifelse(par()$usr[1]<0,0.7,1.2)), par()$usr[3]*(ifelse(par()$usr[3]<0,0.9,1.1)), stress,bg="white", font=2, cex=1.3)   
    legend(pos.l.2, legend=c("HF", "HC", "OF", "OC"), col=c("tomato","tomato4","skyblue","skyblue4"),pch=c(19,18,19,18), bg="white", cex=1.4)
}

###Associated color swatches function:
col.factor.mf <-  function(ds, eco.df=eco.df){
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

###Associated point type selection function:
col.factor.pch <-  function(ds, eco.df=eco.df){
    ss <- eco.df[eco.df$Sample %in% rownames(ds),]
    sapply(1:nrow(ds), function(i){
        if(ss[i,]$Sampling.Site=="Near"){
            return(19)}
           else if(ss[i,]$Sampling.Site=="Far"){
                    return(18)}
       })
}

####5.3 Adonis analysis
###Use adonis statistic, options:
##sqrt.tr: square root transform data
##pres.abs: use presence absence data with pre.thresh determining the minimum abundance
##eco: provide a ecology reference table

adonise.ph.mtrx <- function(the.matrix, sqr.tr=TRUE, pres.abs=FALSE, pres.thresh=0, eco.df=eco){
    if(pres.abs==TRUE){
        the.matrix <-ifelse((the.matrix>pres.thresh)==TRUE, 1, 0)}
    if(sqr.tr==TRUE){
        the.matrix <- sqrt(the.matrix)}
    adonis(the.matrix~Fish.farm*Sampling.Site,data=eco.df[eco.df$Sample%in%rownames(the.matrix),] )[[1]]}


#####################################################################################
########################
### 6: Misc
### 6.1: Rarefaction curve analysis found on the internet
quickRareCurve <- function (x, step = 1, sample, xlab = "Sample Size",
                            ylab = "Species", label = TRUE, col, lty, max.cores = T, nCores = 150, ...){
    require(parallel)
    x <- as.matrix(x)
    if (!identical(all.equal(x, round(x)), TRUE))
        stop("function accepts only integers (counts)")
    if (missing(col))
        col <- par("col")
    if (missing(lty))
        lty <- par("lty")
    tot <- rowSums(x) # calculates library sizes
    S <- specnumber(x) # calculates n species for each sample
    if (any(S <= 0)) {
        message("empty rows removed")
        x <- x[S > 0, , drop = FALSE]
        tot <- tot[S > 0]
        S <- S[S > 0]
    } # removes any empty rows
    nr <- nrow(x) # number of samples
    col <- rep(col, length.out = nr)
    lty <- rep(lty, length.out = nr)
    # parallel mclapply
    # set number of cores
    mc <- getOption("mc.cores", ifelse(max.cores, detectCores(), nCores))
    message(paste("Using ", mc, " cores"))
    out <- mclapply(seq_len(nr), mc.cores = mc, function(i) {
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i])
            n <- c(n, tot[i])
        drop(rarefy(x[i, ], n))
    })
    Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
    Smax <- sapply(out, max)
     plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
       type = "n", ...)
    if (!missing(sample)) {
      abline(v = sample)
      rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),
         y = z, xout = sample, rule = 1)$y)
      abline(h = rare, lwd = 0.5)
      }
    for (ln in seq_along(out)) {
      N <- attr(out[[ln]], "Subsample")
      lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
    }
    if (label) {
        ordilabel(cbind(tot, S), labels = rownames(x), ...)
    }
    invisible(out)
}

## 6.2 Getting sized for objects in the R environment: useful when working in an environment for too long and, since there is no R studio fancy interface, I do not know which object takes up the most space
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
