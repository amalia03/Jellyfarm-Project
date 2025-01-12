#source("tax_viz_functions.R")

###Set a list of metazoans
mtz <- sort(c("Mollusca","Echinodermata","Arthropoda", "Cnidaria", "Annelida", "Nematoda", "Nemertea", "Platyhelminthes", "Priapulida", "Brachiopoda", "Xenacoelomorpha","Hemichordata","Bryozoa", "Porifera","Sipuncula","Chaetognatha", "Tardigradia","Chordata", "Rotifera", "Gastrotricha", "Kinorhyncha", "Ctenophora", "Entoprocta"))
factored <- fact.indx(1:14, eco.s)

######
###First, set up the morphological data
morpho <- read.delim("jelly_morpho.tsv", header=T,sep="\t", stringsAsFactors=FALSE, quote= "" )
colnames(morpho) <- c("station", "sample.loc", "phylum","species", "family", "count", "biomass")
morpho.a <- morpho[!is.na(morpho$count),]

###Fix some issies with the data..
morpho.fams <- morpho.a[, "family"]
morpho.fams <- gsub("[^A-Za-z]+", "", morpho.fams)
morpho.a$family <- morpho.fams

morpho.a <- cbind(Genus=unlist(lapply(strsplit(morpho.a$species, " "), function(x){paste(x[1])})), morpho.a)

##Remove the sample that has no mol caounterpart (180)
morpho.a <- morpho.a[morpho.a$sample.loc!="180",]
morpho.a <- morpho.a[morpho.a$sample.loc!="195",]
morpho.a <- morpho.a[morpho.a$sample.loc!="190",]

##Fix some mislabelled groups
morpho.a[morpho.a$Genus=="Hetromastus","Genus"] <- "Heteromastus"
morpho.a[morpho.a$Genus=="Aphaelochaeta","Genus"] <- "Aphelochaeta"
morpho.a[morpho.a$Genus=="Terebllides","Genus"] <- "Terebellides"
morpho.a[morpho.a$Genus=="Terrebelides","Genus"] <- "Terebellides"

length(unique(morpho.a$sample.loc))

sh.agg.f <-function(dat, phylum, indx="shannon")
    lapply(unique(dat$sample.loc), function(x){
        dat.sample <- dat[dat$sample.loc==x,]
        dat.group <- dat.sample[dat.sample$phylum==phylum,]
###        dat.agg.fam <- aggregate(dat.group$count, by=list(Category=dat.group$family), FUN=sum)
###        colnames(dat.agg.fam) <- c("family", "fam.sum")
        dat.agg.fam <- aggregate(dat.group$count, by=list(Category=dat.group$Genus), FUN=sum)
        dat.agg.fam <- aggregate(dat.group$biomass, by=list(Category=dat.group$Genus), FUN=sum)
        colnames(dat.agg.fam) <- c("Genus", "genus.sum")     
        dat.agg.fam
    })

agg.pol <- sh.agg.f(morpho.a, "Polychaeta", indx="shannon")
names(agg.pol) <- unique(morpho.a$sample.loc)


##Make a family abundance matrix
###all.pol <- unique(sort(unlist(lapply(agg.pol, function(x){x$family}))))
all.pol <- unique(sort(unlist(lapply(agg.pol, function(x){x$Genus}))))
pol.matrix <- matrix(0, nrow=length(all.pol), ncol=nrow(eco.s))
rownames(pol.matrix) <- all.pol
colnames(pol.matrix) <- eco.s$Sample

for(i in 1:length( agg.pol)){
    for(j in all.pol){
        if(which(rownames(pol.matrix)==j)>0){
###            a <- as.numeric( agg.pol[[i]][agg.pol[[i]]$family==j,"fam.sum"])
            a <- as.numeric( agg.pol[[i]][agg.pol[[i]]$Genus==j,"genus.sum"])
            pol.matrix[which(rownames(pol.matrix)==j),i]  <- ifelse(length(a)>0,a,0)
        }
    }}

pol.m.f <- do.call(cbind,lapply(factored, function(x){rowSums(pol.matrix[,x])}))

##FORAMS
####Add the morphological Foraminifera data
forams <- read.delim("foram_morpho_updated.csv", sep=",")
forams <- forams[1:147, 1:18]
forams[is.na(forams)] <- 0
forams <- cbind(Genus=unlist(lapply(strsplit(forams$Species, " "), function(x){paste(x[1])})), forams)

forams.agg <- sapply(unique(forams$Family), function(x){
    colSums(forams[forams$Family==x,4:length(forams)])})

##Since there are some entries that are not exactly genera, since they are all lowercase, I remove them)
forams <- forams[grepl("^[A-Z]", forams$Genus),]
forams.agg <- sapply(unique(forams$Genus), function(x){
   colSums(forams[forams$Genus==x,4:length(forams)])})



for.l <- strsplit(rownames(forams.agg), "X")
rownames(forams.agg) <- unlist(lapply(for.l, function(x){x[[2]]}))
##Remove the sampling location  that we are not using
forams.agg <- forams.agg[rownames(forams.agg)!="180",]
forams.agg <- forams.agg[rownames(forams.agg)!="190",]
###Remove the null values
forams.agg <- forams.agg[,colnames(forams.agg)!="NULL"]
##Put the samples in a numeric order
forams.agg <- forams.agg[order(rownames(forams.agg)),]
###Remove columns with 0 values
foram.matrix <- t(forams.agg[, colSums(forams.agg)!=0])
foram.m.f <- do.call(cbind, lapply(factored, function(x){rowSums(foram.matrix[,x])}))

##################################################################################################################
####Molecular data
#####################################################
##Import the salmon quantifiedvalues
sal.files <- list.files("/home/genomics/ama/data/nextseq/trinity_sans_humsalm/salmon_out","quant.sf", recursive=T)
sal.l <- lapply(1:length(sal.files), function(x){
    read.delim(paste0('/home/genomics/ama/data/nextseq/trinity_sans_humsalm/salmon_out/', sal.files[x]), header=T,sep="\t", stringsAsFactors=FALSE, quote= "" )
})

sample.names <- sub("([0-9A-Z]+).*", "\\1",sal.files)
names(sal.l) <- sample.names

###Import the ecologi cal data sheet
eco<- read.delim("../eco_info.tsv", sep="\t", header=T, stringsAsFactors=FALSE, quote="")
#Make a column that contains only the replicate information
eco$Location <- sub("([0-9]+).*", "\\1", eco$Sample)
oversamples <- sub("([0-9]+).*", "\\1", sample.names)
core.l <- c(2 ,5, 8, 11, 12, 14, 16, 20, 22, 26, 27, 29, 33, 36)
eco.s <- eco[core.l, ]
colnames(eco.s) <- c("Replicate", "Sampling.Site", "Depth", "Fish.farm", "Sample")

overeco <- eco[!duplicated(eco$Location),]

##And make a variable with the factor names
factor.names <- c("HF", "HC", "OF", "OC")

##Now import the taxonomic information and pick the contig match with the highest score
system("ls /home/genomics/ama/data/nextseq/trinity_sans_humsalm/blast_sanstrin/")
tax.info <- read.delim("/home/genomics/ama/data/nextseq/trinity_sans_humsalm/blast_sanstrin/trinity_vs_nt_euk_tax.tsv", sep="\t", header=F, stringsAsFactors=FALSE, quote="")

##Name the columns
colnames(tax.info) <- c('query', 'subject', 'title', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'ssend','evalue', 'score', 'length', 'pident', 'nident', 'gapopen', 'gaps', 'qcovs', 'junk','family', 'phylum')

###Remove Family entries without taxonomic information (but keep the phyla)
tax.info.na <- tax.info[tax.info$family!="NULL",]
tax.info.na <- tax.info.na[tax.info.na$family != "", ]
1-nrow(tax.info.na)/nrow(tax.info)

null.ann <- read.csv("Null_annotation.csv", header=F)
null.ann <- unique(null.ann)
colnames(null.ann) <- c("family", "phylum")

for (i in null.ann$family){
    tax.info.na$phylum[tax.info.na$family==i] <- null.ann[null.ann$family==i,"phylum"]
}

#Now set a list of non-metazoans
not.mtz <- unique(tax.info.na$phylum[tax.info.na$phylum%in%mtz==FALSE])

pid.thresh <- 95
alilength <- 100
pid.thresh.2 <- 95
count.mode="read.nmbr"
filter=TRUE
read.nmbr.thr <- 1

####Import taxonomic annotations. 
tax.info.hp <- tax.info.na[tax.info.na$pid>=pid.thresh & tax.info.na$length>=alilength,]
m.ph <- sapply(unique(tax.info.hp$query), function(x){length(unique(tax.info.hp[tax.info.hp$query==x,"phylum"]))})
print(paste0("Proportion of single phyla per contig; ", round((sum(m.ph[m.ph==1])/length(m.ph))*100, 2), "%"))
s.ph <- m.ph[m.ph==1]
tax.info.sph <- tax.info.hp[tax.info.hp$query %in% names(s.ph), ]

###Order based on score and only keep the first value (this needs a bit of a redoing as I think it makes for a very absolute arbitration)
tax.sc.sph <- tapply(1:nrow(tax.info.sph),tax.info.sph$query, function(i){
        i[ order(tax.info.sph[i, 'score'], decreasing=T) ]
    })

tax.t.sc.sph <- sapply(tax.sc.sph, function(x){ x[1] })
tax.ts.sph <- tax.info.sph[tax.t.sc.sph,]
tax.ts.sph <- tax.ts.sph[grep("chromosome|genome",tax.ts.sph$title, invert=TRUE),]
##Number of unique high quality queries before and after multiple phylum filter and chromosome/genome filter
nrow(tax.ts.sph)/length(unique(tax.info.hp$query))

l.tx.sph <- lapply(names(sal.l),function(x){
    colnames(sal.l[[x]]) <- c("query", "length", "eff.length","tpm","read.nmbr")
    merge(sal.l[[x]], tax.ts.sph[,c("query","subject", "title", "family","pident","score","phylum")],by="query")
})

list.len <- length(sample.names)
l.tx.p <- lapply(l.tx.sph, function(x){x[x$pident >= pid.thresh.2, ]})
l.tx.p <- lapply(l.tx.p, function(x){x[x$family!="NULL"|x$phylum!="NULL",]})

##Make a tpm normalization step here 
l.tx.p<- lapply(1:list.len, function(x){
    tpm.read <- (l.tx.p[[x]][,"read.nmbr"]*1000000)/sum(l.tx.p[[x]][,"read.nmbr"])
    cbind(l.tx.p[[x]], tpm.read)
})

l.tx.p <- lapply(l.tx.p, function(x){x[x$read.nmbr>=read.nmbr.thr,]})
##Remove groups that are connected to model organisms that weren't already removed prior
l.tx.p <- lapply(l.tx.p, function(x){x[x$phylum!="Basidiomycota"&
                                       x$family!="Schistosomatidae"&
                                       x$family!="Rhabditidae"&
                                       x$family!="Aplysiidae"&
                                       x$family!="Malasseziaceae",]})

###Clean up some messy genera entries
l.tx.p <- lapply(l.tx.p, function(y){
    a <- cbind(y, species=unlist(lapply(strsplit(y[, "title"], " "), function(x){ 
        if(x[1]%in% c("PREDICTED:", "TPA:")){
            paste(x[2], x[3])}
        else{paste(x[1], x[2])}})),
        genus=unlist(lapply(strsplit(y[, "title"], " "), function(x){ 
            if(x[1]%in% c("PREDICTED:", "TPA:")){
                paste(x[2])}
            else{paste(x[1])}})))
    a
})

### Aggregate read numbers based on phylum
agg.ph <- lapply(1:list.len,function(x){
    agg.ph <- aggregate(l.tx.p[[x]][,count.mode], by=list(Category=l.tx.p[[x]]$phylum), FUN=sum)
    colnames(agg.ph) <- c("phylum", "ph.sum")
    agg.ph
                                        #        l.t[l.t$ph.sum>10,]
})

### Aggregate read numbers based on famly
agg.f <- lapply(1:list.len,function(x){
    agg.f <- aggregate(l.tx.p[[x]][,count.mode], by=list(Category=l.tx.p[[x]]$family), FUN=sum)
    colnames(agg.f) <- c("family", "fam.sum")
    agg.f
})

### Aggregate read numbers based on genus
agg.gen <- lapply(1:list.len,function(x){
    agg.g <- aggregate(l.tx.p[[x]][,count.mode], by=list(Category=l.tx.p[[x]]$genus), FUN=sum)
    colnames(agg.g) <- c("genus", "gen.sum")
    agg.g
})

##Make a family abundance matrix
all.fams <- unique(sort(unlist(lapply(agg.f, function(x){x$family}))))
fam.matrix <- matrix(0, nrow=length(all.fams), ncol=length(sample.names))
rownames(fam.matrix) <- all.fams
colnames(fam.matrix) <- c(sample.names)

for(i in 1:length( agg.f)){
    for(j in all.fams){
        if(which(rownames(fam.matrix)==j)>0){
            a <- as.numeric( agg.f[[i]][agg.f[[i]]$family==j,"fam.sum"])
            fam.matrix[which(rownames(fam.matrix)==j),i]  <- ifelse(length(a)>0,a,0)
        }}
}

###Make a phylum abundance matrix
all.ph <- unique(sort(unlist(lapply(agg.ph, function(x){x$phylum}))))
ph.matrix <- matrix(0, nrow=length(all.ph), ncol=length(sample.names))
rownames(ph.matrix) <- all.ph
colnames(ph.matrix) <- c(sample.names)

for(i in 1:length( agg.ph)){
    for(j in all.ph){
        if(which(rownames(ph.matrix)==j)>0){
            a <- as.numeric( agg.ph[[i]][agg.ph[[i]]$phylum==j,"ph.sum"])
            ph.matrix[which(rownames(ph.matrix)==j),i]  <- ifelse(length(a)>0,a,0)
        }}
}

###Make a genus abundance matrix
all.gen <- unique(sort(unlist(lapply(agg.gen, function(x){x$genus}))))
gen.matrix <- matrix(0, nrow=length(all.gen), ncol=length(sample.names))
rownames(gen.matrix) <- all.gen
colnames(gen.matrix) <- c(sample.names)

for(i in 1:length( agg.gen)){
    for(j in all.gen){
        if(which(rownames(gen.matrix)==j)>0){
            a <- as.numeric( agg.gen[[i]][agg.gen[[i]]$genus==j,"gen.sum"])
            gen.matrix[which(rownames(gen.matrix)==j),i]  <- ifelse(length(a)>0,a,0)
        }}
}

##Remove any genera that are either not correctly inputed, or undefined (such as families, "dae")
dim(gen.matrix)
gen.matrix <- gen.matrix[rownames(gen.matrix)[!grepl("\\.",rownames(gen.matrix))],]
dim(gen.matrix)
gen.matrix <- gen.matrix[rownames(gen.matrix)[!grepl("dae$",rownames(gen.matrix))],]
dim(gen.matrix)

## Remove all the lower core samples and keep only the upper one using core. l
l.tx.p <- lapply(core.l, function(x){l.tx.p[[x]]})

fam.matrix <- fam.matrix[,core.l]
colnames(fam.matrix) <- eco.s$Sample
###Make sure that no tax entries remain empty as some samples might only had some taxa in the lower core layers
fam.matrix <- fam.matrix[rowSums(fam.matrix)>0,]

ph.matrix <- ph.matrix[,core.l]
colnames(ph.matrix) <- eco.s$Sample
ph.matrix <- ph.matrix[rowSums(ph.matrix)>0,]

gen.matrix <- gen.matrix[,core.l]
colnames(gen.matrix) <- eco.s$Sample
gen.matrix <- gen.matrix[rowSums(gen.matrix)>0,]

a <- fact.indx(1:14, eco.s)
ph.m.f <- do.call(cbind,lapply(a, function(x){rowSums(ph.matrix[,x])}))
fam.m.f <- do.call(cbind,lapply(a, function(x){rowSums(fam.matrix[,x])}))
gen.m.f <- do.call(cbind,lapply(a, function(x){rowSums(gen.matrix[,x])}))

####Create a taxonomy reference table of 
ph.fam <- do.call(rbind, l.tx.p)[,c("family","phylum", "genus")]
ph.fam <- ph.fam[!duplicated(ph.fam),]
