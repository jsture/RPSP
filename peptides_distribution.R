## Script to generate peptide distribution plots and circos plot 
# Author: Niels Banhos Danneskiold-Samsoee

setwd("./PeptidePredictor")

library(readxl)
library(openxlsx)
library(dplyr)
library(purrr)
library(networkD3)
library(tidyr)
library(dittoSeq)
library(circlize)
library(stringr)
library(scales)
library(MASS)
library(ggpubr)
library(ggridges)
library(ggplot2)
library(tidyverse)
library(plyr)

# Loading proteins, peptides and tissues
peptides.xls <-lapply(excel_sheets("./matches_35.xlsx"), 
                      read_excel, path = "matches_35.xlsx")
tissues <- getSheetNames("matches_35.xlsx")
names(peptides.xls) <- tissues
head(peptides.xls)
peptides <- do.call(rbind.data.frame, peptides.xls)
head(peptides) 

# Counting number of proteins per tissue
tissue.times <- peptides.xls %>%
  map_int(nrow)
# Set tissue variable
peptides$tissue <- rep(names(tissue.times),times=tissue.times)
head(peptides)
table(peptides$id)
length(peptides$id)
length(unique(peptides$id))

# Get peptides
pept <- str_split(peptides$`peptide sequences`,pattern=",")
head(pept)
names(pept) <- peptides$tissue 
protein.times <- sapply(pept, function(x) length(x))
pept <- gsub("\n","",unlist(pept))
head(pept)
pept_data <- data.frame(tissue=names(pept),peptide=unname(pept),protein=rep(peptides$id,times=protein.times))
head(pept_data)
pept_data <- pept_data %>% filter(nchar(peptide)>0)
pept_data$protein
pept_data$tissue <- gsub("\\d+$","",pept_data$tissue)
head(pept_data)
pept_data$peplength <- nchar(pept_data$peptide)
head(pept_data)

# remove trailing lysine and arginines 
pept_data[grep("[KR]+$",pept_data$peptide),]
pept_data$peptide <- gsub("[KR]+$","",pept_data$peptide)

# get unique peptides
dim(pept_data)[1] # number of peptides generated in different tissues
dim(pept_data %>% distinct(protein,peptide, .keep_all=TRUE)) # number of peptides assuming same modifications across tissues
dim(pept_data %>% distinct(peptide, .keep_all=TRUE)) # number of unique peptide sequences
pept_data %>% distinct(protein,peptide, .keep_all=TRUE) %>% group_by(peptide) %>% 
  filter(n()>1)

# Get summary statistics on predicted peptide lengths
quantile(pept_data$peplength, probs=c(.25, .75))
pept_data %>% distinct(protein,peptide, .keep_all=TRUE) %>% 
  summarize(mean=mean(peplength),
            sd=sd(peplength),
            median=median(peplength))

# Get summary statistics on predicted peptide lengths by tissue
quantile(pept_data$peplength, probs=c(.25, .75))
pept_data %>% distinct(protein,peptide, .keep_all=TRUE) %>% 
  group_by(tissue) %>%
  summarize(mean=mean(peplength),
            sd=sd(peplength),
            median=median(peplength)) %>% print(n=36)

# Get all peptides, and peptide lengths matched to parent protein 
dim(pept_data %>% distinct(protein,peptide, .keep_all=TRUE))
head(pept_data %>% distinct(protein,peptide, .keep_all=TRUE))
pept_data %>% distinct(protein,peptide, .keep_all=TRUE)

# Plot length distribution of peptides
ggplot(data = pept_data %>% distinct(protein,peptide, .keep_all=TRUE), aes(x =peplength)) +
  geom_histogram(binwidth = 2, color="#E69F00",fill="#E69F00",linewidth=0.1) + 
  xlab("peptide length (AA)") +
  theme_bw()+
  theme(text = element_text(size=9),
        axis.text = element_text(size=9))
#ggsave("./peptide_dist_linear.emf",width=2.5,height=2)

# Load and filter known peptides from Foster et al. 
# https://www.cell.com/cms/10.1016/j.cell.2019.10.010/attachment/dd67b8ea-768d-48d2-a00e-16e5cbb983d3/mmc1.xlsx
known_GPCR_ligands <- read_xlsx("./foster_et_al_Suppl_Table_1.xlsx","HumanGPCRLigands")
head(known_GPCR_ligands)
known_GPCR_ligands <- known_GPCR_ligands %>% filter(type=="Peptide") # Remove small molecules
known_GPCR_ligands <- known_GPCR_ligands %>% filter(!str_detect(Family, "Chemokine")) # Remove chemokines
known_GPCR_ligands <- known_GPCR_ligands %>% filter(PreAA!="" | AfterAA!="") # remove peptides without annotated start or ending
dim(known_GPCR_ligands)

# Remove duplicate peptides
known_GPCR_ligands %>% group_by(Uniprot,PepStartLoc,LigandLength) %>% 
  filter(n()>1)
dim(known_GPCR_ligands %>% distinct(Uniprot,PepStartLoc,LigandLength, .keep_all=TRUE))

# Combine predicted and known peptides
combined_pep <-  data.frame(peplength=pept_data %>% distinct(protein,peptide, .keep_all=TRUE) %>% pull(peplength),type="predicted")
combined_pep <- rbind(combined_pep,data.frame(peplength=known_GPCR_ligands %>% distinct(Uniprot,PepStartLoc,LigandLength, .keep_all=TRUE) %>% pull(LigandLength),type="known"))
tail(combined_pep)

# Investigate distribution of peptide lengths
ggqqplot(log(pept_data$peplength))
ggqqplot(sqrt(pept_data$peplength))

# Peptide length is not log distributed
shapiro.test(log(pept_data %>% distinct(protein,peptide, .keep_all=TRUE) %>% pull(peplength)))

# Plot length distribution of known and predicted peptides
combined_pep$type <- as.factor(combined_pep$type)
type_levels <- levels(combined_pep$type)
type_levels_reversed <- rev(type_levels)
combined_pep$type <- factor(combined_pep$type, levels = type_levels_reversed)
ggplot(data = combined_pep, aes(x =log(peplength),color=type,fill=type)) +
  geom_histogram(aes(y=..count..),binwidth = 0.1,linewidth=0.1, alpha=0.6,position='identity') + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  xlab("log(peptide length) (AA)") +
  ylab("# peptides") +
  #coord_trans(x='log') +
  stat_function(fun = function(x) 270 * dnorm(x, mean = 4.2, sd = 1), color = "#56B4E9", linewidth = 0.3) +
  annotate("rect", xmin=log(5), xmax=log(25), ymin=-10, ymax=0, alpha=0.8, fill="#3CBB75FF") + 
  theme(text = element_text(size=12),
        axis.text = element_text(size=12)) +
  theme_bw()
#ggsave("./peptide_dist.emf",width=4,height=2.5)
#ggsave("./peptide_dist_count.pdf",width=4,height=2.5)

ggplot(data = combined_pep, aes(x = peplength, color = type, fill = type)) +
  geom_histogram(aes(y = ..count..), binwidth = 0.1, linewidth = 0.1, alpha = 0.6, position = 'identity') + 
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  scale_x_continuous(trans = "log10", breaks = c(5, 10, 25, 50, 100,200,500,1000)) +  # Set breaks as desired
  xlab("Peptide length (AA)") +
  ylab("# peptides") +
  stat_function(fun = function(x) 550 * dnorm(log(x), mean = 4.2, sd = 1), color = "#56B4E9", linewidth = 0.3) +
  annotate("rect", xmin = 5, xmax = 25, ymin = -10, ymax = 0, alpha = 0.8, fill = "#3CBB75FF") + 
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  theme_bw()
#ggsave("./peptide_dist_count_linear_x.pdf",width=4,height=2.5)

# Plot length distribution of known and predicted peptides
type_levels <- levels(combined_pep$type)
type_levels_reversed <- rev(type_levels)
combined_pep$type <- factor(combined_pep$type, levels = type_levels_reversed)
ggplot(data = combined_pep, aes(x =log(peplength),color=type,fill=type)) +
  geom_histogram(aes(y=2*(..density..)/sum(..density..)),binwidth = 0.1,linewidth=0.1, alpha=0.5,position='identity') + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  xlab("log(peptide length) (AA)") +
  ylab("Percent of peptides\nfor each peptide type") +
  scale_y_continuous(labels=percent_format()) +
  stat_function(fun = function(x) 0.09 * dnorm(x, mean = 4.2, sd = 1), color = "#56B4E9", linewidth = 0.5) +
  annotate("rect", xmin=log(5), xmax=log(25), ymin=-0.005, ymax=0, alpha=0.8, fill="#3CBB75FF") + 
  theme(text = element_text(size=12),
        axis.text = element_text(size=12)) +
  theme_bw()
#ggsave("./peptide_dist.emf",width=4,height=2.5)
#ggsave("./peptide_dist_percent.pdf",width=4,height=2.5)

# Annotate peptides within lengths used for screen
pept_data_distinct <- pept_data %>% distinct(protein,peptide, .keep_all=TRUE) %>% 
  mutate(inrange=if_else(peplength>=5 & peplength<=25,"yes","no")) 

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

# Make pie diagram illustrating percent within lengths included in screen
inrange <- pept_data_distinct %>% 
            group_by(inrange) %>% 
            tally() %>% 
            mutate(percent = n / sum(n) * 100) %>%
            ungroup()
bp<- ggplot(inrange, aes(x="", y=percent, fill=inrange))+
  geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0)
pie
pie + scale_fill_manual(values=c("#E69F00", "#3CBB75FF")) + 
  blank_theme +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = percent/2 + c(0, cumsum(percent)[-length(percent)]), 
                label = percent(percent/100)), size=5)
#ggsave("./peptides_pie_included.emf",width=2.5,height=2.5)

# Count number of peptides in screen range
inrange <- pept_data_distinct %>% 
  group_by(inrange) %>% 
  tally() %>% 
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()
tail(inrange)
min(inrange$n)
max(inrange$n)
inrange

# Plot general length distribtuon of peptides across tissues
ggplot(pept_data, aes(x = peplength, y = factor(tissue))) +
  geom_density_ridges(alpha = 0.5,jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0,yoffset = .2),
                      point_shape = '|', point_size = 1.5, point_alpha = 2, alpha = 0,color="gray30") + 
  xlab("peptide length (AA)") +
  xlim(c(0,1500)) +
  theme_bw() +
  theme(axis.title.y = element_blank())
ggplot(pept_data, aes(x = log(peplength), y = factor(tissue))) +
  geom_density_ridges(alpha = 0.5) + 
  xlab("peptide length log(AA)") +
  theme_bw() +
  theme(axis.title.y = element_blank())

pept_data <- pept_data  %>% 
  mutate(inrange=if_else(peplength>=5 & peplength<=25,"yes","no"))

inrange <- pept_data %>% 
  group_by(tissue,inrange) %>% 
  tally() %>% 
  mutate(percent = n / sum(n) * 100) %>%
  ungroup() %>%
  filter(inrange=="yes") %>%
  arrange(desc(percent))

pept_data <- left_join(pept_data,inrange %>% dplyr::select(tissue,percent),by="tissue")

# No difference in peptide length between tissues within peptides < 24 & > 5 amino acids in length
inrange_all <- pept_data %>% filter(peplength<=25 & peplength>=5) %>% arrange(desc(percent))
res.aov <- aov(peplength ~ tissue, data = inrange_all)
summary(res.aov)
ggplot(inrange_all) +
  geom_density_ridges(aes(x = peplength, y = factor(tissue)),alpha = 0.5,jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0,yoffset = .2),
                      point_shape = '-', point_size = 2, point_alpha = 2, alpha = 0,color="gray30",fill="#3CBB75FF") + 
  coord_flip() +
  xlab("peptide length (AA)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))
#ggsave("./peptide_dist_inrange_pertissue.pdf",width=8,height=3)

# Plot length distribution of peptides < 24 & > 5 amino acids in length by tissue
inrange_all <- inrange_all %>% mutate(tissue = fct_relevel(tissue, levels = inrange$tissue))
head(inrange_all)
inrange_all <- inrange_all %>% arrange(desc(percent))
inrange_all$tissue <- factor(inrange_all$tissue, levels = inrange$tissue, ordered = TRUE)
ggplot(inrange_all) +
  geom_density_ridges(data=inrange_all,aes(x = peplength, y = tissue),alpha = 0.5,jittered_points = TRUE,
                    position = position_points_jitter(width = 0.05, height = 0,yoffset = .2),
                    point_shape = '-', point_size = 2, point_alpha = 2, alpha = 0,color="gray30",fill="#3CBB75FF") + 
  coord_flip() +
  xlab("peptide length (AA)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))
#ggsave("./peptide_dist_inrange_pertissue.pdf",width=8,height=3)

write.table(inrange_all,"./peptides_5_to_25AA.tsv",sep="\t")

# Plot percentage of total peptides in < 24 & > 5 amino acid length range by tissue
inrange$tissue <- factor(inrange$tissue, levels = inrange$tissue, ordered = TRUE)
ggplot(data=inrange,aes(tissue,percent)) +
  geom_col(color="gray30",fill="#3CBB75FF",alpha=0.5) + 
  theme_bw() + 
  ylab("Percent of predicted\npeptides ≥ 5 & ≤ 25 (AA)") +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
#ggsave("./peptide_inrange_pecent_pertissue.pdf",width=8,height=3)

pept_data <- pept_data %>% mutate(tissue = fct_relevel(tissue, levels = pept_data$tissue))
head(pept_data)
pept_data <- pept_data %>% arrange(desc(percent))
pept_data$tissue <- factor(pept_data$tissue, levels = pept_data$tissue, ordered = TRUE)
ggplot(pept_data) +
  geom_density_ridges(data=pept_data,aes(x = peplength, y = tissue),alpha = 0.5,jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0,yoffset = .2),
                      point_shape = '-', point_size = 2, point_alpha = 2, alpha = 0,color="gray30",fill="#3CBB75FF") + 
  #geom_density_ridges(aes(x = log(peplength), y = Species) +
  coord_flip() +
  xlim(c(0,100))+
  xlab("peptide length (AA)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))
#ggsave("./peptide_dist_inrange_pertissue.pdf",width=8,height=3)

# List percentage of peptides included by size per tissue
pept_data %>% group_by(tissue) %>%
  mutate(small=peplength<=25 & peplength>=5) %>%
  group_by(tissue,small) %>%
  count() %>% group_by(tissue) %>%
  mutate(fraction = 100*n / sum(n)) %>%
  ungroup() %>% filter(small==TRUE) %>%
  arrange(fraction) %>% print(n=70)
pept_data %>% group_by(tissue) %>%
  filter(peplength<=25& peplength>=5) %>%
  group_by(tissue) %>%
  count() %>% group_by(tissue) %>%
  mutate(fraction = 100*n / sum(n)) %>%
  ungroup() %>% 
  arrange(fraction) %>% print(n=70)
  
pepdist <- pept_data %>% group_by(tissue) %>%
  mutate(small=peplength<=25) %>%
  group_by(tissue,small) %>%
  count() %>% group_by(tissue) %>%
  mutate(fraction = 100*n / sum(n)) %>%
  ungroup() %>% filter(small==TRUE) %>%
  arrange(fraction)

### Make circus plot

## Make links
peptides['num_peptides'] <- peptides[, c(4)] + 1 #number of peptides will be number of cleavage sites+1
peptides <- peptides[order(peptides$num_peptides),] #sort dataframe by number of peptides
links <- peptides %>%
  dplyr::arrange(tissue) %>% dplyr::group_by(tissue) %>%
  dplyr::mutate(idnr = cumsum(id != lag(id, default=""))) %>%
  dplyr::arrange(id) %>% dplyr::group_by(id) %>% 
  dplyr::summarise(id, tissue, `number of cleavage sites`,idnr) %>%
  dplyr::mutate(target=tissue, idnr1=idnr) %>%
  expand(tissue,target) %>%
  dplyr::filter(tissue!=target) %>%
  dplyr::filter(!duplicated(paste0(pmax(tissue, target), pmin(tissue, target))))

# Assign id number to each prohormone by number of peptides
idnr <- peptides %>%
  dplyr::arrange(tissue,num_peptides) %>% dplyr::group_by(tissue) %>% #replace with id
  dplyr::mutate(idnr = cumsum(id != lag(id, default=""))) %>%
  dplyr::arrange(id) %>% group_by(id) %>% 
  dplyr::summarise(tissue, `number of cleavage sites`,idnr)


# Populate links from each prohormone's idnr 
links$idtissue <- 0
links$idtarget <- 0
for (i in 1:length(links$id)){
  links$idtissue[i] <-
    idnr[idnr$id == links[i,]$id & idnr$tissue==links[i,]$tissue,]$idnr
  links$idtarget[i] <-
    idnr[idnr$id == links[i,]$id & idnr$tissue==links[i,]$target,]$idnr
}
links <- links %>% arrange(tissue,target)
links %>% print(n=10)

# Order tissue
tisorder <- c("adipose tissue","breast","intestine","blood","lymphoid tissue",
              "testis","brain","placenta","cervix, uterine","endometrium 1",
              "esophagus","tongue","skin 1","lung","kidney","liver",
              "adrenal gland","pituitary gland","bone marrow",
              "ductus deferens","epididymis","fallopian tube","gallbladder",
              "heart muscle","no HPA annotation","ovary","pancreas",
              "parathyroid gland","thyroid gland","prostate","retina",
              "salivary gland","seminal vesicle","skeletal muscle",
              "smooth muscle","urinary bladder","vagina","wide")
tisorder <- tissue.times[match(tisorder,names(tissue.times))]

# If ordering in alphabetical order
tisorder <- tissue.times
total_peptides <- idnr %>% group_by(tissue) %>% 
  dplyr::summarise(Freq = sum(`number of cleavage sites`))

# Drawing circular plot
idnr <- idnr %>% arrange(tissue, idnr)
m <- data.frame(xlim1=0,xlim2=as.numeric(unname(tisorder)))
m <- as.matrix(m)
rownames(m) <- names(tisorder)
par(mar=rep(0,4))

circos.clear()
circos.par(track.margin=c(0,0))
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0.0,0.18), start.degree = 90,
           gap.degree =2.4)
circos.initialize(sectors=names(tisorder),xlim=m,
                  sector.width=unname(tisorder))
circos.trackPlotRegion(sectors = names(tisorder), y=unname(tisorder), 
                       bg.col = dittoColors()[1:38] , bg.border = "gray",
                       bg.lwd=0.5,
                       track.height = 0.2,
                       panel.fun = function(x, y) {
                         #Tissue labels
                         circos.axis(labels=get.current.sector.index(),
                                     major.at = m[grep(paste("^",get.current.sector.index(),"$", sep=""),rownames(m)),2]/2,
                                     labels.facing="clockwise", labels.cex = 1)
                         circos.axis(h = "bottom",labels=NULL, #labels=c("","10","","30","","50","","","","",""),
                                     major.at = c(0,10,20,30,40,50,60,70,80,90,100),
                                     direction ="inside",
                                     labels.facing="clockwise",
                                     major.tick.length=30)
                         #numbers 
                         circos.axis(labels=table(peptides$tissue)[get.current.sector.index()],
                                     major.at = m[grep(paste("^",get.current.sector.index(),"$", sep=""),rownames(m)),2]/2,
                                     major.tick=FALSE,
                                     labels.facing="clockwise",
                                     h = "bottom",
                                     direction ="inside",major.tick.length=35,
                                     labels.cex = 1)
                       })

#Add bar graph to each tissue in order of number of cleavage sites
circos.track(track.index=1,ylim = c(0, max(idnr$`number of cleavage sites`)+2),sectors = names(tisorder), y=unname(tisorder), panel.fun = function(x, y) {
  for(pos in seq(1, m[grep(paste("^",get.current.sector.index(),"$", sep=""),rownames(m)),2], by = 1)) {
    value=sort(idnr$`number of cleavage sites`[grep(get.current.sector.index(),idnr$tissue)])[pos]
    circos.barplot(value, pos,bar_width=0.0001)
  }
})

#draw links
#cou=1
tisnr <- 999

for(i in 1:length(links$id)) {
  if (i<(length(links$id)-1) & links$tissue[i]==links$tissue[i+1] & 
      links$target[i]==links$target[i+1] &
      links$idtissue[i]+1==links$idtissue[i+1] & 
      links$idtarget[i]+1==links$idtarget[i+1]){
    if (links$idtissue[i]<tisnr){
      tisnr <- links$idtissue[i]
      tarnr <- links$idtarget[i]
    }
  } else {
    if (tisnr==999){
      #print(cou)
      circos.link(sector.index1=links$tissue[i], links$idtissue[i], 
                  sector.index2=links$target[i], links$idtarget[i],
                  rou=0.48,h.ratio = 0.9,col="grey50")
    }
    if (tisnr!=999){
      tisnr
      tarnr
      circos.link(sector.index1=links$tissue[i], point1=c(tisnr,links$idtissue[i]), 
                  sector.index2=links$target[i], point2=c(tarnr,links$idtarget[i]),
                  rou=0.48,h.ratio = 0.9,col="grey50")
      tisnr <- 999
      tarnr <- 999
    }
  }
  #cou <- cou +1
}

#export to pdf 
#dev.copy2pdf(file=paste0("circular_plot",paste0(".pdf")), width=16, height=16)
#dev.off()
