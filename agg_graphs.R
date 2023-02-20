# Script for making graphs

options(stringsAsFactors = FALSE)
options(warn = 1)
globalTextSize <- 7

aggData_rmsd_0.5 <- read.table("IRON_DATA/combined_agg_ligand_summary_info_with_novelty_rmsd_0.5.txt", 
                               sep = "", header = F, stringsAsFactors = F)

aggData_rmsd_0.5$V21 <- gsub('.{4}$', '', aggData_rmsd_0.5$V21)

uniprot_data <- read.table("IRON_DATA/uniprot_writeout.txt",
                           sep = "\t", header = T, stringsAsFactors = F, fill = T)
uniprot_data <- uniprot_data[which(uniprot_data$composition == "4Fe-4S" |
                                       uniprot_data$composition == "3Fe-4S" |
                                       uniprot_data$composition == "2Fe-2S" |
                                       uniprot_data$composition == "Fe" |
                                       uniprot_data$composition == "Zn"),]
uniprot_data$composition <- ifelse(uniprot_data$composition == "Fe",
                                   "Fe/Zn",
                                   uniprot_data$composition)
uniprot_data$composition <- ifelse(uniprot_data$composition == "Zn",
                                   "Fe/Zn",
                                   uniprot_data$composition)


library(ggplot2)
library(viridis)
library(gridExtra)
library(cowplot)
library(scales)
#library(svglite)

integer_breaks <- function(n = 5, ...) {
    fxn <- function(x) {
        breaks <- floor(pretty(x, n, ...))
        names(breaks) <- attr(breaks, "labels")
        breaks
    }
    return(fxn)
}


ligandOrderVector <- c("4Fe4S_4CYS", 
                       "4Fe4S_3CYS",     
                       "3Fe4S_3CYS",
                       "2Fe2S_4CYS",     
                       "2Fe2S_2CYS_2HIS",
                       "2Fe2S_3CYS_1ASP",
                       "1FeZn_4CYS",     
                       "1ZINC_3CYS_1HIS",
                       "1ZINC_2CYS_2HIS",
                       "1ZINC_1CYS_3HIS",
                       "1ZINC_4HIS",     
                       "1ZINC_3HIS")

ligandOrderVectorWriteOut <- c("4Fe-4S 4CYS", 
                               "4Fe4S_3CYS",     
                               "3Fe4S_3CYS",
                               "2Fe2S_4CYS",     
                               "2Fe2S_2CYS_2HIS",
                               "2Fe2S_3CYS_1ASP",
                               "1FeZn_4CYS",     
                               "1ZINC_3CYS_1HIS",
                               "1ZINC_2CYS_2HIS",
                               "1ZINC_1CYS_3HIS",
                               "1ZINC_4HIS",     
                               "1ZINC_3HIS")

ligandCommonComparisonVector <- c("4Fe-4S",
                                  "4Fe-4S",
                                  "3Fe-4S",
                                  "2Fe-2S",
                                  "2Fe-2S",
                                  "2Fe-2S",
                                  "Fe/Zn",
                                  "Fe/Zn",
                                  "Fe/Zn",
                                  "Fe/Zn",
                                  "Fe/Zn",
                                  "Fe/Zn")

ligandUniprotName <- c("4Fe4S_4CYS", 
                       "4Fe4S_3CYS",    
                       "2Fe2S_4CYS",    
                       "2Fe2S_2CYS_2HIS",
                       "2Fe2S_3CYS_1ASP",
                       "1FeZn_4CYS",    
                       "1FeZn_4CYS",
                       "3Fe4S_3CYS",
                       "1ZINC_3CYS_1HIS",
                       "1ZINC_2CYS_2HIS",
                       "1ZINC_1CYS_3HIS",
                       "1ZINC_4HIS",    
                       "1ZINC_3HIS")


ligandCommonOrderVector <- c("4Fe-4S",
                             "3Fe-4S",
                             "2Fe-2S",
                             "Zn")


orgOrderVector <- c("ORG1",
                    "ORG2",
                    "ORG3",
                    "ORG4",
                    "ORG5",
                    "ORG6",
                    "ORG7",
                    "ORG8",
                    "ORG9",
                    "ORG10",
                    "ORG11",
                    "ORG12",
                    "ORG13",
                    "ORG14",
                    "ORG15",
                    "ORG16",
                    "ORG17",
                    "ORG18",
                    "ORG19",
                    "ORG20",
                    "ORG21")

orgNameVector <- c("A. thaliana",          
                   "C. elegans",        
                   "C. albicans",              
                   "D. rerio",                   
                   "D. discoideum",      
                   "D. melanogaster",       
                   "E. coli",              
                   "G. max",                   
                   "H. sapiens",                  
                   "L. infantum",          
                   "M. jannaschii",
                   "M. musculus",                 
                   "M. tuberculosis",   
                   "O. sativa",                 
                   "P. falciparum",        
                   "R. norvegicus",            
                   "S. cerevisiae",     
                   "S. pombe",    
                   "S. aureus",        
                   "T. cruzi",            
                   "Z. mays")

orgNameOrderVector <- c("E. coli",
                        "S. aureus",
                        "M. tuberculosis",
                        "M. jannaschii",
                        "A. thaliana",
                        "G. max",
                        "O. sativa",
                        "Z. mays",
                        "S. pombe",
                        "S. cerevisiae",
                        "C. albicans",
                        "C. elegans",
                        "D. melanogaster",
                        "D. rerio",
                        "H. sapiens",
                        "M. musculus",
                        "R. norvegicus",
                        "L. infantum",
                        "T. cruzi",
                        "P. falciparum",
                        "D. discoideum")


update_names <- function(inData, inOrgNames, outOrgNames) {
    
    for (i in 1:length(inOrgNames)) {
        print(i)
        for (j in 1:length(inData)) {
            if (inData[j] == inOrgNames[i]) {
                inData[j] <- outOrgNames[i]
            }
        }
    }
    return(inData)
    
}

add_general_classification <- function(inData, inLigandNames, outCommonNames) {
    V26 <- character(length = nrow(inData))
    for (i in 1:length(inLigandNames)) {
        print(i)
        for (j in 1:nrow(inData)) {
            if (inData$V21[j] == inLigandNames[i]) {
                V26[j] <- outCommonNames[i]
            }
        }
    }
    inData <- cbind.data.frame(inData, V26)
    return(inData)
}

add_ligand_names_to_uniprot <- function(inData,
                                        inVector) {
    ligandName <- character(length = nrow(inData))
    
    for (i in 1:nrow(inData)) {
        iVal <- inData$matchIndex[i]
        if (iVal == 0) {
            ligandName[i] <- "NOTHING"
        }
        if (iVal == -1) {
            ligandName[i] <- "REGION"
        }
        if (iVal > 0) {
            ligandName[i] <- inVector[as.integer(iVal)]
        }
        
    }
    
    inData$ligandName <- ligandName
    return(inData)
    
}


aggData_rmsd_0.5$V25 <- update_names(aggData_rmsd_0.5$V25, orgOrderVector, orgNameVector)
uniprot_data$composition <- update_names(uniprot_data$composition, ligandOrderVector, ligandCommonComparisonVector)
aggData_rmsd_0.5 <- add_general_classification(aggData_rmsd_0.5, ligandOrderVector, ligandCommonComparisonVector)
uniprot_data[,7] <- update_names(uniprot_data[,7], orgOrderVector, orgNameVector)
uniprot_data <- add_ligand_names_to_uniprot(uniprot_data, ligandUniprotName)
uniprot_data <- uniprot_data[which(uniprot_data$ligandName != "NOTHING"),]
update_names(aggData_rmsd_0.5$V21, orgOrderVector, orgNameVector)


p1 <- ggplot(data = aggData_rmsd_0.5, aes(x = factor(V21, level = ligandOrderVector))) +
    geom_bar(aes(fill = factor(V21, level = orgNameOrderVector)), color = "black") +
    facet_wrap(~factor(V25, levels = orgOrderVector), scales = "free", nrow = 3) +
    theme_bw() +
    xlab("") +
    ylab("Counts") +
    ggtitle("RMSD = 0.5") +
    scale_fill_discrete(drop = FALSE) +
    guides(fill=guide_legend(title="Ligand")) +
    theme(axis.text.x=element_blank())

ggsave(filename = "../IRON_DATA/barchart_of_ligand_fequencies_all_orgs.bmp",
       plot = p1,
       units = "in", dpi = 300, height = 15, width = 20)

p2 <- ggplot(data = aggData_rmsd_0.75, aes(x = factor(V21, level = ligandOrderVector))) +
    geom_bar(aes(fill = factor(V21, level = ligandOrderVector)), color = "black") +
    facet_wrap(~factor(V25, levels = orgOrderVector), scales = "free", nrow = 3) +
    theme_bw() +
    xlab("") +
    ylab("Counts") +
    ggtitle("RMSD = 0.75") +
    scale_fill_discrete(drop = FALSE) +
    guides(fill=guide_legend(title="Ligand")) +
    theme(axis.text.x=element_blank())



# Stacked bar chart for ligand type
p1 <- ggplot(data=aggData_rmsd_0.5, aes(factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V26, level = ligandCommonOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("a") +
    #guides(fill=guide_legend(title="Ligand")) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    scale_fill_brewer(palette="Spectral") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

p2 <- ggplot(data=uniprot_data, aes(factor(uniprot_data[,7], levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(composition, level = ligandCommonOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("b") +
    guides(fill=guide_legend(title="Ligand")) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    scale_fill_brewer(palette="Spectral") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

p3 <- ggplot(data=aggData_rmsd_0.5, aes(factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("c") +
    guides(fill=guide_legend(title="Ligand")) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

tmpData <- uniprot_data[grep("-", uniprot_data$resid),]
p4 <- ggplot(data=subset(tmpData, as.integer(matchIndex) > 0), aes(factor(rep.iOrg..nrow.iDF.., levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(ligandName, level = ligandOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("d") +
    guides(fill=guide_legend(title="Ligand")) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

pOut <- grid.arrange(p1,p2,p3,p4, nrow = 2, ncol = 2)
ggsave(filename = "../IRON_DATA/stacked_barchart_comparison.bmp",
       plot = pOut,
       units = "in", dpi = 300, height = 10, width = 12)
ggsave(filename = "../IRON_DATA/stacked_barchart_comparison.svg",
       plot = pOut,
       units = "in", dpi = 300, height = 10, width = 12)


# Stacked bar chart for ligand type
# Figure 5
# Revamped panel c
p1 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70), aes(factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V26, level = ligandCommonOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts of identified ligands (this work)") +
    #guides(fill=guide_legend(title="Ligand")) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    scale_fill_brewer(palette="Spectral") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=8)) +
    theme(legend.position = "none")

p2 <- ggplot(data=uniprot_data, aes(factor(uniprot_data[,7], levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(composition, level = ligandCommonOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts of ligands found in UniProt") +
    guides(fill=guide_legend(title="Ligand")) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    scale_fill_brewer(palette="Spectral") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=8)) +
    theme(legend.position = "none")

p3 <- ggplot(data=subset(aggData_rmsd_0.5, V23 == "KNOWN_UNIPROT_MB" & V12 > 70), aes(factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts of identified ligands (this work)\n with known UniProt metal binding sites") +
    guides(fill=guide_legend(title="Ligand")) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=8)) +
    theme(legend.position = "none")

tmpData <- uniprot_data[grep("-", uniprot_data$resid),]
p4 <- ggplot(data=subset(tmpData, as.integer(matchIndex) > 0), aes(factor(rep.iOrg..nrow.iDF.., levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(ligandName, level = ligandOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts of ligands found in\n UniProt with known metal binding sites") +
    guides(fill=guide_legend(title="Ligand")) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=8)) +
    theme(legend.position = "none")

pOut <- grid.arrange(p1,p2,p3,p4, nrow = 2, ncol = 2)
ggsave(filename = "../IRON_DATA/stacked_barchart_comparison_2.bmp",
       plot = pOut,
       units = "in", dpi = 300, height = 8, width = 8)




# Barcharts of: x-axis - organism, y-axis stacked counts
# Main text

p1 <- ggplot(data=aggData_rmsd_0.5, aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="stack") +
    ylab("Counts") +
    xlab("Organism") +
    ggtitle("Absolute counts, pLDDT > 0") +
    ylim(c(0, 15000)) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

p2 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="stack") +
    ylab("Counts") +
    xlab("Organism") +
    ggtitle("Absolute counts, pLDDT > 70") +
    ylim(c(0, 15000)) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

p3 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 90), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="stack") +
    ylab("Counts") +
    xlab("Organism") +
    ggtitle("Absolute counts, pLDDT > 90") +
    ylim(c(0, 15000)) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")
pOut1 <- grid.arrange(p1,p2,p3, nrow = 1, ncol = 3)
ggsave(filename = "../IRON_DATA/stacked_barchart_abs_counts_by_bfac.bmp",
       plot = pOut1,
       units = "in", dpi = 300, height = 5, width = 12)
ggsave(filename = "../IRON_DATA/stacked_barchart_abs_counts_by_bfac.svg",
       plot = pOut1,
       units = "in", dpi = 300, height = 5, width = 12)



p1 <- ggplot(data=aggData_rmsd_0.5, aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="fill", ) +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts, pLDDT > 0") +
    #ylim(c(0, 15000)) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

p2 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts, pLDDT > 70") +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

p3 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 90), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts, pLDDT > 90") +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")
pOut2 <- grid.arrange(p1,p2,p3, nrow = 1, ncol = 3)
ggsave(filename = "../IRON_DATA/stacked_barchart_rel_counts_by_bfac.bmp",
       plot = pOut2,
       units = "in", dpi = 300, height = 5, width = 12)
ggsave(filename = "../IRON_DATA/stacked_barchart_rel_counts_by_bfac.svg",
       plot = pOut2,
       units = "in", dpi = 300, height = 5, width = 12)


pOut3 <- grid.arrange(pOut1, pOut2, nrow = 2, ncol = 1)
ggsave(filename = "../IRON_DATA/stacked_barchart_abs_and_rel_counts_by_bfac.bmp",
       plot = pOut3,
       units = "in", dpi = 300, height = 8, width = 9)
# ggsave(filename = "../IRON_DATA/stacked_barchart_abs_and_rel_counts_by_bfac.bmp",
#        plot = pOut3,
#        units = "in", dpi = 300, height = 6, width = 6.5)


# Stacked bar charts of only truely novel hits
# Main text
p1 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70 & V23 == "NONE"), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="stack") +
    ylab("Counts") +
    xlab("Organism") +
    ggtitle("Absolute counts of identified ligands (this work)\n not found in UniProt") +
    ylim(c(0, 8500)) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

p2 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70 & V23 == "NONE"), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts of identified ligands (this work)\n not found in UniProt") +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")
pOut1 <- grid.arrange(p1, p2, nrow = 1, ncol = 2)
ggsave(filename = "../IRON_DATA/stacked_barchart_abs_and_rel_counts_by_bfac_novel_by_uniprot.bmp",
       plot = pOut1,
       units = "in", dpi = 300, height = 6, width = 9)



p1 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70 & V23 == "NONE" & V22 == "COMPLETELY_NEW"), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="stack") +
    ylab("Counts") +
    xlab("Organism") +
    ggtitle("Absolute counts of identified ligands (this work)\n not found in UniProt or by hhsearch") +
    ylim(c(0, 8500)) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

p2 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70 & V23 == "NONE" & V22 == "COMPLETELY_NEW"), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts of identified ligands (this work)\n not found in UniProt or by hhsearch") +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")
pOut2 <- grid.arrange(p1, p2, nrow = 1, ncol = 2)
ggsave(filename = "../IRON_DATA/stacked_barchart_abs_and_rel_counts_by_bfac_novel.bmp",
       plot = pOut2,
       units = "in", dpi = 300, height = 6, width = 9)

pOut3 <- grid.arrange(pOut1, pOut2, nrow = 2, ncol = 1)
ggsave(filename = "../IRON_DATA/stacked_barchart_abs_and_rel_counts_by_bfac_novel_both.bmp",
       plot = pOut3,
       units = "in", dpi = 300, height = 7.5, width = 7.5)

# Reorder the data to allow for normalization
# divideByVector numbers are taken directly from the Alphafold website on # of predicted structures
divideByVector <- c(27434,
                    19694,
                    5974,
                    24664,
                    12622,
                    13458,
                    4363,
                    55799,
                    20505, # was 21457
                    7924,
                    1773,
                    21615,
                    3988,
                    43649,
                    5187,
                    21272,
                    6040,
                    5128,
                    2888,
                    19036,
                    39199)
normalize_counts_df <- function(inData, 
                                pLDDTCutoff, 
                                orgVector,
                                ligandVector,
                                divideByVector) {
    
    myLigand <- character()
    myOrg <- character()
    myCountsRaw1 <- integer()
    myCountsRaw2 <- integer()
    myCountsAdjusted1 <- numeric()
    myCountsAdjusted2 <- numeric()
    
    myPosition <- 1
    for (i in 1:length(orgVector)) {
        iOrg <- orgVector[i]
        iDivide <- divideByVector[i]
        for (j in 1:length(ligandVector)) {
            jLigand <- ligandVector[j]
            jIndex <- which(inData$V25 == iOrg & 
                                inData$V12 > pLDDTCutoff &
                                inData$V21 == jLigand)
            
            myLigand[myPosition] <- jLigand
            myOrg[myPosition] <- iOrg
            myCountsRaw1[myPosition] <- length(jIndex)
            myCountsAdjusted1[myPosition] <- myCountsRaw1[myPosition] / iDivide
            jIndex <- unique(inData[jIndex,18])
            myCountsRaw2[myPosition] <- length(jIndex)
            myCountsAdjusted2[myPosition] <- myCountsRaw2[myPosition] / iDivide
            myPosition <- myPosition + 1
        }
    }
    
    return(data.frame(myLigand,myOrg, myCountsRaw1, myCountsAdjusted1,myCountsRaw2, myCountsAdjusted2))
}
normDF <- normalize_counts_df(aggData_rmsd_0.5, 70, orgNameVector, ligandOrderVector, divideByVector)


# Normalized
# SI
p1 <- ggplot(data=normDF, aes(x = factor(myOrg, levels = orgNameOrderVector),
                              y = myCountsAdjusted1)) +
    geom_bar(aes(fill=factor(myLigand, level = ligandOrderVector)), position="stack", stat="identity") +
    ylab("Number of ligand binding sites /\n number of proteins in proteome") +
    xlab("Organism") +
    ggtitle("Identified number of ligand binding\n sites per protein, pLDDT > 70") +
    ylim(c(0, 0.55)) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=7)) +
    theme(legend.position = "none")
p2 <- ggplot(data=normDF, aes(x = factor(myOrg, levels = orgNameOrderVector),
                              y = myCountsAdjusted2)) +
    geom_bar(aes(fill=factor(myLigand, level = ligandOrderVector)), position="stack", stat="identity") +
    ylab("Number of ligand binding proteins /\n number of proteins in proteome") +
    xlab("Organism") +
    ggtitle("Identified fraction of ligand binding\n proteins, pLDDT > 70") +
    ylim(c(0, 0.55)) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=7)) +
    theme(legend.position = "none")
pOut3 <- grid.arrange(p1, p2, nrow = 1, ncol = 2)
ggsave(filename = "../IRON_DATA/stacked_bar_chart_normalized_nProteins.bmp", plot = pOut3, 
       width = 6.5, height = 3.25, units = "in", dpi = 300)





# Histograms of # ligs per protein 
# Iron Sulfur seperate from Fe/Zn
# Main text
find_max_number_of_ligs_per_uniprot_predicted_data <- function(inData,
                                                               feSVector,
                                                               feZnVector) {
    uniqueUniprots <- unique(inData[,18])
    myMaxFeS <- numeric(length = length(uniqueUniprots))
    myMaxFeZn <- numeric(length = length(uniqueUniprots))
    for (i in 1:length(uniqueUniprots)) {
        print(i)
        iDF <- inData[which(inData[,18] == uniqueUniprots[i]),]
        myMaxFeS[i] <- length(iDF[which(iDF[,21] %in% feSVector),20])
        myMaxFeZn[i] <- length(iDF[which(iDF[,21] %in% feZnVector),20])
    }
    
    typeVector <- c(rep("FeS", length(uniqueUniprots)), rep("Fe/Zn", length(uniqueUniprots)))
    valueVector <- c(myMaxFeS, myMaxFeZn)
    uniqueUniprots <- c(uniqueUniprots, uniqueUniprots)
    return(data.frame(uniqueUniprots, valueVector, typeVector))
}

feSVector <- c("4Fe4S_4CYS", 
               "4Fe4S_3CYS",     
               "3Fe4S_3CYS",
               "2Fe2S_4CYS",     
               "2Fe2S_2CYS_2HIS",
               "2Fe2S_3CYS_1ASP")

feZnVector <- c("1FeZn_4CYS",     
                "1ZINC_3CYS_1HIS",
                "1ZINC_2CYS_2HIS",
                "1ZINC_1CYS_3HIS",
                "1ZINC_4HIS",     
                "1ZINC_3HIS")

nLigsPerProteinDF <- find_max_number_of_ligs_per_uniprot_predicted_data(aggData_rmsd_0.5,
                                                                        feSVector = feSVector,
                                                                        feZnVector = feZnVector)
library(ggprism)
p1 <- ggplot(nLigsPerProteinDF) +
    geom_histogram(aes(x = valueVector,
                       fill = factor(typeVector, levels = c("FeS", "Fe/Zn"))), binwidth = 1, color = "black", position = "dodge") +
    xlab("Number of ligands per Protein") +
    ylab("Counts") +
    ggtitle("Distribution of the number of ligands per protein") +
    theme_bw() +
    #theme_prism() +
    #xlim(1,36) + 
    scale_y_continuous(trans="log1p", breaks = c(1,10,100,1000, 10000)) +
    scale_x_continuous(guide = guide_prism_minor(), limits = c(1,37)) +
    theme(text = element_text(size=12)) +
    #theme(plot.title = element_text(face = "bold")) +
    guides(fill=guide_legend(title="Ligand")) +
    theme(plot.title = element_text(hjust = 0.5))
#theme(legend.position = "none")
p1
ggsave(filename = "../IRON_DATA/ligands_per_uniprot_histogram.bmp", plot = p1, 
       width = 5.5, height = 5.5, units = "in", dpi = 300)

# p1 <- ggplot(nLigsPerProteinDF) +
#     geom_histogram(aes(x = myMaxFeS), binwidth = 1, color = "black", fill = "blue", alpha = 0.3) +
#     geom_histogram(aes(x = myMaxFeZn), binwidth = 1, color = "black", fill = "red", alpha = 0.3) +
#     xlab("Number of Ligands Per UniProt") +
#     ylab("Counts") +
#     ggtitle("a") +
#     theme_bw() +
#     scale_y_continuous(trans="log1p") +
#     theme(plot.title = element_text(face = "bold")) +
#     #theme(plot.title = element_text(hjust = 0.5))
#     theme(legend.position = "none")
# p1





# RMSD distribution by ligand type
# LOOK HERE
p1 <- ggplot(subset(aggData_rmsd_0.5, V12 > 70), aes(x = V15, 
                                                     fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free") +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("Distributions of RMSDs by Ligand, pLDDT > 70") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    geom_vline(xintercept = 0.5, linetype = 2, color = "black", size = 0.7) +
    theme_bw() +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")
# ZJW
ggsave(filename = "../IRON_DATA/rmsd_0.5_small_grid_bfac_70.bmp", plot = p1, 
       width = 6.5, height = 8.5/2, units = "in", dpi = 300)


p1 <- ggplot(aggData_rmsd_0.5, aes(x = V15, 
                                   fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free") +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("a") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold")) +
    #theme(plot.title = element_text(hjust = 0.5))
    theme(legend.position = "none")

p3 <- ggplot(subset(aggData_rmsd_0.5,V12 > 90), aes(x = V15, 
                                                    fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free") +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("b") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold")) +
    #theme(plot.title = element_text(hjust = 0.5))
    theme(legend.position = "none")

p5 <- ggplot(subset(subset(aggData_rmsd_0.5,V12 > 90), V23 == "KNOWN_UNIPROT_MB"), aes(x = V15, 
                                                                                       fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("c") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold")) +
    theme(legend.position = "none")

p8 <- ggplot(subset(subset(aggData_rmsd_0.5, V22 == "COMPLETELY_NEW"), V12 > 90),
             aes(x = V15, 
                 fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("d") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold")) +
    theme(legend.position = "none")

#pOut <- grid.arrange(p7, p8, nrow = 1, ncol = 2)



pOut <- grid.arrange(p1, p3, p5, p8, nrow = 2, ncol = 2)

ggsave(filename = "../IRON_DATA/rmsd_0.5_big_grid.bmp", plot = pOut, 
       width = 15, height = 15, units = "in", dpi = 300)



# Plots of high confidence

p2 <- ggplot(subset(subset(aggData_rmsd_0.5,V12 > 90), V23 == "KNOWN_UNIPROT_MB"), aes(x = V15, 
                                                                                       fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n found in known UniProt metal binding sites, pLDDT > 90") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=globalTextSize)) +
    theme(legend.position = "none")


p3 <- ggplot(subset(subset(aggData_rmsd_0.5,V12 > 90), V23 != "KNOWN_UNIPROT_MB"), aes(x = V15, 
                                                                                       fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n not found in known UniProt metal binding sites, pLDDT > 90") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=globalTextSize)) +
    theme(legend.position = "none")



pOut <- grid.arrange(p2, p3, nrow = 2, ncol = 1)
ggsave(filename = "../IRON_DATA/rmsd_0.5_big_grid_bfac_90_comparison.bmp", plot = pOut, 
       width = 15, height = 7.5, units = "in", dpi = 300)


p2 <- ggplot(subset(subset(aggData_rmsd_0.5,V12 > 70), V23 != "NONE"), aes(x = V15, 
                                                                           fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n found in UniProt, pLDDT > 70") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=7)) +
    theme(legend.position = "none")


p3 <- ggplot(subset(subset(aggData_rmsd_0.5,V12 > 70), V23 == "NONE"), aes(x = V15, 
                                                                           fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n not found in UniProt, pLDDT > 70") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=7)) +
    theme(legend.position = "none")



pOut <- grid.arrange(p2, p3, nrow = 2, ncol = 1)
#ggsave(filename = "../IRON_DATA/rmsd_0.5_big_grid_bfac_70_comparison.bmp", plot = pOut, 
#       width = 15, height = 7.5, units = "in", dpi = 300)
ggsave(filename = "../IRON_DATA/rmsd_0.5_big_grid_bfac_70_comparison.bmp", plot = pOut, 
       width = 6.5, height = 8.5, units = "in", dpi = 300)

p2 <- ggplot(subset(subset(aggData_rmsd_0.5,V12 > 0), V23 == "KNOWN_UNIPROT_MB"), aes(x = V15, 
                                                                                      fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n found in known UniProt metal binding sites, pLDDT > 0") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=globalTextSize)) +
    theme(legend.position = "none")


p3 <- ggplot(subset(subset(aggData_rmsd_0.5,V12 > 0), V23 != "KNOWN_UNIPROT_MB"), aes(x = V15, 
                                                                                      fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n not found in known UniProt metal binding sites, pLDDT > 0") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=globalTextSize)) +
    theme(legend.position = "none")



pOut <- grid.arrange(p2, p3, nrow = 1, ncol = 2)
ggsave(filename = "../IRON_DATA/rmsd_0.5_big_grid_bfac_0_comparison.bmp", plot = pOut, 
       width = 15, height = 7.5, units = "in", dpi = 300)




# Plots hhsearch/hhblits

p1 <- ggplot(subset(subset(aggData_rmsd_0.5, V22 == "NOT_NEW" | V23 != "NONE"), V12 > 90),
             aes(x = V15, 
                 fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n found in known UniProt metal binding sites or by hhsearch, pLDDT > 90") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=globalTextSize)) +
    theme(legend.position = "none")


p2 <- ggplot(subset(subset(subset(subset(aggData_rmsd_0.5, V22 == "COMPLETELY_NEW"), V12 > 90), V23 != "KNOWN_UNIPROT_MB"), V23 != "KNOWN_UNIPROT_REGION"),
             aes(x = V15, 
                 fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n not found in known UniProt metal binding sites or by hhsearch, pLDDT > 90") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=globalTextSize)) +
    theme(legend.position = "none")


pOut <- grid.arrange(p1, p2, nrow = 2, ncol = 1)

ggsave(filename = "../IRON_DATA/big_grid_not_new_vs_completely_new_and_not_uniprot_bfac_90.bmp", plot = pOut, 
       width = 15, height = 7.5, units = "in", dpi = 300)


p1 <- ggplot(subset(subset(aggData_rmsd_0.5, V22 == "NOT_NEW" | V23 != "NONE"), V12 > 70),
             aes(x = V15, 
                 fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n found in UniProt or by hhsearch, pLDDT > 70") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=globalTextSize)) +
    theme(legend.position = "none")


p2 <- ggplot(subset(subset(subset(subset(aggData_rmsd_0.5, V22 == "COMPLETELY_NEW"), V12 > 70), V23 == "NONE"), V23 != "KNOWN_UNIPROT_REGION"),
             aes(x = V15, 
                 fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n not found in UniProt or by hhsearch, pLDDT > 70") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=globalTextSize)) +
    theme(legend.position = "none")


pOut <- grid.arrange(p1, p2, nrow = 2, ncol = 1)

ggsave(filename = "../IRON_DATA/big_grid_not_new_vs_completely_new_and_not_uniprot_bfac_70.bmp", plot = pOut, 
       width = 6.5, height = 8.5, units = "in", dpi = 300)

p1 <- ggplot(subset(subset(aggData_rmsd_0.5, V22 == "NOT_NEW" | V23 != "NONE"), V12 > 0),
             aes(x = V15, 
                 fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n found in UniProt metal binding sites or by hhsearch, pLDDT > 0") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=globalTextSize)) +
    theme(legend.position = "none")


p2 <- ggplot(subset(subset(subset(subset(aggData_rmsd_0.5, V22 == "COMPLETELY_NEW"), V12 > 0), V23 != "KNOWN_UNIPROT_MB"), V23 != "KNOWN_UNIPROT_REGION"),
             aes(x = V15, 
                 fill = factor(V21, level = ligandOrderVector))) +
    geom_histogram(bins = 30, color = "black") +
    facet_wrap(~factor(V21, level = ligandOrderVector), nrow = 3, scales = "free",
               drop = FALSE) +
    xlim(c(0,0.75)) +
    xlab("RMSD") +
    ylab("Counts") +
    ggtitle("RMSD distributions of identified ligands (this work)\n not found in known UniProt metal binding sites or by hhsearch, pLDDT > 0") +
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(breaks= integer_breaks()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size=globalTextSize)) +
    theme(legend.position = "none")


pOut <- grid.arrange(p1, p2, nrow = 1, ncol = 2)

ggsave(filename = "../IRON_DATA/big_grid_not_new_vs_completely_new_and_not_uniprot_bfac_0.bmp", plot = pOut, 
       width = 15, height = 7.5, units = "in", dpi = 300)



# 2D histogram

ggplot(aggData_rmsd_0.5, 
       aes(x = V15,
           y = V12
       )) +
    geom_bin2d(aes(fill = ..ndensity..), bins = 40) +
    facet_wrap(~factor(V21, levels = ligandOrderVector), nrow = 3) +
    ylim(c(0,100)) +
    xlim(c(0,0.75)) +
    ggtitle("a") +
    xlab("RMSD") +
    ylab("Bfac") +
    theme(plot.title = element_text(face = "bold")) +
    scale_fill_viridis() +
    theme_bw()



# Ratio of gains plots

build_relative_novelty_df <- function(inData,
                                      inUniProtData,
                                      pLDDTCutoffVector) {
    myOrg <- character()
    mypLDDT <- numeric()
    myNewLigs <- integer()
    myOldLigs <- integer()
    myLigRatio <- numeric()
    myNewUniprots <- integer()
    myOldUniprots <- integer()
    myUniprotRatio <- numeric()
    
    uniqueOrgs <- unique(inData$V25)
    
    currentPosition <- 1
    
    for (i in 1:length(uniqueOrgs)) {
        iOrg <- uniqueOrgs[i]
        for (k in 1:length(pLDDTCutoffVector)) {
            # N ligs
            myNewLigs[currentPosition] <- nrow(inData[which(
                inData$V23 == "NONE" &
                    inData$V12 > pLDDTCutoffVector[k] &
                    inData$V25 == iOrg
            ), ])
            myOldLigs[currentPosition] <- nrow(inData[which(
                inUniProtData$ligandName != "NOTHING" &
                    inUniProtData$rep.iOrg..nrow.iDF.. == iOrg
            ), ])
            myLigRatio[currentPosition] <- (myNewLigs[currentPosition] + myOldLigs[currentPosition]) / myOldLigs[currentPosition]
            
            # Get the number new uniprots
            h1 <- unique(inData[which(inData$V23 == "NONE" &
                                          inData$V12 > pLDDTCutoffVector[k] &
                                          inData$V25 == iOrg),18])
            #length(h1)
            h2 <- unique(inUniProtData[which(inUniProtData$ligandName != "NOTHING" &
                                                 inUniProtData$rep.iOrg..nrow.iDF.. == iOrg),1])
            myNewUniprots[currentPosition] <- length(setdiff(h1,h2))
            myOldUniprots[currentPosition] <- length(h2)
            myUniprotRatio[currentPosition] <- (myNewUniprots[currentPosition] + myOldUniprots[currentPosition]) / (myOldUniprots[currentPosition])
            
            
            myOrg[currentPosition] <- iOrg
            mypLDDT[currentPosition] <- pLDDTCutoffVector[k]
            currentPosition <- currentPosition + 1
        }
        
    }
    
    return(data.frame(myOrg,
                      mypLDDT,
                      myNewLigs,
                      myOldLigs,
                      myLigRatio,
                      myNewUniprots,
                      myOldUniprots,
                      myUniprotRatio))
    
}

relDF <- build_relative_novelty_df(aggData_rmsd_0.5,
                                   uniprot_data,
                                   c(0,70, 90))

p1 <- ggplot(subset(relDF, mypLDDT == 70), aes(x = factor(myOrg, levels = orgNameOrderVector),
                                               y = myUniprotRatio)) +
    geom_bar(stat = "identity",
             color = "black",
             fill = "#009E73") +
    ylim(c(0,8.5)) +
    geom_hline(yintercept = 1, color = 'black', linetype = "dashed") +
    xlab("Organism") +
    ylab("Novel protein ratio") +
    ggtitle("Novel protein ratio via UniProt") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=10)) +
    theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(subset(relDF, mypLDDT == 70), aes(x = factor(myOrg, levels = orgNameOrderVector),
                                               y = myLigRatio)) +
    geom_bar(stat = "identity",
             color = "black",
             fill = "#CC79A7") +
    ylim(c(0,8.5)) +
    geom_hline(yintercept = 1, color = 'black', linetype = "dashed") +
    xlab("Organism") +
    ylab("Novel ligand ratio") +
    ggtitle("Novel ligand ratio via UniProt") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=globalTextSize)) +
    theme(plot.title = element_text(hjust = 0.5))

pOut1 <- grid.arrange(p1, p2, nrow = 1, ncol = 2)
ggsave(filename = "../IRON_DATA/novelty_ratio_barcharts_by_uniprot.bmp", 
       plot = pOut1, 
       width = 10, height = 5, units = "in", dpi = 300)




build_relative_novelty_df_hhsearch <- function(inData,
                                               inUniProtData,
                                               pLDDTCutoffVector) {
    myOrg <- character()
    mypLDDT <- numeric()
    myNewLigs <- integer()
    myOldLigs <- integer()
    myLigRatio <- numeric()
    myNewUniprots <- integer()
    myOldUniprots <- integer()
    myUniprotRatio <- numeric()
    
    uniqueOrgs <- unique(inData$V25)
    
    currentPosition <- 1
    
    for (i in 1:length(uniqueOrgs)) {
        iOrg <- uniqueOrgs[i]
        for (k in 1:length(pLDDTCutoffVector)) {
            # N ligs
            myNewLigs[currentPosition] <- nrow(inData[which(
                inData$V23 == "NONE" &
                    inData$V12 > pLDDTCutoffVector[k] &
                    inData$V22 == "COMPLETELY_NEW" &
                    inData$V25 == iOrg
            ), ])
            myOldLigs[currentPosition] <- nrow(inData[which(
                inUniProtData$ligandName != "NOTHING" &
                    inUniProtData$rep.iOrg..nrow.iDF.. == iOrg
            ), ])
            myLigRatio[currentPosition] <- (myNewLigs[currentPosition] + myOldLigs[currentPosition]) / myOldLigs[currentPosition]
            
            # Get the number new uniprots
            h1 <- unique(inData[which(inData$V23 == "NONE" &
                                          inData$V12 > pLDDTCutoffVector[k] &
                                          inData$V22 == "COMPLETELY_NEW" &
                                          inData$V25 == iOrg),18])
            #length(h1)
            h2 <- unique(inUniProtData[which(inUniProtData$ligandName != "NOTHING" &
                                                 inUniProtData$rep.iOrg..nrow.iDF.. == iOrg),1])
            myNewUniprots[currentPosition] <- length(setdiff(h1,h2))
            myOldUniprots[currentPosition] <- length(h2)
            myUniprotRatio[currentPosition] <- (myNewUniprots[currentPosition] + myOldUniprots[currentPosition]) / (myOldUniprots[currentPosition])
            
            
            myOrg[currentPosition] <- iOrg
            mypLDDT[currentPosition] <- pLDDTCutoffVector[k]
            currentPosition <- currentPosition + 1
        }
        
    }
    
    return(data.frame(myOrg,
                      mypLDDT,
                      myNewLigs,
                      myOldLigs,
                      myLigRatio,
                      myNewUniprots,
                      myOldUniprots,
                      myUniprotRatio))
    
}

relDF2 <- build_relative_novelty_df_hhsearch(aggData_rmsd_0.5,
                                             uniprot_data,
                                             c(0,70, 90))

p1 <- ggplot(subset(relDF2, mypLDDT == 70), aes(x = factor(myOrg, levels = orgNameOrderVector),
                                                y = myUniprotRatio)) +
    geom_bar(stat = "identity",
             color = "black",
             fill = "#009E73") +
    ylim(c(0,8.5)) +
    geom_hline(yintercept = 1, color = 'black', linetype = "dashed") +
    xlab("Organism") +
    ylab("Novel protein ratio") +
    ggtitle("Novel protein ratio via\n UniProt and hhsearch") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=10)) +
    theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(subset(relDF2, mypLDDT == 70), aes(x = factor(myOrg, levels = orgNameOrderVector),
                                                y = myLigRatio)) +
    geom_bar(stat = "identity",
             color = "black",
             fill = "#CC79A7") +
    ylim(c(0,8.5)) +
    geom_hline(yintercept = 1, color = 'black', linetype = "dashed") +
    xlab("Organism") +
    ylab("Novel ligand ratio") +
    ggtitle("Novel ligand ratio via\n UniProt and hhsearch") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=10)) +
    theme(plot.title = element_text(hjust = 0.5))

pOut2 <- grid.arrange(p1, p2, nrow = 1, ncol = 2)
ggsave(filename = "../IRON_DATA/novelty_ratio_barcharts_by_uniprot_and_hhsearch.bmp", 
       plot = pOut2, 
       width = 10, height = 5, units = "in", dpi = 300)

pOut3 <- grid.arrange(pOut1, pOut2, nrow = 2, ncol = 1)
ggsave(filename = "../IRON_DATA/novelty_ratio_barcharts_by_uniprot_and_hhsearch_total.bmp", 
       plot = pOut3, 
       width = 7.5, height = 7.5, units = "in", dpi = 300)









p1 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70 & V23 == "NONE"), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="stack") +
    ylab("Counts") +
    xlab("Organism") +
    ggtitle("Absolute counts of identified ligands (this work)\n not annotated in UniProt") +
    ylim(c(0, 8500)) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=8)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

p2 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70 & V23 == "NONE"), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts of identified ligands (this work)\n not annotated in UniProt") +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=8)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

relDF <- build_relative_novelty_df(aggData_rmsd_0.5,
                                   uniprot_data,
                                   c(0,70, 90))

p3 <- ggplot(subset(relDF, mypLDDT == 70), aes(x = factor(myOrg, levels = orgNameOrderVector),
                                               y = myUniprotRatio)) +
    geom_bar(stat = "identity",
             color = "black",
             fill = "#009E73") +
    ylim(c(0,8.5)) +
    geom_hline(yintercept = 1, color = 'black', linetype = "dashed") +
    xlab("Organism") +
    ylab("Novel protein ratio") +
    ggtitle("Novel protein ratio relative to UniProt") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=8)) +
    theme(plot.title = element_text(hjust = 0.5))

p4 <- ggplot(subset(relDF, mypLDDT == 70), aes(x = factor(myOrg, levels = orgNameOrderVector),
                                               y = myLigRatio)) +
    geom_bar(stat = "identity",
             color = "black",
             fill = "#CC79A7") +
    ylim(c(0,8.5)) +
    geom_hline(yintercept = 1, color = 'black', linetype = "dashed") +
    xlab("Organism") +
    ylab("Novel ligand binding site ratio") +
    ggtitle("Novel ligand binding site ratio relative to UniProt") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=8)) +
    theme(plot.title = element_text(hjust = 0.5))

pOut1 <- grid.arrange(p4, p3, p1, p2, nrow = 2, ncol = 2)

ggsave(filename = "../IRON_DATA/counts_and_novelty_ratio_via_uniprot.bmp", 
       plot = pOut1, 
       width = 6.5, height = 7.5, units = "in", dpi = 300)







p1 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70 & V23 == "NONE" & V22 == "COMPLETELY_NEW"), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="stack") +
    ylab("Counts") +
    xlab("Organism") +
    ggtitle("Absolute counts of identified ligands (this work)\n not annotated in UniProt or by hhsearch") +
    ylim(c(0, 8500)) +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=8)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

p2 <- ggplot(data=subset(aggData_rmsd_0.5, V12 > 70 & V23 == "NONE" & V22 == "COMPLETELY_NEW"), aes(x = factor(V25, levels = orgNameOrderVector))) +
    geom_bar(aes(fill=factor(V21, level = ligandOrderVector)), position="fill") +
    ylab("Fraction") +
    xlab("Organism") +
    ggtitle("Relative counts of identified ligands (this work)\n not annotated in UniProt or by hhsearch") +
    theme_bw() +
    scale_fill_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    #theme(plot.title = element_text(face = "bold")) +
    theme(text = element_text(size=8)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")

relDF2 <- build_relative_novelty_df_hhsearch(aggData_rmsd_0.5,
                                             uniprot_data,
                                             c(0,70, 90))

p3 <- ggplot(subset(relDF2, mypLDDT == 70), aes(x = factor(myOrg, levels = orgNameOrderVector),
                                                y = myUniprotRatio)) +
    geom_bar(stat = "identity",
             color = "black",
             fill = "#009E73") +
    ylim(c(0,8.5)) +
    geom_hline(yintercept = 1, color = 'black', linetype = "dashed") +
    xlab("Organism") +
    ylab("Novel protein ratio") +
    ggtitle("Novel protein ratio relative to\n UniProt and hhsearch") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=8)) +
    theme(plot.title = element_text(hjust = 0.5))

p4 <- ggplot(subset(relDF2, mypLDDT == 70), aes(x = factor(myOrg, levels = orgNameOrderVector),
                                                y = myLigRatio)) +
    geom_bar(stat = "identity",
             color = "black",
             fill = "#CC79A7") +
    ylim(c(0,8.5)) +
    geom_hline(yintercept = 1, color = 'black', linetype = "dashed") +
    xlab("Organism") +
    ylab("Novel ligand binding site ratio") +
    ggtitle("Novel ligand binding site ratio relative to\n UniProt and hhsearch") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.text.x=element_text(face="italic")) +
    theme(text = element_text(size=8)) +
    theme(plot.title = element_text(hjust = 0.5))
pOut1 <- grid.arrange(p4, p3, p1, p2, nrow = 2, ncol = 2)

ggsave(filename = "../IRON_DATA/counts_and_novelty_ratio_via_uniprot_and_hhsearch.bmp", 
       plot = pOut1, 
       width = 6.5, height = 7.5, units = "in", dpi = 300)

