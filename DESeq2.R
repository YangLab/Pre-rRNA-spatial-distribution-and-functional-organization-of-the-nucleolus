library(tidyverse)
library(DESeq2)

DFfloat2int <- function(df){
    rowNames <- rownames(df)
    df <- apply(df,2,as.integer) %>% as.data.frame
    rownames(df) <- rowNames
    return(df)
}

Count2DEG <- function(df,col="EnsemblID"){
    df <- df %>% column_to_rownames(col) %>% DFfloat2int
    return(df)
}

get_coldata <- function(sample,sample_field){
    
    coldata <- data.frame(V1=colnames(sample))
    coldata <- coldata %>% mutate(V2=V1)%>% 
                           column_to_rownames('V2')%>% 
                            separate(V1,into=sample_field)
                            
    return(coldata)
    
}

DEG <- function(Counts,pos,colData,formula){

    VStem <- Counts[,pos]
#     n1 <- length(pos1)
#     n2 <- length(pos2)
#     set1 <- rep(control,n1)
#     set2 <- rep(treatment,n2)
    
    

#     # case1: IntS6 versus DL1
#     conditiontem <- factor(c(set1,set2),levels = c(control,treatment))

#     coldatatem <- data.frame(row.names=colnames(VStem),
#                            condition=conditiontem)
    
    # batch <- substitute(batch)
    # print(typeof(quote(batch)))
    

    ddstem <- DESeqDataSetFromMatrix(VStem, colData, design = as.formula(formula))
    ddstem <- DESeq(ddstem)


    # Get differential expression results
    restem <- DESeq2::results(ddstem)
    # table(res$padj<0.001)
    ## Order by adjusted p-value
    restem <- restem[order(restem$padj), ]
    ## Merge with normalized count data
    resdatatem <- merge(as.data.frame(restem), as.data.frame(counts(ddstem, normalized=TRUE)), by="row.names", sort=FALSE) %>% arrange(padj)
    # names(resdatatem)[1] <- "EnsemblID"
    return(resdatatem)
    # return(restem)
}

# fill null of output of DEG function
DESeq2fillNA <- function(res){
    data <- res
    
    if (sum(is.na(data)) > 0){
        data$padj[is.na(data$padj)] <- 1
        data$pvalue[is.na(data$pvalue)] <- 1
        data$log2FoldChange[is.na(data$log2FoldChange)] <- 0
        data$lfcSE[is.na(data$lfcSE)] <- 0
        data$stat[is.na(data$stat)] <- 0        
    }

    return(data)
}

# get DEG gene dataframe
# differential expression gene
UpDownGene <- function(dds,FC = 1.5, P = 0.05, key = "padj"){
    if (!is.na(key)) {
        
        up <-  dds %>% dplyr::filter(log2FoldChange > log2(FC) & 
                                     .data[[key]] < P)%>%
                        arrange(.data[[key]])
        down <- dds %>% dplyr::filter(log2FoldChange < -log2(FC) & 
                                     .data[[key]] < P)%>%
                        arrange(.data[[key]])
    }
    else {
        
        up <-  dds %>% dplyr::filter(log2FoldChange > log2(FC)) %>%
                        arrange(log2FoldChange)
        down <- dds %>% dplyr::filter(log2FoldChange < -log2(FC))%>%
                        arrange(log2FoldChange)
    }

    return(list('up'=up,'down'=down))
}

addUpDowntag <- function (dds, FC = 1.5, Padj = 0.001) {
  print(FC)
  print(Padj)
    data <- dds %>% mutate(tag = case_when(log2FoldChange >= 
        log2({
            {
                FC
            }
        }) & padj <= {
        {
            Padj
        }
    } ~ "Up", log2FoldChange <= -log2({
        {
            FC
        }
    }) & padj <= {
        {
            Padj
        }
    } ~ "Down", TRUE ~ "Nochange"))
    data$tag <- factor(data$tag, levels = c("Up", "Down", "Nochange"))
    data <- data %>% arrange(padj)
    data %>% group_by(tag) %>% summarise(N = n()) %>% print
    cat("\n")

    return(data)
}

# Volcano plot
DEseq2volcanoPlot <- function(data,FC = 1.5, P = 0.05,key = "padj"){
    
    up <- UpDownGene(dds = data,FC = FC,P = P,key = key)$up
    down <- UpDownGene(dds = data,FC = FC,P = P,key = key)$down
    # library(ggplot2)
    cat("Total:",nrow(data),"\n")
    cat("Up-ragulated:",nrow(up),"\n")
    cat("Down-ragulated:",nrow(down),"\n")
    cat("Other:",nrow(data) - nrow(up) - nrow(down),"\n\n\n")
    
    p <- ggplot()+
    geom_point(data = data,
               aes(x = log2FoldChange,y= -log10(.data[[key]] + (0.1)**100)),alpha=0.3,shape=20,size=3,
              color = "#7E807F")+
    
    geom_point(data = up,
               aes(x = log2FoldChange,y= -log10(.data[[key]] + (0.1)**100),alpha=0.5),
               shape = 20,size=3,
               color = "#DB2143")+
    
    # geom_text_repel(data = data[data$log2FoldChange > log2(FC) & data$padj < P,],
    #                 aes(x=log2FoldChange,y=-log10(padj + (0.1)**100),label=Symbol),
    #                 size=5,color="black",point.size = 5)+
    
    geom_point(data = down,
               aes(x = log2FoldChange,y= -log10(.data[[key]] + (0.1)**100),alpha=0.5),
               shape = 20,size=3,
               color = "#13D7ED")+
    
    # geom_text_repel(data = data[data$log2FoldChange < -log2(FC) & data$padj < P,],
    #            aes(x = log2FoldChange,y= -log10(padj + (0.1)**100),label=Symbol),
    #                size=5,color="black",point.size = 5)+
                
    geom_hline(yintercept = -log10(P),linetype = "dashed",color = "#496982")+
    geom_vline(xintercept = -log(FC),linetype = "dashed",color = "#496982")+
    geom_vline(xintercept = log(FC),linetype = "dashed",color = "#496982")+
    theme_prism() + scale_y_continuous(expand=c(0.01,0))
    
    return(p)
}

DESeq2MAplot <- function(data,FC = 1.5,
                         control = c(8:9)){
    
    data <- data %>% mutate(control_mean = 
                    log10(rowSums(data[,control])/length(control)+1))
    up <- UpDownGene(dds = data,FC = FC,key = NA)$up
    down <- UpDownGene(dds = data,FC = FC,key = NA)$down
    
    cat("Total:",nrow(data),"\n")
    cat("Up-ragulated:",nrow(up),"\n")
    cat("Down-ragulated:",nrow(down),"\n")
    cat("Other:",nrow(data) - nrow(up) - nrow(down),"\n\n\n")
    
    p <- ggplot()+
        geom_point(data = data,
                  aes(x = control_mean, y = log2FoldChange),
                  alpha=0.3,shape=21,size=5,
                  fill="#7E807F",color = "#7E807F")+
        geom_point(data = up,
               aes(x = control_mean,y= log2FoldChange,alpha=0.5),
               shape = 21,size=5,
               fill="#DB2143",color = "#7E807F")+
        geom_point(data = down,
               aes(x = control_mean,y= log2FoldChange,alpha=0.5),
               shape = 21,size=5,
               fill="#13D7ED",color = "#7E807F")+
        geom_hline(yintercept = log2(FC),
                   linetype = "dashed",color = "#496982")+
        geom_hline(yintercept = -log2(FC),
                   linetype = "dashed",color = "#496982")+
        theme_prism()
        
    return(p)
}
