library(survival)
library(survminer)
source("functions.r")

#markers-----
p <- FormsurviResult(TCGA.list,
                     group.by = "LIF",
                     auto = T,
                     survival.message="os")
ggsave("/data/khaoxian/project/Stagecrc/output/survival/LIF_os.pdf", p, width = 21, height = 10, units = "cm")

#marker lists-----
igg.list <- list()
igg.list[["IgG_signature"]] <- c( "IGHG1", "IGHG2", "IGHG3", "IGHG4")

iga.list <- list()
iga.list[["IgA_signature"]] <- c( "IGHA1", "IGHA2", "JCHAIN")

cohort.DEList <- Addsignature2DEList(cohort.DEList,
                                     signature.list = igg.list)

cohort.DEList <- Addsignature2DEList(cohort.DEList,
                                     signature.list = iga.list)
p <- FormsurviResult(cohort.DEListx,
                     group.by = "IgG_signature",
                     auto = T,
                     survival.message="os")

p <- FormsurviResult(cohort.DEListx,
                     group.by = "IgA_signature",
                     auto = T,
                     survival.message="os")
#cell types-----
#add coad cibersortx results
coad.ciber <- read.table("/data/khaoxian/project/Stagecrc/data/cibersortx/webout/coad.cibersortx.result.txt", sep = "\t",
                         header = T,row.names = 1)
coad.ciber <- as.data.frame(t(coad.ciber))
cluster <- rownames(coad.ciber)
all(colnames(cohort.DEList$COAD$counts)==colnames(coad.ciber))
cohort.DEList$COAD$counts <- rbind(cohort.DEList$COAD$counts, coad.ciber)

#add read cibersortx results
read.ciber <- read.table("/data/khaoxian/project/Stagecrc/data/cibersortx/webout/read.cibersortx.result.txt", sep = "\t",
                         header = T,row.names = 1)
read.ciber <- as.data.frame(t(read.ciber))
all(colnames(cohort.DEList$READ$counts)==colnames(read.ciber))
cohort.DEList$READ$counts <- rbind(cohort.DEList$READ$counts, read.ciber)

#add gse17536 cibersortx results
GSE17536.ciber <- read.table("/data/khaoxian/project/Stagecrc/data/cibersortx/webout/GSE17536.cibersortx.result.txt", sep = "\t",
                             header = T,row.names = 1)
GSE17536.ciber <- as.data.frame(t(GSE17536.ciber))
all(colnames(cohort.DEList$GSE17536$counts)==colnames(GSE17536.ciber))
cohort.DEList$GSE17536$counts <- rbind(cohort.DEList$GSE17536$counts, GSE17536.ciber)

#add gse39582 cibersortx results
GSE39582.ciber <- read.table("/data/khaoxian/project/Stagecrc/data/cibersortx/webout/GSE39582.cibersortx.result.txt", sep = "\t",
                             header = T,row.names = 1)
GSE39582.ciber <- as.data.frame(t(GSE39582.ciber))
all(colnames(cohort.DEList$GSE39582$counts)==colnames(GSE39582.ciber))
cohort.DEList$GSE39582$counts <- rbind(cohort.DEList$GSE39582$counts, GSE39582.ciber)

#add gse17537 cibersortx results
GSE17537.ciber <- read.table("/data/khaoxian/project/Stagecrc/data/cibersortx/webout/GSE17537.cibersortx.result.txt", sep = "\t",
                             header = T,row.names = 1)
GSE17537.ciber <- as.data.frame(t(GSE17537.ciber))
all(colnames(cohort.DEList$GSE17537$counts)==colnames(GSE17537.ciber))
cohort.DEList$GSE17537$counts <- rbind(cohort.DEList$GSE17537$counts, GSE17537.ciber)

#plot
plot.list <- purrr::map(cluster, function(x){
  tryCatch({
    p <- FormsurviResult(cohort.DEList,
                         group.by = x,
                         survival.message="os",
                         auto=T)}, 
    warning = function(w){p}, 
    error = function(e){NA},
    finally = {p})
})

names(plot.list) <- cluster
saveRDS(plot.list,"/data/khaoxian/project/Stagecrc/output/survival/cibersortx.plot.list.rds")

