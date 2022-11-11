# Average Number of Genes Discovered - Graph B

library(readr)
library(data.table)
library(ggplot2)

sim_101622 <- read_rds("Data/simulationresults102222_5.rds")

ranked_res <- lapply(sim_101622, function (x) {
  PPI_num <- c(length(intersect(x$PPI$gene[1:50], x$Known)), 
               length(intersect(x$PPI$gene[1:100], x$Known)), 
               length(intersect(x$PPI$gene[1:150], x$Known)),
               length(intersect(x$PPI$gene[1:200], x$Known)),
               length(intersect(x$PPI$gene[1:250], x$Known)),
               length(intersect(x$PPI$gene[1:300], x$Known)),
               length(intersect(x$PPI$gene[1:350], x$Known)),
               length(intersect(x$PPI$gene[1:400], x$Known)),
               length(intersect(x$PPI$gene[1:450], x$Known)),
               length(intersect(x$PPI$gene[1:500], x$Known)))
  
  novel_num <- c(length(intersect(x$Novel$gene[1:50], x$Known)), 
                 length(intersect(x$Novel$gene[1:100], x$Known)), 
                 length(intersect(x$Novel$gene[1:150], x$Known)),
                 length(intersect(x$Novel$gene[1:200], x$Known)),
                 length(intersect(x$Novel$gene[1:250], x$Known)),
                 length(intersect(x$Novel$gene[1:300], x$Known)),
                 length(intersect(x$Novel$gene[1:350], x$Known)),
                 length(intersect(x$Novel$gene[1:400], x$Known)),
                 length(intersect(x$Novel$gene[1:450], x$Known)),
                 length(intersect(x$Novel$gene[1:500], x$Known)))
  
  KS_num <- c(length(intersect(x$KS$gene[1:50], x$Known)), 
              length(intersect(x$KS$gene[1:100], x$Known)), 
              length(intersect(x$KS$gene[1:150], x$Known)),
              length(intersect(x$KS$gene[1:200], x$Known)),
              length(intersect(x$KS$gene[1:250], x$Known)),
              length(intersect(x$KS$gene[1:300], x$Known)),
              length(intersect(x$KS$gene[1:350], x$Known)),
              length(intersect(x$KS$gene[1:400], x$Known)),
              length(intersect(x$KS$gene[1:450], x$Known)),
              length(intersect(x$KS$gene[1:500], x$Known)))
  
  list("PPI_num" = PPI_num,
       "PPI_genes" = intersect(x$PPI$gene, x$Known), 
       "PPI_rank" = which(x$PPI$gene%in%x$Known),
       "novel_num" = novel_num,
       "novel_genes" = intersect(x$Novel$gene, x$Known), 
       "novel_rank" = which(x$Novel$gene%in%x$Known),
       "KS_num" = KS_num,
       "KS_genes" = intersect(x$KS$gene, x$Known), 
       "KS_rank" = which(x$KS$gene%in%x$Known),
       "Known" = x$Known)
})

PPI_num <- c(0)
novel_num <- c(0)
KS_num <- c(0)

for (i in 1:10) {
  PPI_pct <- mean(unlist(lapply(ranked_res, function (x) {
    length(which(x$PPI_rank <= 50*i))
  })))
  PPI_num <- c(PPI_num, PPI_pct)
  
  novel_pct <- mean(unlist(lapply(ranked_res, function (x) {
    length(which(x$novel_rank <= 50*i))
  })))
  novel_num <- c(novel_num, novel_pct)
  
  KS_pct <- mean(unlist(lapply(ranked_res, function (x) {
    length(which(x$KS_rank <= 50*i))
  })))
  KS_num <- c(KS_num, KS_pct)
}


df_line  <- data.frame(test = rep(c("Novel", "PPI", "KS"), each = 11),
                       top = rep(c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), 3),
                       num = c(novel_num, PPI_num, KS_num))


png(paste0("Output/sim_102222_accuracy.png"), width=1500, height=1200, res=250)
ggplot(df_line, aes(x=top, y=num, group=test)) +
  geom_line(aes(linetype=test))+
  geom_point(aes(shape=test)) +
  ggtitle("b. Average percent of true signal genes out of top ranked m genes") +
  scale_x_continuous(name="Top m ranked genes", limits=c(0, 500)) +
  scale_y_continuous(name="Average number of true signal genes", limits=c(0, 50)) +
  theme_bw()
dev.off()



