setwd("E:/MBA@Brandeis/Syllabus/193HS-256F Healthcare Data Analytics and Data Mining/Final")
getwd()

library(dplyr)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(fpc)
library(amap)
library(clusterSim)
library(ggplot2)
library(plotly)

drg_pccr <- read.csv("PCCR_DRG_DX.csv",header = T,fileEncoding = 'UTF-8-BOM')
summary(drg_pccr)
mdc <- read.csv("DRG_MDC.csv", header = T)
#Normalize
nmp <- drg_pccr$X3700.4000.PCCR_OR_and_Anesth_Costs
nmp_m <- mean(nmp)
nmp_sd <- sd(nmp)

drg_z <- scale(nmp, center = nmp_m, scale = nmp_sd)
drg_z

#Calculate distance
dist_eu <- get_dist(drg_z, method = "euclidean")
fviz_dist(dist_eu, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# Scree plot
wss <- (nrow(drg_z) -1) * sum(apply(drg_z, 2, var))
# var
for (i in 2:10) wss[i] <- sum(kmeans(drg_z, centers = i)$withinss)
plot(1:10,wss,
     type ='b',
     xlab = "Number of Clusters",
     ylab = "within Cluster Sum of Squares")

# Choosing K ....
set.seed(200)
k <- list()
for (i in 1:5){
  k[[i]] <- kmeans(drg_z, i,nstart =15)
}

between_totss <- list()
for (i in 1:5) {
  between_totss[[i]] <- k[[i]]$betweenss/k[[i]]$totss
} 

plot(1:5,
     between_totss,
     type = 'b',
     ylab = "Between ss/ Total SS",
     xlab = "Cluster #")

plot(drg_z, col = k[[2]]$cluster)
# K-means clustering
## k=2
#fviz_cluster(k[[2]], data = drg_z)

index.G1(drg_z,k[[2]]$cluster,d=NULL,centrotypes="centroids")

f_stat <- list()
for (i in 2:5) {
  f_stat[i] <-  round(calinhara(drg_z,k[[i]]$cluster),digits=2) ## f-stat
}

f_stat_cluster <- cbind(cluster = 2:5, f_value = f_stat[2:5])
f_stat_cluster <- as.data.frame(f_stat_cluster)
f_stat_cluster$cluster <- as.character(f_stat_cluster$cluster)


f_stat_cluster %>% ggplot(aes(x= cluster, y = f_value)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = f_value))

## k=3
col_pl <- c('green','blue','red')
plot(drg_z, col = k[[3]]$cluster)
legend("right",
       legend = levels(as.factor(k[[3]]$cluster)),
       fill = col_pl)

##
sil <- silhouette(k[[3]]$cluster, dist(drg_z))
fviz_silhouette(sil)

agg_m <- aggregate(drg_pccr[,57],list(k[[3]]$cluster), mean)

drg.cluster <- cbind(drg_pccr[,c(1:2,57)],cluster = k[[3]]$cluster)


plot(X3700.4000.PCCR_OR_and_Anesth_Costs ~ MSDRG,data = drg.cluster, 
     col = col_pl[ k[[3]]$cluster],
     las =2)
legend("topleft",
       legend = levels(as.factor(k[[3]]$cluster)),
       fill = col_pl)

drg.cluster$DRG_Desc <- factor(drg.cluster$DRG_Desc, 
                     levels = drg.cluster$DRG_Desc[order(drg.cluster$X3700.4000.PCCR_OR_and_Anesth_Costs)])
drg.cluster$cluster <- as.character(drg.cluster$cluster)
q <- drg.cluster %>% ggplot(aes(x= DRG_Desc, y =X3700.4000.PCCR_OR_and_Anesth_Costs)) +
  geom_point(mapping = aes(shape = cluster,col = cluster)) +
  scale_y_continuous(name = 'PCCR_OR_and_Anesth_Costs', breaks = seq(0,65000,10000)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_text("MSDRG_DESC")) 

ggplotly(q)

drg_mdc.cluster <- merge(drg.cluster,mdc, by.x = "MSDRG", by.y = "DRG" )

drg_mdc.cluster %>% group_by(cluster, MDC,MDC_CAT_NAME) %>% summarise(n = n())

drg_mdc.cluster %>% group_by(cluster, MDC,MDC_CAT_NAME) %>% summarise(n = n()) %>% 
  ggplot(aes(x = MDC_CAT_NAME, y = n, fill = cluster)) +
  geom_bar(stat = "identity", ) +
  theme(legend.position="right",
        axis.text.x = element_text(angle = 90, hjust = 1))

# SOM
library(kohonen)
set.seed(222)

drg_g <- somgrid(xdim = 1, ydim = 5, topo = "rectangular")
map <- som(drg_z,
           grid = drg_g,
           alpha = c(0.05,0.01),
           radius = 1)
map$codes
plot(map,type = "changes")

plot(map,
     type = 'codes',
     palette.name = rainbow,
     main = "2 by 3  Mapping of Phamar")
plot(map, type = "count")
plot(map, type = 'dist.neighbours')

plot(map, type = 'mapping')
plot(map, type = "quality")

#### option 2
set.seed(100)

#kmeans
clust <- kmeans(map$codes[[1]], 3)

round(calinhara(drg_z,clust$cluster[map$unit.classif]),digits=2) 
index.G1(drg_z,k[[3]]$cluster,d=NULL,centrotypes="centroids")
index.G1(drg_z,clust$cluster[map$unit.classif],d=NULL,centrotypes="centroids")


between_totss[[3]]

plot(map, type = "codes",
     palette.name = rainbow,
     bgcol = col_pl[clust$cluster], 
     main = "SOM Cluster Map"
)
add.cluster.boundaries(map, clust$cluster)
legend("right",
       legend = levels(as.factor(clust$cluster)),
       fill = col_pl)

plot(map, type = 'mapping')


drg_som <-data.frame(drg_pccr[,c(1:2,57)], cluster = clust$cluster[map$unit.classif])


drg_som$DRG_Desc <- factor(drg_som$DRG_Desc, 
                               levels = drg_som$DRG_Desc[order(drg_som$X3700.4000.PCCR_OR_and_Anesth_Costs)])
drg_som$cluster <- as.character(drg_som$cluster)
p <- drg_som %>% ggplot(aes(x= DRG_Desc, y =X3700.4000.PCCR_OR_and_Anesth_Costs)) +
  geom_point(aes(shape = cluster,color = cluster))+
  scale_y_continuous(name = 'PCCR_OR_and_Anesth_Costs', breaks = seq(0,65000,10000)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_text("MSDRG_DESC"),
        legend.position = "left") 

ggplotly(p)


drgsom_mdc.cluster <- merge(drg_som,mdc, by.x = "MSDRG", by.y = "DRG" )

drgsom_mdc.cluster %>% group_by(cluster, MDC,MDC_CAT_NAME) %>% summarise(n = n()) %>% 
  ggplot(aes(x = MDC_CAT_NAME, y = n, fill = cluster)) +
  geom_bar(stat = "identity", ) +
  theme(legend.position="right",
        axis.text.x = element_text(angle = 90, hjust = 1))




library(GGally)
ggparcoord(data = drg.cluster, columns = c(2:3), 
           groupColumn = "cluster",scale = "globalminmax") + 
  labs(x = "PCCR", y = "value", title = "Clustering")+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 90, hjust = 1))

drg.cluster$X3700.Operating.Room
drg.cluster$cluster <- as.character(drg.cluster$cluster)



plot(X5600.Drugs.Charged.to.Patients ~ X3700.Operating.Room, 
     data =drg.cluster,
     col = drg.cluster$cluster,
     main = "Drugs.Charged.to.Patient Vs Operating.Room"
)
text(X5600.Drugs.Charged.to.Patients ~ X3700.Operating.Room, 
     labels=DRG,
     data=drg.cluster,
     pos =2,
     cex=0.6, 
     font=1)


clusplot(z, drg.cluster$cluster, 
         color=TRUE, 
         shade=TRUE,
         labels= drg.cluster$DRG,
         lines=0)

drg.cluster %>% select(DRG,cluster,X5600.Drugs.Charged.to.Patients) %>% 
ggplot(aes(x = DRG, y = X5600.Drugs.Charged.to.Patients)) +
  geom_point(mapping = aes(colour = cluster))
