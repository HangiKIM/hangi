library(dplyr)
library(stringr)
library(combinat)
library(tm)
library(qgraph)
library(rJava)
library(KoNLP)
library(RColorBrewer)
library(wordcloud)
library(arules)
library(igraph)
library(combinat)
library(ggplot2)
library(sna)
library(cluster)
library(ggthemes)

setwd("C:/Users/user1/desktop/about_r/result/party_faction/anti_sae")



################################################
########  previous work : data sorting #########
################################################

faction <- read.csv("antisae_faction_2012.csv")              ## import raw data


faction$weight <- gsub("1-10 / ","",faction$weight)          ## removing unnecessary character 
faction$weight <- gsub("1-9 / ","",faction$weight)
faction$weight <- gsub("1-8 / ","",faction$weight)
faction$weight <- gsub("1-7 / ","",faction$weight)
faction$weight <- gsub("1-6 / ","",faction$weight)
faction$weight <- gsub("1-5 / ","",faction$weight)
faction$weight <- gsub("1-4 / ","",faction$weight)
faction$weight <- gsub("1-3 / ","",faction$weight)
faction$weight <- gsub("1-2 / ","",faction$weight)
faction$weight <- gsub("1-1 / ","",faction$weight)
faction$weight <- gsub("°Ç","",faction$weight)
faction$weight <- gsub(",","",faction$weight)


weight_nu <- as.numeric(faction$weight)                     ## transforming weight : character -> numeric
faction <- cbind(faction, weight_nu)                        ## adding trnasformed numeric weight


## too important : below only be used to make MST Matrix or graph 
faction <- cbind(faction, weight_1 = (faction$weight_nu*-1)+max(faction$weight_nu)-mean(faction$weight_nu))


qq <- transform(faction, z_score = log(weight_1))                     # transforming weight : log variable
data <- data.frame(name1=qq$name, name2=qq$name1, z_score=qq$z_score) # arraying score 

   # data  <- data %>% mutate(z_score = ifelse(z_score < 5.991,  0, z_score)) ## removing below 75%


################################################
#############  data-set check  #################
################################################

head(rank(data$z_score))
data %>% head
head(z_score)
dim(data)
table(data$z_score)
data %>% select(z_score) %>% summary


################################################
####### making co_occurance data set  ##########
################################################

uniq = sort(unique(rbind(as.matrix(data[,1]), as.matrix(data[,2]))))

res_mat = matrix(0, length(uniq), length(uniq))
dimnames(res_mat)[[1]] = uniq
dimnames(res_mat)[[2]] = uniq


for(i in 1:dim(data)[[1]])
{
  res_mat[which(dimnames(res_mat)[[1]] == data[i,1]), which(dimnames(res_mat)[[2]] == data[i,2])] = data[i,3]
  res_mat[which(dimnames(res_mat)[[2]] == data[i,2]), which(dimnames(res_mat)[[1]] == data[i,1])] = data[i,3]
}

res_mat = t(res_mat)
co.matrix = data.frame(res_mat)

write.csv(co.matrix, "co_matrix.csv")                                 ## saving co_occurance data set 


################################################
##########  Estimated each Degree ##############
################################################

co.matrix <- as.matrix(co.matrix)

## co.matrix_minor <- as.matrix(co.matrix_minor)


indegree <- pstar(co.matrix,effects=c("indegree"))   ## Network centralization 

degree <- degree(co.matrix)                          ## Degree centrality
betweenness <- betweenness(co.matrix)                ## Betweenness centrality
closeness <- closeness(co.matrix)                    ## closeness centrality
cor_degree <- cor(co.matrix)                         ## cor-relation by each congressmen

index <- data.frame(uniq, degree, betweenness, cowork=betweenness/degree)      ## saving centrality
write.csv(index,file="index.csv")



################################################
#########  Plot by each Centrality   ###########
################################################

ggplot(index_minor, aes(degree, betweenness, label=uniq)) + geom_text(aes(color=commit), size=6, family="HYwulM") +
  geom_smooth(size=0.5, color=1,se=F) +
  theme_economist(base_size=13, base_family="UnPilgi") +
  theme(legend.position='none',
        axis.title.x = element_text(size=14, face="bold", lineheight=1),
        axis.title.y = element_text(size=14, face="bold"),
        plot.title = element_text(lineheight=50, size = 24, face="UnPilgi", colour = "#006633") )



################################################
#############  cluster analysis   ##############
################################################

fit <- kmeans(co.matrix, 3)
aggregate(co.matrix,by=list(fit$cluster),FUN=mean)
mydata <- data.frame(co.matrix, fit$cluster)

d <- dist(mydata, method = "euclidean")
fit <- hclust(d, method="ward.D") 
clusplot(pam(d, 4), labels=5, lines=5, color=T)

clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

dd <- c(d)


################################################
########    SNA Graph Drawing : MST   ##########
################################################ 

# mst <- minimum.spanning.tree(g, algorithm="unweighted")

pal3 <- brewer.pal(9, "Blues")                                    # whole color in graph 

set.seed(3952)                                                    # set seed to make the layout reproducible
g <- graph.adjacency(co.matrix, weighted=T,mode=c("directed"))    # build a graph from the above matrix
g <- mst(g)                                                       # remove all of minimum number
V(g)$label <- V(g)$name                                           # set labels and degrees of vertices
V(g)$degree <-  apply(co.matrix,1,sum)+apply(co.matrix,2,sum)     
V(g)$color <- "steel blue"                                        ## vertex color
V(g)$size <- V(g)$degree*0.0009                                   ## vertex size
V(g)$label.cex <- 0.85                                            ## vertex label size
V(g)$label.color <- "grey10"                                      ## vertex label color
E(g)$color <- pal3                                                ## edge color
egam <- log(E(g)$weight+.4)*0.5                                   ## edge width step 1
E(g)$width <- egam                                                ## edge width step 2

aa <- tkplot(g, layout=layout.fruchterman.reingold)               ## graph drawing: fruchterman.reingold
aa <- tkplot(g, layout=layout.kamada.kawai)                       ## graph drawing: kamada.kawai



################################################
#########  SNA Graph Drawing : normal ##########
################################################

pal3 <- brewer.pal(9, "Blues")                                    # whole color in graph 

set.seed(3952)                                                    # set seed to make the layout reproducible
g <- graph.adjacency(co.matrix, weighted=T,mode=c("directed"))    # build a graph from the above matrix
g <- simplify(g)                                                  # remove loops
V(g)$label <- V(g)$name                                           # set labels and degrees of vertices
V(g)$degree <-  apply(co.matrix,1,sum)+apply(co.matrix,2,sum)
V(g)$color <- "steel blue"                                        ## vertex color
V(g)$size <- V(g)$degree*0.0009                                   ## vertex size
V(g)$label.cex <- 0.85                                            ## vertex label size
V(g)$label.color <- "grey10"                                      ## vertex label color
E(g)$color <- pal3                                                ## edge color
egam <- log(E(g)$weight+.4)*0.5                                   ## edge width step 1
E(g)$width <- egam                                                ## edge width step 2

aa <- tkplot(g, layout=layout.fruchterman.reingold)               ## graph drawing: fruchterman.reingold
aa <- tkplot(g, layout=layout.kamada.kawai)                       ## graph drawing: kamada.kawai




###################################################
##### committee, priodical, and regional gap ######
################################################### 

## below materials are manually written by Excel


# committee
dd <- index %>% 
  group_by(commit) %>% 
  summarise(mean.degree = mean(degree), mean.betweenness = mean(betweenness), mean.cowork=mean(cowork))

write.csv(dd, "commit.csv")


# priodical
dd <- index %>% 
  group_by(number) %>% 
  summarise(mean.degree = mean(degree), mean.betweenness = mean(betweenness), mean.cowork=mean(cowork))

write.csv(dd, "priodical.csv")


# regional
dd <- index %>% 
  group_by(region) %>% 
  summarise(mean.degree = mean(degree), mean.betweenness = mean(betweenness), mean.cowork=mean(cowork))

write.csv(dd, "regional.csv")



