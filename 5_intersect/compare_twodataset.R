library(readr)
deg17025_cluster_1<- read_csv("~deg17025_cluster_1.csv")
deg17025_cluster_2<- read_csv("~deg17025_cluster_2.csv")
deg17025_cluster<- union(deg17025_cluster_1[,7],deg17025_cluster_2[,7]) #union all shared name
deg17025_cluster= as.data.frame(deg17025_cluster)

deg17025_tur<- read_csv("~deg17025_module_ turquoise .csv")

deg17025_inter<- intersect(deg17025_cluster,deg17025_tur)
deg17025_inter= as.data.frame(deg17025_inter)

deg39099_cluster_1<- read_csv("~deg39099_cluster_1.csv")
deg39099_cluster_2<- read_csv("~deg39099_cluster_2.csv")
deg39099_cluster<- union(deg39099_cluster_1[,20],deg39099_cluster_2[,20]) #union all shared name
deg39099_cluster= as.data.frame(deg39099_cluster)

all_inter<- intersect(deg17025_inter,deg39099_cluster)