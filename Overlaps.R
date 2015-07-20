setwd("//ad.ucl.ac.uk/slms/home3/rmhvlrr/Desktop/Microarray results/New Files")

f1 <- "differential ECxENaB2000.csv"
f2 <- "differential PGxPGNaB2000.csv"
f3 <- "differential UNxNaB2000.csv"


f1r <- read.csv(f1)
f2r <- read.csv(f2)
f3r <- read.csv(f3)

names(f1r) = paste("ECxENaB2000",names(f1r),sep="_")
names(f2r) = paste("PGxPGNaB2000",names(f2r),sep="_")
names(f3r) = paste("UNxNaB2000",names(f3r),sep="_")

f12r = merge(f1r,f2r,by.x="ECxENaB2000_ID",by.y="PGxPGNaB2000_ID")
f123r = merge(f12r, f3r, by.x="ECxENaB2000_ID", by.y="UNxNaB2000_ID")
f13r = merge(f1r,f3r,by.x="ECxENaB2000_ID",by.y="UNxNaB2000_ID")
f23r = merge(f2r,f3r,by.x="PGxPGNaB2000_ID",by.y="UNxNaB2000_ID")

f123r
f12r
f1r
which(f123r$UNxNaB2000_adj.P.Val<0.05 & f123r$PGxPGNaB2000_adj.P.Val<0.05)

write.csv(f123r,file="All 3 together.csv",row.names=F)
head(f123r)

g1 = as.character(f1r$ECxENaB2000_ID)
g2 = as.character(f2r$PGxPGNaB2000_ID)
g3 = as.character(f3r$UNxNaB2000_ID)

all = unique(c(g1,g2,g3))
all
g1only = all[which(all %in% g1 & !all %in% g2 & !all %in% g3)]

g12only = all[which(all %in% g3 & !all %in% g1 & !all %in% g2)]
#g1only = all[which(all %in% g1 & !all %in% g2 & !all %in% g3)]

#g2only = all[which(all %in% g2 & !all %in% g1 & !all %in% g3)]
#g3only = all[which(all %in% g3 & !all %in% g1 & !all %in% g2)]
#g12only
#f1r

#h1 = XXXX [which(XXXX $ XXXXUNxNaB2000_ID %in% g12only),]

h1 = f3r[which(f3r$UNxNaB2000_ID %in% g12only),]

write.csv(h1,file="UN only.csv",row.names=F)




