GBC <- read.table(file = "Z37_GBC_all_raw_real.txt", header = T)

GBC2 <- GBC[,4:ncol(GBC)]
rownames(GBC2) = GBC[,2]


# Calculate guide proportions
table_props = round(t(apply(GBC2,1,function(x) x/sum(x))),4)

# Table of guide proportions, ordered by most common 
props_order = t(apply(table_props,1,function(x) x[order(as.matrix(x),decreasing=T)]))

# Table of guide identities ordered by most common
guide_id_order = t(apply(GBC2,1,function(x) c(1:34)[order(as.matrix(x),decreasing=T)]))



# Find range separating single and double guides by % 2nd most common guide
# Hist of 2nd most common guide (for cells passing low 3rd guide and total read filters)
# Maximum 3rd most common guide?

library(dplyr)

props_order2 = apply(props_order,2,as.numeric)
rownames(props_order2) = rownames(props_order)
props_order2 = as.data.frame(props_order2)
filt_2rd = dplyr::filter(props_order2,props_order2[,2] < 0.01) 
filt_3rd = dplyr::filter(props_order2,props_order2[,3] < 0.01) 
filt_4rd = dplyr::filter(props_order2,props_order2[,4] < 0.01)
filt_2rd_list = props_order2[,2]<0.01
filt_3rd_list = props_order2[,3]<0.01
filt_4rd_list = props_order2[,4]<0.01

sum = apply(GBC2,1,function(x) sum(x))
sum = as.data.frame(sum)
props_order3=cbind(sum$sum,props_order2)
props_order3 = dplyr::filter(props_order3,props_order3[,1] > 100)
props_order3 = props_order3[,-1]

# Track cells passing depth + 3rd guide, plus 2nd guide range to determine single vs. double hit
filt_cells = vector(mode="logical",length=nrow(GBC2))
filt_single = vector(mode="logical",length=nrow(GBC2))
filt_double = vector(mode="logical",length=nrow(GBC2))
filt_trible = vector(mode="logical",length=nrow(GBC2))


filt_cells[filt_4rd_list] = TRUE # cells passing depth and 3rd guide
filt_single[props_order2[,2]<0.01 & filt_4rd_list] = TRUE # single hits
filt_double[props_order2[,2]>0.01 & props_order2[,3]<0.01 & filt_4rd_list] = TRUE  # double hits
filt_trible[props_order2[,3]>0.01 & filt_4rd_list]= TRUE


# Iterate through cells and call guide hits

cell_hit = matrix(0,nrow=dim(table_props)[1],ncol=3)
rownames(cell_hit) = rownames(GBC2)

for (cell in c(1:dim(table_props)[1])){
  if (cell %in% which(filt_single)){
    cell_hit[cell,1] = guide_id_order[cell,1]
  }
  if (cell %in% which(filt_double)){
    cell_hit[cell,1] = guide_id_order[cell,1]
    cell_hit[cell,2] = guide_id_order[cell,2]
  }
  if (cell %in% which(filt_trible)){
    cell_hit[cell,1] = guide_id_order[cell,1]
    cell_hit[cell,2] = guide_id_order[cell,2]
    cell_hit[cell,3] = guide_id_order[cell,3]
  }
}


# How many cells have calls?
length(which(rowSums(cell_hit)>0))  # single or double
length(which(cell_hit[,1]>0 & cell_hit[,2]==0)) # single 
length(which(cell_hit[,1]>0 & cell_hit[,2]>0 & cell_hit[,3]==0)) # double
length(which(cell_hit[,1]>0 & cell_hit[,2]>0 & cell_hit[,3]>0)) # trible


# Print QC stats and calls to text
guide_names = c("NA","7SK","LINC00189","LINC00511","LINC00578","LINC00657","LINC01011","LINC01012","LINC01133","LINC01578","LINC02693","BRD4","CASC15","FAM66c","FENDR","GAS5","NEAT1","HOTAIRM1","JPX","LRRC75a","MIAT","MALAT1","MEG3","MIR222HG","MIR155HG","PURPL","RMRP","SNHG3","SNHG1","XIST","SNHG16","LINC01615","LINC01638","nontarget-1","nontarget-2")
cell_hit_names = matrix(0,nrow=nrow(cell_hit),ncol=3)
for (i in c(1:nrow(cell_hit))){
  cell_hit_names[i,1] = guide_names[cell_hit[i,1]+1]
  cell_hit_names[i,2] = guide_names[cell_hit[i,2]+1]
  cell_hit_names[i,3] = guide_names[cell_hit[i,3]+1]
}

# Print cell hit matix, guide counts, guide proportions to text
summary_table = cbind(GBC2,table_props,cell_hit_names)
colnames(summary_table)[c(67,68,69)] = c("Call_1","Call_2","Call_3")
options(scipen=10)
write.table(summary_table,"Z37_Guides_calls_hitnames.txt",row.names=T,col.names=T,quote=F,sep="\t")


##############################################
##############################################

# Find counts of each guide hit and display heat map

hit_mat = matrix(,nrow=34,ncol=34)
for (i in c(1:34)){
  for (j in c(1:34)){
    check1 = cell_hit[,1]== i & cell_hit[,2] == j
    #check2 = cell_hit[,1]== j & cell_hit[,2] == i
    #hit_mat[i,j] = length(which(check1 | check2))
    hit_mat[i,j] = length(which(check1))
  }
}

# check3 = cell_hit[,1]== 1 & cell_hit[,2] == 1
# Correct diagonal
for (i in c(1:34)){
  hit_mat[i,i] = length(which(cell_hit[,1]==i & cell_hit[,2]==0))
}


# Heat map of counts for hits
colors = c(seq(0,25,length=22)) #max = 268
my_palette <- colorRampPalette(brewer.pal(6,"YlOrRd"))(n = 21) 


library(gplots)

labRow = c("7SK","LINC00189","LINC00511","LINC00578","LINC00657","LINC01011","LINC01012","LINC01133","LINC01578","LINC02693","BRD4","CASC15","FAM66c","FENDR","GAS5","NEAT1","HOTAIRM1","JPX","LRRC75a","MIAT","MALAT1","MEG3","MIR222HG","MIR155HG","PURPL","RMRP","SNHG3","SNHG1","XIST","SNHG16","LINC01615","LINC01638","nontarget-1","nontarget-2")
labCol = c("7SK","LINC00189","LINC00511","LINC00578","LINC00657","LINC01011","LINC01012","LINC01133","LINC01578","LINC02693","BRD4","CASC15","FAM66c","FENDR","GAS5","NEAT1","HOTAIRM1","JPX","LRRC75a","MIAT","MALAT1","MEG3","MIR222HG","MIR155HG","PURPL","RMRP","SNHG3","SNHG1","XIST","SNHG16","LINC01615","LINC01638","nontarget-1","nontarget-2")

hm2_call = heatmap.2(as.matrix(hit_mat),col=my_palette,breaks=colors,
                     key = T,keysize = 1.5,
                     key.title = "Cell number",key.xlab = "",
                     #xlab="second guide", ylab ="first guide",
                     density.info="none",trace="none",Rowv=F,Colv=F,dendrogram="both",
                     symm=F,labRow= labRow,labCol= labCol,margins=c(4,4),
                     offsetRow = 0.1, offsetCol = 0.1,
                     scale="none",cexRow=0.7,cexCol=0.7,colsep=c(0:33),rowsep=c(0:33),
                     sepcolor="white",sepwidth=c(0.005,0.005),
                     cellnote=hit_mat,notecol="black",
                     border = T,
                     key.xtickfun =  function(){breaks <- parent.frame()$breaks
                     return(list(
                       at=parent.frame()$scale01(c(breaks[1],breaks[12],
                                                   breaks[length(breaks)])),
                       labels=c(as.character(breaks[1]),"12",
                                ">25"))
                     )})

library(ComplexHeatmap)
library(circlize)
display.brewer.pal(11,"Spectral") 
brewer.pal(11,"Spectral")
display.brewer.pal(11,"PuOr") 
brewer.pal(11,"PuOr")

col_fun <- colorRamp2(
  c(0,10,20,30), 
  c("light yellow", "#FDAE61", "#F46D43","#D53E4F")
)

x = as.matrix(hit_mat)

Heatmap(x,name = "Cell Number",col=col_fun,
        cluster_columns=F,cluster_rows = F,row_labels =labRow,row_names_gp = gpar(fontsize = 8,fontface = "bold"),
        column_labels =labCol,column_names_gp = gpar(fontsize = 8,fontface = "bold"),border = T, border_gp=gpar(col = "black"),
        rect_gp = gpar(col = "white", lwd = 0.6),
        row_title = "First Guide",row_title_side = c("left"),row_title_gp = gpar(fontsize = 12, fontface ="bold"),
        column_title = "Second Guide",column_title_side = c("top"),column_title_gp = gpar(fontsize = 12, fontface ="bold"),
        #cell_fun = function(j, i, x, y, w, h, fill) {grid.text(hit_mat[i, j], x, y,gp = gpar(col="black",fontsize = 6)) }
          )

##############################################
##############################################


# Heat map of props for all guides based on hit selection 
colors = c(seq(0,1,length=22))
my_palette <- colorRampPalette(brewer.pal(6,"Blues"))(n = 21) 

g1=3
g2=7

set1 = cell_hit[,1] == g1 & cell_hit[,2]==0
set2 = cell_hit[,1] == g2 & cell_hit[,2]==0
set3 = (cell_hit[,1] == g1 & cell_hit[,2]==g2) | (cell_hit[,1] == g2 & cell_hit[,2]==g1)
set4 = cell_hit[,1] >0 & cell_hit[,2]==0

input_mat = table_props[c(which(set1),which(set2),which(set3)),]
#input_mat = table_props[which(set3),]
#input_mat = table_nlog[which(set3),]
input_mat = table_props[set4,][order(table_props[set4,1]),]

# Ordered heat map to display proportions for all single hits
input_mat = rbind(table_props[cell_hit[,1]==1 & set4,],table_props[cell_hit[,1]==2 & set4,],table_props[cell_hit[,1]==3 & set4,],table_props[cell_hit[,1]==4 & set4,],table_props[cell_hit[,1]==5 & set4,],table_props[cell_hit[,1]==6 & set4,],table_props[cell_hit[,1]==7 & set4,])



#select_singles = which(cell_hit[,1]>0 & cell_hit[,2]==0)
#select_doubles = which(cell_hit[,1]>0 & cell_hit[,2]>0)

# Ordered heat map to display proportions for all double hits
input_mat = matrix(,ncol=7,nrow=0)
for (i in c(1:7)){
  for (j in c(1:7)){
    if (i>j){
      temp_hits = which((cell_hit[,1]==i & cell_hit[,2] ==j) | (cell_hit[,1]==j & cell_hit[,2] ==i))
      if (length(temp_hits)>0){
        input_mat = rbind(input_mat,table_props[temp_hits,])
      }
    }
  }
}


hm2_call = heatmap.2(as.matrix(input_mat),col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=F,Colv=F,dendrogram="both",symm=F,labRow=c(""),labCol=colnames(table),margins=c(4,4),scale="none",cexRow=0.6,cexCol=0.6,colsep=c(0:10),sepcolor="gray",sepwidth=c(0.001,0.01),cellnote=round(input_mat,2),notecol="black",notecex=0.5)


hm2_call = heatmap.2(as.matrix(input_mat),col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=F,Colv=F,dendrogram="both",symm=F,labRow=c(""),labCol=colnames(table),margins=c(4,4),scale="none",cexRow=0.6,cexCol=0.6,colsep=c(0:10),sepcolor="gray",sepwidth=c(0.001,0.01))


##############################################
##############################################

# Scatter proportions for hit categories

par(mfrow=(c(7,7)),cex.axis=0.8,mar=c(4,4,2,2))
pt_cex=1
for (g1 in c(1:7)){
  for (g2 in c(1:7)){
    set1 = cell_hit[,1] == g1 & cell_hit[,2]==0
    set2 = cell_hit[,1] == g2 & cell_hit[,2]==0
    set3 = (cell_hit[,1] == g1 & cell_hit[,2]==g2) | (cell_hit[,1] == g2 & cell_hit[,2]==g1)
    if (g1 ==g2){
      plot(table_props[,c(g1,g2)],pch=19,cex=0.8,col=rgb(0,0,0,00))
      points(table_props[set1,g1], table_props[set1,g2],pch=19,cex= pt_cex,col=rgb(1,0,0,0.8))
      points(table_props[set2,g1], table_props[set2,g2],pch=19,cex= pt_cex,col=rgb(0,0,1,0.8))
      abline(h=0.95,lty=2,col="gray")
      abline(v=0.95,lty=2,col="gray")
    }
    else{
      plot(table_props[,c(g1,g2)],pch=19,cex=0.8,col=rgb(0,0,0,00),xlab=g1,ylab=g2)
      points(table_props[set3,g1], table_props[set3,g2],pch=19,cex= pt_cex,col=rgb(1,0,1,0.9))
      points(table_props[set1,g1], table_props[set1,g2],pch=19,cex= pt_cex,col=rgb(1,0,0,0.8))
      points(table_props[set2,g1], table_props[set2,g2],pch=19,cex= pt_cex,col=rgb(0,0,1,0.8))
      abline(h=0.95,lty=2,col="gray")
      abline(v=0.95,lty=2,col="gray")
    }
  }
}

cell_hit_names2 <- cbind(rownames(cell_hit),cell_hit_names)

