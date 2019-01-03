library(png)
#setwd("~/dss-data-cat/")

lexicon <- c("one","two","three","four","five","six","seven","eight","nine")
path = "dss-data-cat/train"
dirs = list.files(path)
hold_images <- matrix(rep(0,28*28),28,28)
hold_labels <- c()
hold_names <- c()

for(name in dirs){
  nmz <- match(name, lexicon)
  pngs <- list.files(paste(path,name,sep = "/"))
  for(ds_png in pngs){
    inter_step <- readPNG(paste(paste(path,name,sep = "/"),ds_png,sep = "/"))
    hold_labels <- append(hold_labels,nmz)
    hold_names <- append(hold_names,ds_png)
    if (length(dim(inter_step)) == 3){ hold_images <- cbind(hold_images,inter_step[,,2])}
    else {hold_images <- cbind(hold_images,inter_step)}}
  }
}

hold_images <- hold_images[,-(1:28)]
hold_names <- hold_names[-length(hold_names)]
hold_labels <- hold_labels[-length(hold_labels)]

uniqLabels <- unique(hold_labels)
uniqLabels <- match(hold_labels,uniqLabels)

resizedImages <- array(hold_images,c(28,28,length(hold_labels)))

for(i in 1:length(hold_labels)){
	inter_matrix <- resizedImages[,,i]
	inter_matrix <- inter_matrix*255
  inter_matrix[inter_matrix != 0] <- inter_matrix[inter_matrix != 0] + 50
  inter_matrix[inter_matrix == 0] <- inter_matrix[inter_matrix == 0] + rnorm(length(inter_matrix[inter_matrix == 0]))
  resizedImages[,,i] <- inter_matrix
}





