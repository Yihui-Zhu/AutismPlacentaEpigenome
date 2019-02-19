## CpG Features Extraction
## Extract CpG island, shores, shelves and open sea from CpG Islands

#################################################
######### Load the packages #####################
#################################################

# Download CpG island bed file from USCS genome blowser
# Use https://gist.github.com/decodebiology/30cf37a451ce411a4da698e38cc4e782 to generate CpG shores, shelves

# Generate open sea 

# Generate a function
analysis3 <- function(dataset,index){
  subdata <- dataset[which(dataset[,1]==paste("chr",toString(index), sep="")),] 
  subdata = subdata[order(subdata$V2),]
  m = nrow(subdata)
  V = integer(m*2)
  C2 = as.numeric(as.character(subdata[,2]))
  C3 = as.numeric(as.character(subdata[,3]))

if (C2[1] > 4000){
    V[1] = C2[1] - 4000
}else{
  V[1] = 2
}
if (abs(C3[1]-C2[2]) > 8000){
  V[2] = C3[1] + 4000
}else{
  V[2] = C3[1]
}

for (i in 2:(m-1)){
  if (abs(C2[i] - C3[i-1]) > 8000){
    V[2*i-1] = C2[i] - 4000
  }else{
    V[2*i-1] = C3[i-1] + 1
  }
  if (abs(C3[i]-C2[i+1]) > 8000){
    V[2*i] = C3[i] + 4000
  }else{
    V[2*i] = C3[i]
  }
}

if (abs(C2[m] - C3[m-1]) > 8000){
  V[2*m-1] = C2[m] - 4000
}else{
  V[2*m-1] = C3[m-1] + 1
}
  return(V)
} # end function



Q = analysis3(Island_CpG,1)
P = shift_array(Q)
U = reconstruction(P,1)
for (i in 2:22){
  Q = analysis3(Island_CpG,i)
  P = shift_array(Q)
  U = rbind(U,reconstruction(P,i))
} 
Open_sea <- U
write.table(Open_sea, "Open_sea.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

