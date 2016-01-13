library("caTools")
resultsDir = "/home/kathryn/ChasteClean/testoutput/"
window = 8

directories = c("TESTINGTestBuskeDistal_eps0.0002_d0.001_k1000_")
dirIndices = c(1)
PS = FALSE
if(PS){
   postscript("OutputComparison.eps", horizontal = TRUE, onefile = TRUE, paper = "special", width=8, height=8)
}

colors = c(rgb(200,0,150,maxColorValue=255), rgb(255,0,71,maxColorValue=255), rgb(255,73,0,maxColorValue=255),
           rgb(255,198,0,maxColorValue=255), rgb(210,255,0,maxColorValue=255), rgb(3,163,0,maxColorValue=255),
           rgb(0,255,121,maxColorValue=255), rgb(0,222,255,maxColorValue=255), rgb(48,147,252,maxColorValue=255),
           rgb(13,80,131,maxColorValue=255), rgb(47,23,125,maxColorValue=255), rgb(0,0,0,maxColorValue=255))

for(d in dirIndices){

	data = read.table( paste(resultsDir, directories[d], "/PlaneDeathsData.txt", sep=""),
	                   as.is=TRUE, header=FALSE)
	
	timesMin = data[1,1];
	timesMax = data[nrow(data),1];
	boxWidth = 1;

	boxBoundaries = seq(timesMin, timesMax+boxWidth, boxWidth); 
	boxMidpoints = 0.5*(boxBoundaries[1:(length(boxBoundaries)-1)] 
		              + boxBoundaries[2:length(boxBoundaries)])

	counts = rep(0, length(boxMidpoints))
	for(i in seq(1,length(counts))){
		cellsOutput = data[(data[,1]>boxBoundaries[i]) & (data[,1]<boxBoundaries[i+1]), 2]
		counts[i] = length(cellsOutput)
	}
	
	if(d==dirIndices[1]){
		plot(boxMidpoints, runmean(counts, window, endrule="mean"), lwd=2, type='l', xlab="Time (hours)", 
			 ylab="Niche output (n cells)", col=colors[d])
	}else{
	    points(boxMidpoints, runmean(counts, window, endrule="mean"), lwd=2, type='l', col=colors[d])
	}	
}

legend("topleft", fill=colors[dirIndices], legend=directories[dirIndices])

if(PS){
	dev.off()
}