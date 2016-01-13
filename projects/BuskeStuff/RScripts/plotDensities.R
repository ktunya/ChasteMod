library("caTools")
resultsDir = "/home/kathryn/ChasteClean/testoutput/"
window = 2

directories = c("TESTINGTestBuskeDistal_eps0.0002_d0.001_k1000_")
dirIndices = c(1)
PS = FALSE
if(PS){
   postscript("DensityComparison.eps", horizontal = TRUE, onefile = TRUE, paper = "special", width=8, height=8)
}

third = 'L' #'F' = first, 'M' = mid, 'L' = last

colors = c(rgb(200,0,150,maxColorValue=255), rgb(255,0,71,maxColorValue=255), rgb(255,73,0,maxColorValue=255),
           rgb(255,198,0,maxColorValue=255), rgb(210,255,0,maxColorValue=255), rgb(3,163,0,maxColorValue=255),
           rgb(0,255,121,maxColorValue=255), rgb(0,222,255,maxColorValue=255), rgb(48,147,252,maxColorValue=255),
           rgb(13,80,131,maxColorValue=255), rgb(47,23,125,maxColorValue=255), rgb(0,0,0,maxColorValue=255))

for(d in dirIndices){

	data = read.table( paste(resultsDir, directories[d], "/results_from_time_0/TrackingData.txt", 
		               sep=""), as.is=TRUE, header=FALSE)
	
	times = unique(data[,1])
	armlength = max(data[,3])
	thirds = c(armlength/3, 2*armlength/3)
	nCells = rep(0, length(times))

	for( i in seq(1, length(nCells)) ){

		if(third == 'F'){
			cells = data[ (data[,1]==times[i]) & (data[,3]<thirds[1]), 2]
		}else if(third == 'M'){
			cells = data[ (data[,1]==times[i]) & (data[,3]>thirds[1]) & (data[,3]<thirds[2]), 2]
		}else if(third == 'L'){
			cells = data[ (data[,1]==times[i]) & (data[,3]>thirds[2]), 2]
		}

		nCells[i] = nCells[i] + length(cells)
	}
	
	if(d==dirIndices[1]){
		plot(times, runmean(nCells, window, endrule="mean"), lwd=2, type='l', xlab="Time (hours)", 
			 ylab=paste("Number of cells in third ",third," of the arm.", sep=""), 
			 col=colors[d])
	}else{
	    points(times, runmean(nCells, window, endrule="mean"), lwd=2, type='l', col=colors[d])
	}	
}

legend("topleft", fill=colors[dirIndices], legend=directories[dirIndices])

if(PS){
	dev.off()
}