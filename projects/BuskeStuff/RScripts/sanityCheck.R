library("caTools")
resultsDir = "/home/kathryn/ChasteClean/testoutput/"

directories = c("TestBuskeDistal_d0.001_ArtificialBoundary",
				"TestBuskeDistal_eps0.0002_d0.001_ArtificialBoundaryTruncated",
                "TestBuskeDistal_d0.001_ForceBoundary",
                "TestBuskeDistal_d0.001_k1000_ForceBoundaryTruncated",
                "TestBuskeDistal_eps0.0002_d0.001_ForceBoundaryTruncated",
                "TestBuskeDistal_eps0.0002_d0.001_k1000_ForceBoundaryTruncated",
                "TestBuskeDistal_eps0.0002_d0.001_k1000_varying_ForceBasedTruncated")
dirIndices = c(3)
PS = FALSE
if(PS){
   postscript("SanityCheck.eps", horizontal = TRUE, onefile = TRUE, paper = "special", width=8, height=8)
}

colors = c(rgb(200,0,150,maxColorValue=255), rgb(255,0,71,maxColorValue=255), rgb(255,73,0,maxColorValue=255),
           rgb(255,198,0,maxColorValue=255), rgb(210,255,0,maxColorValue=255), rgb(3,163,0,maxColorValue=255),
           rgb(0,255,121,maxColorValue=255), rgb(0,222,255,maxColorValue=255), rgb(48,147,252,maxColorValue=255),
           rgb(13,80,131,maxColorValue=255), rgb(47,23,125,maxColorValue=255), rgb(0,0,0,maxColorValue=255))

par(mfrow=c(1,3))

t = 100

for(d in dirIndices){

	data = read.table( paste(resultsDir, directories[d], "/results_from_time_0/TrackingData.txt", 
		               sep=""), as.is=TRUE, header=FALSE)

	cellsAtT = data[data[,1]==t, c(3,4,5)]

	plot(cellsAtT[,1],cellsAtT[,2], xlab="x", ylab="y", pch='o')
	plot(cellsAtT[,2],cellsAtT[,3], xlab="y", ylab="z", pch='o')
	plot(cellsAtT[,1],cellsAtT[,3], xlab="x", ylab="z", pch='o')

	armRadius = 11.3
	cellRadius = 2.8
	allCells = data[, c(3,4,5)]
	maxSep = 0
	totalSep = 0
	for(i in seq(1,nrow(allCells))){

		if(allCells[i,1] < 0){
			sep = abs((armRadius-cellRadius) - sqrt(allCells[i,1]^2 + allCells[i,2]^2 + allCells[i,3]^2))
		}

	    if(allCells[i,1] >= 0){
	    	sep = abs((armRadius-cellRadius) - sqrt(allCells[i,2]^2 + allCells[i,3]^2))
	    }
	    
	    if(sep > maxSep){
	    	maxSep = sep
	    }

	    totalSep = totalSep + sep
	}
	print(paste("Max distance from intended midline separation:", maxSep))
	print(paste("Mean distance from intended midline separation:", totalSep/nrow(allCells)))	
}

if(PS){
	dev.off()
}