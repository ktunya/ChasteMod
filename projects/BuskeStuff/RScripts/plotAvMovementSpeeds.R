library("caTools")
resultsDir = "/home/kathryn/ChasteClean/testoutput/"
window = 8

directories = c("TestBuskeDistal_d0.001_ArtificialBoundary",
				"TestBuskeDistal_eps0.0002_d0.001_ArtificialBoundaryTruncated",
                "TestBuskeDistal_d0.001_ForceBoundary",
                "TestBuskeDistal_d0.001_k1000_ForceBoundaryTruncated",
                "TestBuskeDistal_eps0.0002_d0.001_ForceBoundaryTruncated",
                "TestBuskeDistal_eps0.0002_d0.001_k1000_ForceBoundaryTruncated",
                "TestBuskeDistal_eps0.0002_d0.001_k1000_varying_ForceBasedTruncated")
dirIndices = seq(1, 7)
PS = FALSE
if(PS){
   postscript("SpeedComparison.eps", horizontal = TRUE, onefile = TRUE, paper = "special", width=8, height=8)
}

colors = c(rgb(200,0,150,maxColorValue=255), rgb(255,0,71,maxColorValue=255), rgb(255,73,0,maxColorValue=255),
           rgb(255,198,0,maxColorValue=255), rgb(210,255,0,maxColorValue=255), rgb(3,163,0,maxColorValue=255),
           rgb(0,255,121,maxColorValue=255), rgb(0,222,255,maxColorValue=255), rgb(48,147,252,maxColorValue=255),
           rgb(13,80,131,maxColorValue=255), rgb(47,23,125,maxColorValue=255), rgb(0,0,0,maxColorValue=255))

for(d in dirIndices){

	data = read.table( paste(resultsDir, directories[d], "/results_from_time_0/TrackingData.txt", 
		               sep=""), as.is=TRUE, header=FALSE)
	
	times = unique(data[,1])
	timeMidpoints = 0.5*(times[1:(length(times)-1)] + times[2:length(times)])
	avSpeeds = rep(0, length(times)-1)

	for( i in seq(1, length(avSpeeds)) ){

		cellsT1 = data[data[,1]==times[i], c(2,3,4,5)]
		cellsT2 = data[data[,1]==times[i+1], c(2,3,4,5)]
		ids = cellsT1[cellsT1[,1]%in%cellsT2[,1], 1]

		for(c in ids){
			p1 = cellsT1[cellsT1[,1]==c,c(2,3,4)]
			p2 = cellsT2[cellsT2[,1]==c,c(2,3,4)]
			avSpeeds[i] = avSpeeds[i] + sqrt(sum((p2-p1)*(p2-p1)))
		}

		avSpeeds[i] = avSpeeds[i]/length(ids)
	}

	trimLen = window %/% 2
	
	if(d==dirIndices[1]){
		plot(timeMidpoints[trimLen:(length(timeMidpoints)-trimLen)], runmean(avSpeeds, window, endrule="trim"), lwd=2, type='l', xlab="Time (hours)", 
			 ylab=expression(paste("Average speed (",mu,"m per hour)", sep="")), 
			 col=colors[d])
	}else{
	    points(timeMidpoints[trimLen:(length(timeMidpoints)-trimLen)], runmean(avSpeeds, window, endrule="trim"), lwd=2, type='l', col=colors[d])
	}	
}

legend("topleft", fill=colors[dirIndices], legend=directories[dirIndices], bty="n")

if(PS){
	dev.off()
}