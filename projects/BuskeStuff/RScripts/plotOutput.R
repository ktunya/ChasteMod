library("caTools")
resultsDir = "/home/kathryn/ChasteClean/testoutput/"
window = 8

BuskeDirectories = c("Run6/TestBuskeDistal_eps2592_d0.00010288_k12960_varying_",
				"Run6/TestBuskeDistal_eps2592_d0.00010288_k12960_",
                "Run6/TestBuskeDistal_d0.00010288_k12960_",
                "Run6/TestBuskeDistal_eps2592_d0.00010288_",
                "Run6/TestBuskeDistal_d0.00010288_")
SpringDirectories = c("TestBuskeDistal_eps2592_d0.00010288_k12960_",
	                  "Run3/TestSpringDistal_mu15_alpha5_",   #2-9
				      "Run3/TestSpringDistal_mu150_alpha5_",
				      "Run5/TestSpringDistal_mu225_alpha5_",
				      "Run5/TestSpringDistal_mu300_alpha5_",
				      "Run5/TestSpringDistal_mu375_alpha5_",
				      "Run5/TestSpringDistal_mu450_alpha5_",
				      "Run5/TestSpringDistal_mu600_alpha5_",
				      "Run5/TestSpringDistal_mu750_alpha5_",
                      "Run3/TestSpringDistal_mu15_alpha100_", #10-20
                      "Run3/TestSpringDistal_mu75_alpha100_",
                      "Run3/TestSpringDistal_mu150_alpha100_",
                      "Run3/TestSpringDistal_mu225_alpha100_",
                      "Run3/TestSpringDistal_mu300_alpha100_",
                      "Run3/TestSpringDistal_mu375_alpha100_",
                      "Run4/TestSpringDistal_mu450_alpha100_",
                      "Run4/TestSpringDistal_mu600_alpha100_",
                      "Run4/TestSpringDistal_mu750_alpha100_",
                      "Run4/TestSpringDistal_mu900_alpha100_",
                      "Run4/TestSpringDistal_mu1050_alpha100_",
                      "Run4/TestSpringDistal_mu150_alpha60_", #21-28
                      "Run4/TestSpringDistal_mu225_alpha60_",
                      "Run4/TestSpringDistal_mu300_alpha60_",
                      "Run4/TestSpringDistal_mu375_alpha60_",
                      "Run4/TestSpringDistal_mu450_alpha60_",
                      "Run5/TestSpringDistal_mu600_alpha60_",
                      "Run5/TestSpringDistal_mu750_alpha60_",
                      "Run5/TestSpringDistal_mu900_alpha60_")
directories = SpringDirectories
dirIndices = seq(4,9)
PS = TRUE
if(PS){
   postscript("OutputComparison1.eps", horizontal = TRUE, onefile = TRUE, paper = "special", width=8, height=8)
}

runClassSizes = c(8,11,8)
colors = c(rgb(0,0,0,maxColorValue=255))
for(c in seq(1,3)){
	for(val in seq(0, 255, length.out=runClassSizes[c])){
		if(c==1){
			colors = c(colors,rgb(val, 0, 255, maxColorValue=255))
		}
		if(c==2){
			colors = c(colors,rgb(255, val, 0, maxColorValue=255))
		}
		if(c==3){
			colors = c(colors,rgb(0, 255, val, maxColorValue=255))
		}
	}	
}
for(d in dirIndices){

	tryCatch({

	data =  read.table( paste(resultsDir, directories[d], "/PlaneDeathsData.txt", sep=""),
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
	
	trimLen = window%/%2

	if(d == dirIndices[1]){
		plot(boxMidpoints[trimLen:(length(boxMidpoints)-trimLen)], runmean(counts, window, endrule="trim"), lwd=2, type='l', xlab="Time (hours)", 
			 ylab="Niche output (n cells)", col=colors[d], ylim=c(5,55), xlim=c(0,50))
	}else{
	    points(boxMidpoints[trimLen:(length(boxMidpoints)-trimLen)], runmean(counts, window, endrule="trim"), lwd=2, type='l', col=colors[d])
	}	


	}, error = function(e){})
}

legend("topright", fill=colors[dirIndices], legend=directories[dirIndices], bty="n")

if(PS){
	dev.off()
}