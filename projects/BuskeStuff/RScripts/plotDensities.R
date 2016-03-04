library("caTools")
resultsDir = "/home/kathryn/ChasteClean/testoutput/"
window = 8

BuskeDirectoriesFB = c("Run6/TestBuskeDistal_eps2592_d0.00010288_k12960_",
                	   "Run6/TestBuskeDistal_d0.00010288_k12960_",
                       "Run6/TestBuskeDistal_eps2592_d0.00010288_",
                       "Run6/TestBuskeDistal_d0.00010288_")
BuskeDirectoriesA = c("Run8/ABTestBuskeDistal_eps2592_d0.00010288_k12960_",
	                  "Run8/ABTestBuskeDistal_d0.00010288_k12960_",
	                  "Run8/ABTestBuskeDistal_eps2592_d0.00010288_",
	                  "Run8/ABTestBuskeDistal_d0.00010288_")
					
BuskeVaryingDirectories = c("Run6/TestBuskeDistal_eps2592_d0.00010288_k12960_varying_",
	                        "Run7/TestBuskeDistal_eps2592_d0.00010288_k12960_varying_A",
	                        "Run7/TestBuskeDistal_d0.00010288_varying_F",
	                        "Run7/TestBuskeDistal_d0.00010288_k12960_varying_A",
	                        "Run7/TestBuskeDistal_eps2592_d0.00010288_k12960_varying_WF")

SpringDirectories = c("Run6/TestBuskeDistal_eps2592_d0.00010288_k12960_",
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
                      "Ru8n5/TestSpringDistal_mu750_alpha60_",
                      "Run5/TestSpringDistal_mu900_alpha60_")

directories = SpringDirectories
dirIndices = c(1,seq(10,20))
#colors = c(rgb(0,0,0,maxColorValue=255), rgb(206,0,6,maxColorValue=255), 
#	       rgb(26,146,14,maxColorValue=255), rgb(14,90,146,maxColorValue=255))

PS = TRUE
if(PS){
   postscript("Density100L.eps", horizontal = TRUE, onefile = TRUE, paper = "special", width=8, height=8)
}
third = 'L' #'F' = first, 'M' = mid, 'L' = last

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

	trimLen = window %/% 2
	
	if(d==dirIndices[1]){
		plot(times[trimLen:(length(times)-trimLen)], runmean(nCells, window, endrule="trim"), lwd=2, type='l', xlab="Time (hours)", 
			 ylab=paste("Number of cells in third ",third," of the arm.", sep=""), 
			 col=colors[d], ylim=c(150, 600), xlim=c(0,200))
	}else{
	    points(times[trimLen:(length(times)-trimLen)], runmean(nCells, window, endrule="trim"), lwd=2, type='l', col=colors[d])
	}	
}

legend("topright", fill=colors[dirIndices], legend=directories[dirIndices], bty="n")

if(PS){
	dev.off()
}