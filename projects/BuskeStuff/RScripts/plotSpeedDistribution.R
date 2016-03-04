library("caTools")
resultsDir = "/home/kathryn/ChasteClean/testoutput/"
window = 8

BuskeDirectories = c("Run6/TestBuskeDistal_eps2592_d0.00010288_k12960_varying_",
				"Run6/TestBuskeDistal_eps2592_d0.00010288_k12960_",
                "Run6/TestBuskeDistal_d0.00010288_k12960_",
                "Run6/TestBuskeDistal_eps2592_d0.00010288_",
                "Run6/TestBuskeDistal_d0.00010288_")

SpringDirectories = c("Run6/TestBuskeDistal_d0.00010288_",
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
dirIndices = c(1, seq(10,20))
PS = FALSE
if(PS){
   postscript("SpeedComparison2.eps", horizontal = TRUE, onefile = TRUE, paper = "special", width=8, height=8)
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

startT = 75 

for(d in dirIndices){

	data = read.table( paste(resultsDir, directories[d], "/results_from_time_0/TrackingData.txt", 
		               sep=""), as.is=TRUE, header=FALSE)

	times = unique(data[(data[,1]>=startT), 1]) 
	cellData = data[(data[,1]>=startT), c(1,2,3,4,5)]

	if(length(times)!=0){

		print(directories[d])

		speeds = c()

		for( t in times ){

			if(t < times[length(times)]){

				cellsT1 = cellData[cellData[,1]==t, c(2,3,4,5)]
				cellsT2 = cellData[cellData[,1]==(t+1), c(2,3,4,5)]
				ids = cellsT1[cellsT1[,1]%in%cellsT2[,1], 1]
				cellsT1 = cellsT1[cellsT1[,1]%in%ids, ]
				cellsT2 = cellsT2[cellsT2[,1]%in%ids, ]

				cellsT1 = cellsT1[order(cellsT1[,1]), c(2,3,4)]
				cellsT2 = cellsT2[order(cellsT2[,1]), c(2,3,4)]
				
				diff = cellsT2 - cellsT1
				speeds = c(speeds, sqrt(rowSums(diff^2)) )

			}
		}

		maxSpeed = max(speeds)
		minSpeed = min(speeds)
		N = length(speeds)
		h = hist(speeds, breaks=seq(minSpeed, maxSpeed, length.out=15), plot=FALSE)

		if(d==dirIndices[1]){
			plot(h$mids, h$density, 
				 lwd=2,
			     type='l', 
			     ylab="Proportion of recordings", 
				 xlab=expression(paste("Average speed (",mu,"m per hour)", sep="")), 
				 col=colors[d],
				 ylim=c(0,0.2))
		}else{
		    points(h$mids, h$density, lwd=2, type='l', col=colors[d])
		}	

	}
}

legend("topleft", fill=colors[dirIndices], legend=directories[dirIndices], bty="n")

if(PS){
	dev.off()
}