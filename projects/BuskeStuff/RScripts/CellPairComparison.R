Bdirs = c("NoCCFriction", "ElasticOnly","CompIn","RadialIn","Full")
n=4
dir1 = Bdirs[n]
readable1 = "Buske force law"

dir2 = "SpringNewMu"
readable2 = "GLS force law"

SeparationOrRadii = 'R' 
AsPS = TRUE
filename = paste("CellPairRad",n,".eps",sep="")

if(AsPS){
	postscript(filename, horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5)
}

resultsDir = "~/ChasteClean/testoutput/"
data1 = read.table( paste(resultsDir, dir1, "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
times = unique(data1[,1])
data2 = read.table( paste(resultsDir, dir2, "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)

if(SeparationOrRadii == 'S'){

	d1pos1 = data1[ data1[,2]==0, 3:5 ] 
	d1pos2 = data1[ data1[,2]==1, 3:5 ] 

	d2pos1 = data2[ data2[,2]==0, 3:5 ] 
	d2pos2 = data2[ data2[,2]==1, 3:5 ] 

	diff1 = d1pos1 - d1pos2
	diff2 = d2pos1 - d2pos2

	sep1 = sqrt( diff1[,1]*diff1[,1] + diff1[,2]*diff1[,2] + diff1[,3]*diff1[,3] )
    sep2 = sqrt( diff2[,1]*diff2[,1] + diff2[,2]*diff2[,2] + diff2[,3]*diff2[,3] )

	plot(times, sep1, type='l', lwd=2, col="dodgerblue", ylim=c(1,6),
		xlab="Time (hours)", ylab=expression(paste("Separation (",mu,"m)",sep="")))
	points(times, sep2, type='l', lwd=2, col="red3")
	legend("topleft", fill=c("dodgerblue", "red3"), legend=c(readable1, readable2))

	print(paste("L2Norm = ", sqrt(sum((sep1-sep2)*(sep1-sep2))), sep=""))
}

if(SeparationOrRadii == 'R'){
	
	d1rad1 = data1[ data1[,2]==0, 6 ] 
	d1rad2 = data1[ data1[,2]==1, 6 ] 

	d2rad1 = data2[ data2[,2]==0, 6 ] 
	d2rad2 = data2[ data2[,2]==1, 6 ] 

	plot(times, d1rad1, type='l', lwd=2, col="dodgerblue",
		 xlab="Time (hours)", ylab=expression(paste("Cell radii (",mu,"m)",sep="")),)
	points(times, d1rad2, type='l', lwd=2, col="dodgerblue")
	points(times, d2rad1, type='l', lwd=2, col="red3")
	points(times, d2rad2, type='l', lwd=2, col="red3")
	legend("topleft", fill=c("dodgerblue", "red3"), legend=c(dir1, dir2))
}

if(AsPS == TRUE){
	dev.off()
}