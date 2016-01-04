resultsDirectory = "/Users/kathryn/Code/HasteRemote/testoutput/"
#folder = "TestBuskePopulationNoBox"
folder = "TestBuskePopulationInBox"
file = "/results_from_time_0/TrackingData.txt"

fullPath = paste(resultsDirectory,folder,file, sep="")
data = read.table(fullPath, as.is=TRUE, header=FALSE, sep="\t");

allCells = unique(data[,2])
print(allCells)

library("RColorBrewer")
rainbow=c()
for(i in seq(1,length(allCells),by=1)){
	rainbow = c(rainbow, rgb(0,floor(255*(i/length(allCells))),100, maxColorValue=255, alpha=255))
}

par(mfrow=c(1,3))
for(dim in c(1,2,3)){

	plotted = FALSE
	proj = c(1,2,3)
	proj = proj[proj!=dim] 

	counter=1
	for(cell in allCells){

		pos = data[data[,2]==cell,3:5]
		nEntries = dim(pos)[1]
		diffs = pos[1:nEntries-1,]-pos[2:nEntries,]
		print(max(sqrt(apply(diffs*diffs,1,sum))))

		projectedPos = pos[,proj]
		times = data[data[,2]==cell,1]

		if(plotted==FALSE){

			if(dim==1){
				plot(projectedPos[,1], projectedPos[,2], type='l', lwd=1, col=rainbow[counter], xlab="y", ylab="z",
					xlim=c(-6*10^(-6), 6*10^(-6)), ylim=c(-6*10^(-6), 6*10^(-6)))
			}
			if(dim==2){
				plot(projectedPos[,1], projectedPos[,2], type='l', lwd=1, col=rainbow[counter], xlab="x", ylab="z",
					xlim=c(-6*10^(-6), 6*10^(-6)), ylim=c(-6*10^(-6), 6*10^(-6)))
			}
			if(dim==3){
				plot(projectedPos[,1], projectedPos[,2], type='l', lwd=1, col=rainbow[counter], xlab="x", ylab="y",
					xlim=c(-6*10^(-6), 6*10^(-6)), ylim=c(-6*10^(-6), 6*10^(-6)))
			}
			plotted = TRUE
		}else{
			points(projectedPos[,1], projectedPos[,2], type='l', lwd=1, col=rainbow[counter])
		}

		counter=counter+1
	}

	points(c(-5*10^(-6), 5*10^(-6)), c(5*10^(-6), 5*10^(-6)), lty=2, col="red3", type='l', lwd=2)
	points(c(-5*10^(-6), 5*10^(-6)), c(-5*10^(-6), -5*10^(-6)), lty=2, col="red3", type='l', lwd=2)
	points(c(5*10^(-6), 5*10^(-6)), c(-5*10^(-6), 5*10^(-6)), lty=2, col="red3", type='l', lwd=2)
	points(c(-5*10^(-6), -5*10^(-6)), c(-5*10^(-6), 5*10^(-6)), lty=2, col="red3", type='l', lwd=2)
}
