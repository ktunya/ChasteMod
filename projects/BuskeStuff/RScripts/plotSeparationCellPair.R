plotCellPairs = function( vectorOfPlotNumbers ){

resultsDir = "/home/kathryn/ChasteClean/testoutput/"
directories = c("TestBuskeCellPairWithAdhesionNoRadial",
                "TestBuskeCellPairWithAdhesionRadialChangesOnly",
                "TestBuskeCellPairWithElasticityNoRadial",
                "TestBuskeCellPairWithElasticityRadialChangesOnly",
                "TestBuskeCellPairWithCompressionNoRadial",
                "TestBuskeCellPairWithCompressionRadialChangesOnly",
                "TestBuskeCellsInIsolationWithCompressionRadialChangesOnly",
                "TestBuskeVsSpringNoRadial",
                "TestBuskeVsSpringNoRadialwCellCellFriction",
                "TestBuskeVsSpringElasticityOnly",
                "TestBuskeVsSpringElasticityPlusCompression",
                "TestBuskeVsSpringElasticityPlusAdhesionPlusCompression")
par(lwd=2)

for( pn in vectorOfPlotNumbers ){

    if(pn == 1){
        postscript("B1.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        data = read.table( paste(resultsDir, directories[1], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        times = unique(data[,1])
        separation = rep(0, length(times))
        posCell0 = data[data[,2]==0, 3:5]*10^6
        posCell1 = data[data[,2]==1, 3:5]*10^6
        for( t in seq(1,length(times)) ){
            separation[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        r0 = data[data[,2]==0, 6][1]*10^6
        r1 = data[data[,2]==1, 6][1]*10^6
        expectedAsymptote = sqrt(max(r0^2, r1^2)-min(r0^2, r1^2));

        plot(times, separation, lwd=2, xlab="Time (s)", ylab=expression(paste("Separation (",mu,"m)",sep="")),
             , type='l', xlim=c(0,600))#main = "Change in cell pair separation under adhesion"
        abline(a=expectedAsymptote, b=0, lwd=2, lty=2, col="red")
        text(x=40, y=expectedAsymptote+0.01, labels="Expected", col="red")
        dev.off()

    }else if(pn ==2){
        postscript("B2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        data = read.table( paste(resultsDir, directories[2], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        times = unique(data[,1])
        radCell0 = data[data[,2]==0, 6]*10^6
        radCell1 = data[data[,2]==1, 6]*10^6
        
        ylabel = expression(paste("Radius (", mu, "m)"))
        plot(times, radCell0, lwd=2, xlab="Time (s)", ylab=ylabel, #main = "Change in cell radii, fixed cell pair under adhesion",
                type='l',  xlim=c(0,4000), ylim = c(0, max(c(radCell0,radCell1))) , col="dodgerblue")
        points(times, radCell1, lwd=2, col='red3', type='l')
        dev.off()


    }else if(pn == 3){
        postscript("B3.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        data = read.table( paste(resultsDir, directories[3], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        times = unique(data[,1])
        separation = rep(0, length(times))
        posCell0 = data[data[,2]==0, 3:5]*10^6
        posCell1 = data[data[,2]==1, 3:5]*10^6
        for( t in seq(1,length(times)) ){
            separation[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        r0 = data[data[,2]==0, 6][1]*10^6
        r1 = data[data[,2]==1, 6][1]*10^6
        expectedAsymptote = r0 + r1;

        plot(times, separation, lwd=2, xlab="Time (s)", ylab=expression(paste("Separation (",mu,"m)",sep="")),
             , type='l', xlim=c(0,1200))#main = "Change in cell pair separation under elasticity"
        abline(a=expectedAsymptote, b=0, lwd=2, lty=2, col="red")
        text(x=80, y=14.9, labels="Expected", col="red")
        dev.off()

    }else if(pn == 4){
        postscript("B4.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        data = read.table( paste(resultsDir, directories[4], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        times = unique(data[,1])
        radCell0 = data[data[,2]==0, 6]*10^6
        radCell1 = data[data[,2]==1, 6]*10^6
        posCell0 = data[data[,2]==0, 3:5]*10^6
        posCell1 = data[data[,2]==1, 3:5]*10^6
        sumRadii = radCell0 + radCell1
        expectedAsymptote = abs(posCell0[1]-posCell1[1]);

        d01 = rep(0, length(times))
        for( t in seq(1,length(times)) ){
            d01[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        plot(times, sumRadii, lwd=2, xlab="Time (s)", ylab=expression(paste("Sum of radii (",mu,"m)",sep="")),
             , type='l', xlim=c(0,1200), ylim=c(9.9,15.1))#main = "Change in sum of radii for a fixed cell pair under elasticity"
        abline(a=10.0, b=0, lwd=2, lty=2, col="red")
        text(x=80, y=10.0+0.08, labels="Expected", col="red")
        dev.off()

    }else if(pn == 5){
        postscript("B5.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        data = read.table( paste(resultsDir, directories[5], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        times = unique(data[,1])
        radCell0 = data[data[,2]==0, 6]*10^6
        radCell1 = data[data[,2]==1, 6]*10^6
        posCell0 = data[data[,2]==0, 3:5]*10^6
        posCell1 = data[data[,2]==1, 3:5]*10^6
        separation = rep(0, length(times))
        for( t in seq(1,length(times)) ){
            separation[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        expectedAsymptote = abs(radCell0 + radCell1);

        plot(times, separation, lwd=2, xlab="Time (s)", ylab=expression(paste("Separation (",mu,"m)",sep="")),
             , type='l', xlim=c(0,4000),ylim=c(10.5,15.1))#main = "Change in separation of a cell pair under compression"
        abline(a=expectedAsymptote, b=0, lwd=2, lty=2, col="red")
        text(x=250, y=15.10, labels="Expected", col="red")
        dev.off()

    }else if(pn == 6){
        
        postscript("B6.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        data = read.table( paste(resultsDir, directories[6], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        times = unique(data[,1])
        radCell0 = data[data[,2]==0, 6]*10^6
        radCell1 = data[data[,2]==1, 6]*10^6
        posCell0 = data[data[,2]==0, 3:5]*10^6
        posCell1 = data[data[,2]==1, 3:5]*10^6

        d01 = rep(0,length(posCell0))
        for(t in seq(1,length(d01),by=1)){
            d01[t]=sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        x01=(radCell0^2-radCell1^2+d01^2)/(2*d01)
        x10=(radCell1^2-radCell0^2+d01^2)/(2*d01)
        expr = 12*radCell0^2-2*(radCell0-x01)*(2*radCell0-x01)*(1-radCell0/d01) 
               - (radCell0-x01)^2*(2-(radCell0/d01))
        
        VA0 = (4/3)*pi*radCell0^3 - (pi/3)*(radCell0-x01)^2*(2*radCell0-x01);
        VA1 = (4/3)*pi*radCell1^3 - (pi/3)*(radCell1-x10)^2*(2*radCell1-x10);
        VT0 = (4/3)*pi*5^3
        VT1 = (4/3)*pi*10^3

        plot(times, VT0-VA0, lwd=2, xlab="Time (s)", ylab="Difference between target and actual volume",
             , type='l',col='dodgerblue')#main = "Changes in cell volume for a fixed, overlapping cell pair under compression"
        points(times, VT1-VA1, lwd=2, type='l',col='red3')
        legend("topright",legend=c("Cell 1","Cell 2"),fill=c("dodgerblue","red3"),bty='n')
        dev.off() 
    
    }else if(pn == 7){

        postscript("B7.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        data = read.table( paste(resultsDir, directories[7], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        times = unique(data[,1])
        radCell0 = data[data[,2]==0, 6]*10^6
        radCell1 = data[data[,2]==1, 6]*10^6

        plot(times, radCell0, col="dodgerblue", lwd=2, type="l", ylim = c(4,11),
            , xlab="Time (s)", ylab=expression(paste("Radius (",mu,"m)")))#main="Change in cell radii for two isolated, fixed cells under compression"
        points(times, radCell1, col="red3", lwd=2, type="l")

        legend("topright",legend=c(expression(paste("Initial rad. 5",mu,"m, target rad. 10",mu,"m")),
                        expression(paste("Initial rad. 10",mu,"m, target rad. 5",mu,"m"))),
               fill=c("red3","dodgerblue"),bty='n') 
        dev.off()

    #----------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------
    }else if( pn==8 ){

        dataB = read.table( paste(resultsDir, directories[8], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        timesB = unique(dataB[,1])
        posCell0 = dataB[dataB[,2]==0, 3:5]*10^6
        posCell1 = dataB[dataB[,2]==1, 3:5]*10^6
        sepB = rep(0, length(timesB))
        for( t in seq(1,length(timesB)) ){
            sepB[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        radCell0B = dataB[dataB[,2]==0, 6]*10^6
        radCell1B = dataB[dataB[,2]==1, 6]*10^6

        dataS = read.table( paste(resultsDir, "/Spring/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        timesS = unique(dataS[,1])
        posCell0 = dataS[dataS[,2]==0, 3:5]
        posCell1 = dataS[dataS[,2]==1, 3:5]
        sepS = rep(0, length(timesS))
        for( t in seq(1,length(timesS)) ){
            sepS[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        radCell0S = dataS[dataS[,2]==0, 6]
        radCell1S = dataS[dataS[,2]==1, 6]

        #postscript("SpringVsBuske1.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        plot(timesB, sepB, type='l', col="red3", lwd=2, xlab="Time (s)", 
             ylab=expression(paste("Separation (",mu,"m)",sep="")), ylim=c(min(c(sepS,sepB)),max(c(sepS,sepB))) )
        points(timesS, sepS, type='l', col="dodgerblue", lwd=2)
        legend("topleft",legend=c("Buske","Spring"),fill=c("red3","dodgerblue"),bty='n')
        #dev.off()
    

    }else if( pn==9 ){
        #directories[9]
        dataB = read.table( paste(resultsDir, directories[9], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        timesB = unique(dataB[,1])
        posCell0 = dataB[dataB[,2]==0, 3:5]*10^6
        posCell1 = dataB[dataB[,2]==1, 3:5]*10^6
        sepB = rep(0, length(timesB))
        for( t in seq(1,length(timesB)) ){
            sepB[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        radCell0B = dataB[dataB[,2]==0, 6]*10^6
        radCell1B = dataB[dataB[,2]==1, 6]*10^6

        dataS = read.table( paste(resultsDir, "Spring/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        timesS = unique(dataS[,1])
        posCell0 = dataS[dataS[,2]==0, 3:5]
        posCell1 = dataS[dataS[,2]==1, 3:5]
        sepS = rep(0, length(timesS))
        for( t in seq(1,length(timesS)) ){
            sepS[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        radCell0S = dataS[dataS[,2]==0, 6]
        radCell1S = dataS[dataS[,2]==1, 6]

        postscript("SpringVsBuske2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        plot(timesB, sepB, type='l', col="red3", lwd=2, xlab="Time (s)", 
             ylab=expression(paste("Separation (",mu,"m)",sep="")),ylim=c(min(c(sepS,sepB)),max(c(sepS,sepB))) )
        points(timesS, sepS, type='l', col="dodgerblue", lwd=2)
        legend("topleft",legend=c("Buske","Spring"),fill=c("red3","dodgerblue"),bty='n')
        dev.off()
    }

    else if( pn==10 ){
        #directories[9]
        dataB = read.table( paste(resultsDir, directories[10], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        timesB = unique(dataB[,1])
        posCell0 = dataB[dataB[,2]==0, 3:5]*10^6
        posCell1 = dataB[dataB[,2]==1, 3:5]*10^6
        sepB = rep(0, length(timesB))
        for( t in seq(1,length(timesB)) ){
            sepB[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        radCell0B = dataB[dataB[,2]==0, 6]*10^6
        radCell1B = dataB[dataB[,2]==1, 6]*10^6

        dataS = read.table( paste(resultsDir, "Spring/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        timesS = unique(dataS[,1])
        posCell0 = dataS[dataS[,2]==0, 3:5]
        posCell1 = dataS[dataS[,2]==1, 3:5]
        sepS = rep(0, length(timesS))
        for( t in seq(1,length(timesS)) ){
            sepS[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        radCell0S = dataS[dataS[,2]==0, 6]
        radCell1S = dataS[dataS[,2]==1, 6]

        postscript("SpringVsBuskeE1.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        plot(timesB, sepB, type='l', col="red3", lwd=2, xlab="Time (s)", 
             ylab=expression(paste("Separation (",mu,"m)",sep="")),ylim=c(min(c(sepS,sepB)),max(c(sepS,sepB))) )
        points(timesS, sepS, type='l', col="dodgerblue", lwd=2)
        legend("topleft",legend=c("Buske","Spring"),fill=c("red3","dodgerblue"),bty='n')
        dev.off()

        postscript("SpringVsBuskeE2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        plot(timesB, radCell0B, type='l', col="red3", lwd=2, xlab="Time (s)", 
             ylab=expression(paste("Radii (",mu,"m)",sep="")),ylim=c(min(c(sepS,sepB)),max(c(sepS,sepB))) )
        points(timesB, radCell1B, type='l', col="red3", lwd=2)
        points(timesS, radCell0S, type='l', col="dodgerblue", lwd=2)
        points(timesS, radCell1S, type='l', col="dodgerblue", lwd=2)
        legend("topleft",legend=c("Buske","Spring"),fill=c("red3","dodgerblue"),bty='n')
        dev.off()
    }

        else if( pn==11 ){
        #directories[9]
        dataB = read.table( paste(resultsDir, directories[11], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        timesB = unique(dataB[,1])
        posCell0 = dataB[dataB[,2]==0, 3:5]*10^6
        posCell1 = dataB[dataB[,2]==1, 3:5]*10^6
        sepB = rep(0, length(timesB))
        for( t in seq(1,length(timesB)) ){
            sepB[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        radCell0B = dataB[dataB[,2]==0, 6]*10^6
        radCell1B = dataB[dataB[,2]==1, 6]*10^6

        dataS = read.table( paste(resultsDir, "Spring/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        timesS = unique(dataS[,1])
        posCell0 = dataS[dataS[,2]==0, 3:5]
        posCell1 = dataS[dataS[,2]==1, 3:5]
        sepS = rep(0, length(timesS))
        for( t in seq(1,length(timesS)) ){
            sepS[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        radCell0S = dataS[dataS[,2]==0, 6]
        radCell1S = dataS[dataS[,2]==1, 6]

        postscript("SpringVsBuskeEC1.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        plot(timesB, sepB, type='l', col="red3", lwd=2, xlab="Time (s)", 
             ylab=expression(paste("Separation (",mu,"m)",sep="")),ylim=c(min(c(sepS,sepB)),max(c(sepS,sepB))) )
        points(timesS, sepS, type='l', col="dodgerblue", lwd=2)
        legend("topleft",legend=c("Buske","Spring"),fill=c("red3","dodgerblue"),bty='n')
        dev.off()

        postscript("SpringVsBuskeEC2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        plot(timesB, radCell0B, type='l', col="red3", lwd=2, xlab="Time (s)", 
             ylab=expression(paste("Radii (",mu,"m)",sep="")),ylim=c(min(c(sepS,sepB)),max(c(sepS,sepB))) )
        points(timesB, radCell1B, type='l', col="red3", lwd=2)
        points(timesS, radCell0S, type='l', col="dodgerblue", lwd=2)
        points(timesS, radCell1S, type='l', col="dodgerblue", lwd=2)
        legend("topleft",legend=c("Buske","Spring"),fill=c("red3","dodgerblue"),bty='n')
        dev.off()
    }

    else if( pn==12 ){
        #directories[9]
        dataB = read.table( paste(resultsDir, directories[12], "/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        timesB = unique(dataB[,1])
        posCell0 = dataB[dataB[,2]==0, 3:5]*10^6
        posCell1 = dataB[dataB[,2]==1, 3:5]*10^6
        sepB = rep(0, length(timesB))
        for( t in seq(1,length(timesB)) ){
            sepB[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        radCell0B = dataB[dataB[,2]==0, 6]*10^6
        radCell1B = dataB[dataB[,2]==1, 6]*10^6

        dataS = read.table( paste(resultsDir, "Spring/results_from_time_0/TrackingData.txt", sep=""), as.is=TRUE, header=FALSE)
        timesS = unique(dataS[,1])
        posCell0 = dataS[dataS[,2]==0, 3:5]
        posCell1 = dataS[dataS[,2]==1, 3:5]
        sepS = rep(0, length(timesS))
        for( t in seq(1,length(timesS)) ){
            sepS[t] = sqrt( sum( (posCell0[t,]-posCell1[t,])*(posCell0[t,]-posCell1[t,]) ) )
        }
        radCell0S = dataS[dataS[,2]==0, 6]
        radCell1S = dataS[dataS[,2]==1, 6]

        postscript("SpringVsBuskeECA1.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        plot(timesB, sepB, type='l', col="red3", lwd=2, xlab="Time (s)", 
             ylab=expression(paste("Separation (",mu,"m)",sep="")),ylim=c(min(c(sepS,sepB)),max(c(sepS,sepB))) )
        points(timesS, sepS, type='l', col="dodgerblue", lwd=2)
        legend("topleft",legend=c("Buske","Spring"),fill=c("red3","dodgerblue"),bty='n')
        dev.off()

        postscript("SpringVsBuskeECA2.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )
        plot(timesB, radCell0B, type='l', col="red3", lwd=2, xlab="Time (s)", 
             ylab=expression(paste("Radii (",mu,"m)",sep="")),ylim=c(min(c(sepS,sepB)),max(c(sepS,sepB))) )
        points(timesB, radCell1B, type='l', col="red3", lwd=2)
        points(timesS, radCell0S, type='l', col="dodgerblue", lwd=2)
        points(timesS, radCell1S, type='l', col="dodgerblue", lwd=2)
        legend("topleft",legend=c("Buske","Spring"),fill=c("red3","dodgerblue"),bty='n')
        dev.off()
    }

}

}