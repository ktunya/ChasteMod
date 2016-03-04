R1 = 3
R2 = 3
D=1.0288*10^-4

iterations = 10000
interval = c(1000, 15000)
error = c(0,0)

d12 = seq(1, 6, by=10^-3)
BF = -sqrt((R1*R2)/(R1+R2)) * ((R1+R2-d12)^(3/2))/D

mu = interval[1]
SF = mu*(R1+R2)*log(1-(R1+R2-d12)/(R1+R2)) 
error[1] = sum(abs(BF-SF))

mu = interval[2]
SF = mu*(R1+R2)*log(1-(R1+R2-d12)/(R1+R2)) 
error[2] = sum(abs(BF-SF))


for (iter in seq(1,iterations,by=1)){
	
	mu = 0.5*sum(interval)
	SF = mu*(R1+R2)*log(1-(R1+R2-d12)/(R1+R2)) 
	errorCurrent = sum(abs(BF-SF))

	if(errorCurrent<error[1]){
		error[1]= errorCurrent
		interval[1] = mu
	}else if (errorCurrent < error[2]){
		error[2]= errorCurrent
		interval[2] = mu
	}

}

mu = 0.5*sum(interval)
SF = mu*(R1+R2)*log(1-(R1+R2-d12)/(R1+R2)) 
print(paste("R1",R1))
print(paste("R2",R2))
print(paste("mu=",mu))
print(paste("L1 error=",sum(BF-SF)))

postscript("ForceLawComparison.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 5, height = 5 )

plot(d12, BF, type='l', lwd=2, xlab=expression(paste("Separation (", mu, "m)",sep="")), 
	 ylab="Force magnitude", col="red3", ylim=c(min(c(SF,BF)),max(c(BF,SF))))

points(d12, SF, type='l', lwd=2, col="dodgerblue")

legend("topleft",legend=c("Buske","Spring"),fill=c("red3","dodgerblue"),bty='n')

dev.off()