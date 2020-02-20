
#crown radius = 0.5*dbh^0.62; dbh in cm and crown radius in m
#crown area = phi*dbh^theta  (dbh in cm, crown area in m2)
phi = round(pi*.5^2,4); theta = 1.24

PA = 10000 #m2; 1 ha 
deltaT = 5 #years; model timestep
dnot = 1 #recruit initial dbh; 1cm
maxT = 200*deltaT #simulation time (1000 years)

mincohortN = 0.001 #minimum cohort size to keep track of

cutT = deltaT #how often to output statistics of the forest (can be in multiples of deltaT)

##########################################
run = function(mainfolder = "",
	vitalratesfile = "PPA_FG5.txt",
	initialconditionsfile = "PPA_initial_fg5_secondary.txt",
	outputfile = "out_secondary_FG5.txt"){
#runs the simulation given the folder and filename inputs
#mainfolder - the location of the input files and where the output will go
#vitalratesfile - the filename of the file in the mainfolder that includes the vital rates by functional group
#		each row in this file should be a functional group 
#		columns: 
#			fg: functional group number
#			fgname: functional group name
#			G1: Diameter growth rate of trees in layer 1 (full sun); G2: Diameter growth rate of trees in layer 2 (shaded by one layer of trees); G3: Diameter growth rate of trees in layer 3 (shaded by two layers of trees); G4: Diameter growth rate of trees in layer 4 (shaded by three layers of trees).  All in units: cm/year.
#			mu1: mortality rate of trees in layer 1; mu2: mortality rate of trees in layer 2; mu3: mortality rate of trees in layer 3; mu4: mortality rate of trees in layer 4. All in units 1/years.
#			F: new recruitments (at dnot, 1cm dbh). Units individuals/Ha/yr. 
#			wd: wood density (used for biomass calculations, has no effect on simulation). Units g/cm3
#initialconditionsfile - the initital conditions of the forest in the format: 
#		each row in the file is a cohort
#		columns (no names as this will be used as a matrix)
#			(1) diameter in cm
#			(2) number of individuals 
#			(3) functional group number; number should match the vital rates number in vitalratesfile$fg
#outputfile: output destination file of the format: 
#		each row in the file is a time point
#		see calcstats function for explanations of columns 
##########################################
		
	fgdata = read.table(paste(mainfolder,vitalratesfile,sep=""),sep="\t",header=TRUE)
	
	ininitdata = read.table(paste(mainfolder,initialconditionsfile,sep=""),sep="\t",header=FALSE)
	initdata = cbind(ininitdata[,c(1,2)],NaN,ininitdata[,3])
	
	CA = sum(phi*initdata[,1]^theta*initdata[,2]) #total crown area of individuals in the plot
	if(CA<=PA) initdata[,3]=1  #if less than the ground area, no need for CCassign. Everyone is in the canopy. 
	#if greater than the ground area, go through CCassign
	if(CA>PA) initdata = CCassign_manylayers_continuous(initdata,PA,deltaT) 
	#kill plants that fall into a fifth layer immediately
	initdata = initdata[initdata[,3]<5,]
	names(initdata) = NULL; initdata = as.matrix(initdata)

	stats = main_cohorts(fgdata,initdata)
	
	write.table(stats,paste(mainfolder,outputfile,sep=""),sep="\t",row.names=FALSE)
}#end run

#####################
main_cohorts = function(fgdata,initdata=NA){
#runs the simulations, and calculates the statistics for each timestep
#fgdata: dataframe of the vital rates of the functional groups (fg,G1,G2,G3,G4,mu1,mu2,mu3,mu4,F)
#initdata: a matrix of the forest initialization each row is a cohorts, columns: (1) dbh in cm (2) number of individuals (3) crown class (4) functional group
	#if initdata is.na, the simulation will be run from bare ground
#####################

	ndata = matrix(1,nrow=dim(fgdata)[1],ncol=4) #main matrix with columns: (1) cohort diameter (2) number of individuals (3) crown class (4) functional group
	#initialize the matrix with initial cohort of individuals
	ndata[,1] = dnot; ndata[,2] = PA/10000*fgdata$F; ndata[,3] = 1; ndata[,4] = fgdata$fg
	baby1 = ndata; baby1[,1] = baby1[,1]+fgdata$G4*4; baby1[,2] = baby1[,2]*(1-fgdata$mu4)^4
	baby2 = ndata; baby2[,1] = baby2[,1]+fgdata$G4*3; baby2[,2] = baby2[,2]*(1-fgdata$mu4)^3
	baby3 = ndata; baby3[,1] = baby3[,1]+fgdata$G4*2; baby3[,2] = baby3[,2]*(1-fgdata$mu4)^2
	baby4 = ndata; baby4[,1] = baby4[,1]+fgdata$G4*1; baby4[,2] = baby4[,2]*(1-fgdata$mu4)^1
	baby5 = ndata
	babyMatrix = rbind(baby1,baby2,baby3,baby4,baby5)
	babyMatrix[,3] = 4
	
	ndata = babyMatrix
	ndata[,3] = 1

	if(is.na(initdata[[1]])) data = ndata 
	if(!is.na(initdata[[1]])) data = initdata
	
	stats = data.frame(time=0,calcstats(data,fgdata))
	#main loop. Remember data matrix has columns: (1) diameter (2) # of individual (3) crown class (4) functional group
	#crown class = 1 for canopy; crown class = 2 for understory; 3 for 2nd understory; 4 for 3rd understory
	for(t in seq(deltaT,maxT,by=deltaT)){
		for(fg in unique(data[,4])){
			#Step 1: Mortality 
			muV = c(fgdata[fg,seq(7,10)],1) 
			for(i in seq(1,dim(data)[1])[data[,4]==fg]) data[i,2] = data[i,2]*(1-muV[[data[i,3]]])^deltaT 
			data = data[data[,2]>mincohortN,,drop=FALSE]
			#Step 2: Growth by layer
			GV = fgdata[fg,seq(3,6)]
			for(layer in unique(data[data[,4]==fg,,drop=FALSE][,3])) data[data[,4]==fg&data[,3]==layer,1] = data[data[,4]==fg&data[,3]==layer,1]+GV[[layer]]*deltaT
		}
		data = data[data[,2]>mincohortN,,drop=FALSE]
		data = data[data[,1]>0,,drop=FALSE]
		
		#Step 3: Reproduce
		#Note here, reproduction is not a function of this patch itself, but it easily could be with: babies = round(min(PA,sum(data[,1]^theta*phi*data[,2]))*Fnot*deltaT)
		CA = sum(phi*data[,1]^theta*data[,2])
		data = rbind(data,babyMatrix)
		
		#Step 4: Assign crown class
		CA = sum(phi*data[,1]^theta*data[,2]) #total crown area of individuals in the plot
		if(CA<=PA) data[,3]=1  #if less than the ground area, no need for CCassign. Everyone is in the canopy. 
		#if greater than the ground area, go through CCassign
		if(CA>PA) data = CCassign_manylayers_continuous(data,PA,deltaT) 
		#kill plants that fall into a fifth layer immediately
		data = data[data[,3]<5,]

		#Step 5: Record
		if(t>0) if(floor(t/cutT)==t/cutT){  #records data once every cutT timesteps.
			stats = rbind(stats,data.frame(time=t,calcstats(data,fgdata)))
		}
		
	}
	return(stats)
	
}#main_cohorts

######################################
calcstats = function(data,fgdata){
#calculates statistics for the forest for a point in time
#specifically, number of individuals, basal area, and aboveground biomass for each functional group and specific size classes
#Note: assumes specific correspondence of fgname and fg number, does not look it up from fgdata.
#totcarea: total number of canopy layers in the plot
#numcohorts: number of cohorts in the model
#everything else is per hectare. 
#no final number are totals for the functional group. 
#final numbers indicate size classes for sub calculations: 
#	size classes for n (1-6): >1,<5; >5,<20; >=20,<60; >=60; >=5 (cm dbh)
#	size classes for ba and agb (1-3): >5,<20; >=20,<60; >=60 (cm dbh)
######################################	
	out = data.frame(totcarea=NaN,numcohorts=NaN,
		n_Slow=NaN,n_Fast=NaN,n_LLP=NaN,n_SLB=NaN,n_Intermediate=NaN,
		ba_Slow=NaN,ba_Fast=NaN,ba_LLP=NaN,ba_SLB=NaN,ba_Intermediate=NaN,
		agb_Slow=NaN,agb_Fast=NaN,agb_LLP=NaN,agb_SLB=NaN,agb_Intermediate=NaN,
		n_Slow_1=NaN,n_Slow_2=NaN,n_Slow_3=NaN,n_Slow_4=NaN,n_Slow_5=NaN,n_Slow_6=NaN,
		n_Fast_1=NaN,n_Fast_2=NaN,n_Fast_3=NaN,n_Fast_4=NaN,n_Fast_5=NaN,n_Fast_6=NaN,
		n_LLP_1=NaN,n_LLP_2=NaN,n_LLP_3=NaN,n_LLP_4=NaN,n_LLP_5=NaN,n_LLP_6=NaN,
		n_SLB_1=NaN,n_SLB_2=NaN,n_SLB_3=NaN,n_SLB_4=NaN,n_SLB_5=NaN,n_SLB_6=NaN,
		n_Intermediate_1=NaN,n_Intermediate_2=NaN,n_Intermediate_3=NaN,n_Intermediate_4=NaN,n_Intermediate_5=NaN,n_Intermediate_6=NaN,
		ba_Slow_1=NaN,ba_Slow_2=NaN,ba_Slow_3=NaN,
		ba_Fast_1=NaN,ba_Fast_2=NaN,ba_Fast_3=NaN,
		ba_LLP_1=NaN,ba_LLP_2=NaN,ba_LLP_3=NaN,
		ba_SLB_1=NaN,ba_SLB_2=NaN,ba_SLB_3=NaN,
		ba_Intermediate_1=NaN,ba_Intermediate_2=NaN,ba_Intermediate_3=NaN,
		agb_Slow_1=NaN,agb_Slow_2=NaN,agb_Slow_3=NaN,
		agb_Fast_1=NaN,agb_Fast_2=NaN,agb_Fast_3=NaN,
		agb_LLP_1=NaN,agb_LLP_2=NaN,agb_LLP_3=NaN,
		agb_SLB_1=NaN,agb_SLB_2=NaN,agb_SLB_3=NaN,
		agb_Intermediate_1=NaN,agb_Intermediate_2=NaN,agb_Intermediate_3=NaN
	)

	out$totcarea = totcarea = sum(phi*data[,1]^theta*data[,2])/PA 
	data = data[order(data[,1]),]

	out$numcohorts = dim(data)[1]

	bigdata = data[data[,1]>=5,]
	
	for(fg in unique(bigdata[,4])){
		ssdata = bigdata[bigdata[,4]==fg,]
		nbaabg = n_ba_agb(ssdata,fgdata)
		out[,2+fg] = nbaabg$n
		out[,7+fg] = nbaabg$ba
		out[,12+fg] = nbaabg$agb
	}

	for(fg in unique(data[,4])){
		ssdata = data[data[,4]==fg,]
		
		sz1 = ssdata[ssdata[,1]>1,,drop=FALSE]; sz1=sz1[sz1[,1]<5,,drop=FALSE]
		sz2 = ssdata[ssdata[,1]>5,,drop=FALSE]; sz2=sz2[sz2[,1]<20,,drop=FALSE]
		sz3 = ssdata[ssdata[,1]<20,,drop=FALSE]
		sz4 = ssdata[ssdata[,1]>=20,,drop=FALSE]; sz4=sz4[sz4[,1]<60,,drop=FALSE]
		sz5 = ssdata[ssdata[,1]>=60,,drop=FALSE];
		sz6 = ssdata[ssdata[,1]>=5,,drop=FALSE]
		fgstart = seq(1,length(names(out)))[names(out)=="n_Slow_1"]+6*(fg-1)
		out[fgstart] = n_ba_agb(sz1,fgdata)$n
		out[fgstart+1] = n_ba_agb(sz2,fgdata)$n	
		out[fgstart+2] = n_ba_agb(sz3,fgdata)$n
		out[fgstart+3] = n_ba_agb(sz4,fgdata)$n
		out[fgstart+4] = n_ba_agb(sz5,fgdata)$n
		out[fgstart+5] = n_ba_agb(sz6,fgdata)$n

		sz2.1 = ssdata[ssdata[,1]>5,,drop=FALSE]; sz2.1=sz2.1[sz2.1[,1]<20,,drop=FALSE]
		sz2.2 = ssdata[ssdata[,1]>=20,,drop=FALSE]; sz2.2=sz2.2[sz2.2[,1]<60,,drop=FALSE]
		sz2.3 = ssdata[ssdata[,1]>=60,,drop=FALSE];
		
		fgstart = seq(1,length(names(out)))[names(out)=="ba_Slow_1"]+3*(fg-1)
		out[fgstart] = n_ba_agb(sz2.1,fgdata)$ba
		out[fgstart+1] = n_ba_agb(sz2.2,fgdata)$ba	
		out[fgstart+2] = n_ba_agb(sz2.3,fgdata)$ba

		fgstart = seq(1,length(names(out)))[names(out)=="agb_Slow_1"]+3*(fg-1)
		out[fgstart] = n_ba_agb(sz2.1,fgdata)$agb
		out[fgstart+1] = n_ba_agb(sz2.2,fgdata)$agb	
		out[fgstart+2] = n_ba_agb(sz2.3,fgdata)$agb
	}
	return(out)
}# end calcstats

######################################
n_ba_agb = function(godata,fgdata){
#calculates the number, basal area, and aboveground biomass for a given subset of data
#units: 
#		n individuals
#		ba basal area, m2
#		agb aboveground biomass, Mg
######################################	
	godata = godata[!is.na(godata[,4]),,drop=FALSE]

	n = sum(godata[,2])
	ba = sum((godata[,1]/200)^2*pi*godata[,2])
	
	#agb calculation from Chave et al. 2005. Materials and Methods of supplementary text. 
	agb = 0
	if(length(unique(godata[,4]))==1){
		wd = fgdata[fgdata$fg==unique(godata[,4]),]$wd
		agb = sum(wd*exp(-1.499+2.148*log(godata[,1])+0.207*log(godata[,1])^2-0.0281*log(godata[,1])^3)/1000*godata[,2])
	}

	return(list(n=n,ba=ba,agb=agb))
}# end n_ba_agb


######################################
CCassign_manylayers_continuous = function(data,PA,deltaT){
# main_cohorts function 
# assigns crown class based on the PPA assumption
# assumes all individuals have the same crown area and height allometries
# assumes that CAtot>PA 
# works for any number of layers 
######################################

	CAv = phi*data[,1]^theta
	data = data[order(CAv,decreasing=TRUE),]
	CAv = CAv[order(CAv,decreasing=TRUE)]
	cohortCAv = CAv*data[,2]
	cacaV = cumsum(cohortCAv)

	data[,3] = 1
	for(layers in seq(1,floor(sum(cohortCAv)/PA))){
		#make a vector cacaV where the ith entry is the crown area of the i'th cohort plus all cohorts with individuals of greater diameter.
		CAv = phi*data[,1]^theta
		data = data[order(CAv,decreasing=TRUE),]
		CAv = CAv[order(CAv,decreasing=TRUE)]
		cohortCAv = CAv*data[,2]
		cacaV = cumsum(cohortCAv)

		#pull out individuals that are definitely in the canopy and understory
		und = data[cacaV>PA*layers,,drop=FALSE]
		can = data[cacaV<=PA*layers,,drop=FALSE]
		
		#split the first cohort in the understory to fill the leftover open canopy space
		canCA = max(0,sum(phi*can[,1]^theta*can[,2]))
		tosplit = und[1,,drop=FALSE]
		opencan = layers*PA-canCA
		splitind_incan = opencan/(phi*tosplit[1,1]^theta)
		und[1,2] = und[1,2] - splitind_incan
		tosplit[,2] = splitind_incan
		tosplit[,3] = layers
		can = rbind(can,tosplit)
		if(max(can[,3])!=layers) print("error")
		und[,3]=layers+1

		#piece the data back together
		data = rbind(can,und)
		
		data = data[data[,2]>mincohortN,,drop=FALSE]
	}

	return(data)
}# end CCassign_manylayers_continuous