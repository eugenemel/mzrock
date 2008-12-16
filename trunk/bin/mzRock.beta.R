library('randomForest'); #for training and classification
library("waveslim");

#create private env. for mzrock
if (! exists("mzrock")) { mzrock=new.env(); }

#function to load SRM -> Compound spreadsheet
loadCompoundNames <- function() { 
	srmfile = "/mnt/mzRockServer/DataSets/SRM2.csv";
	if(length(grep("linux", paste(version$os)))==0) srmfile = "C:/mzRock/DataSets/SRM2.csv"; 
	compounds=read.csv(srmfile);
	compounds=compounds[!is.na(compounds$precursorIntensity),];
	return(compounds);
}

#load mzCSV files from a given project directory
loadSamples <- function(projectDir=".") {

	csvfiles <- list.files(paste(projectDir),'.mzCSV$');
	mzrock$samples = data.frame(); 
	len=length(csvfiles);

	for (i in 1:len ) { 
		filename = paste(projectDir, "/" , csvfiles[i],sep='');
		sample = read.csv(filename); 

		if(nrow(sample) == 0) {
			cat("Ignoring empty input file", filename, "\n");
			next;
		}
	
		#sample name from filename
		sampleName = filename; 
		sampleName=gsub(".*/", "",sampleName); 
		sampleName=gsub(".mzCSV", "", sampleName);

		#make sure data is clean
		sample$rt = as.numeric(paste(sample$rt));
		sample$intensity = as.numeric(paste(sample$intensity));
		sample = sample[ !is.na(sample$rt) & !is.na(sample$intensity), ];	
		cat(sampleName,nrow(sample),"\n");

		#within each scan, if there are multiple fragments, keep only one fragments with highest intensity
		delRows=c();
		for(i in which(diff(sample$rt)==0)) { 
			if (sample$intensity[i] < sample$intensity[i+1] ) { r=i; } else { r=i+1; }
			delRows=c(delRows,r);
			#cat(i,i+1,r,sample$rt[i],sample$rt[i+1],sample$intensity[i],sample$intensity[i+1],"\n");
		}
		
		#save sample into mzrock env.
		sample$sample <- sampleName;
		if (!is.null(delRows)) {
			mzrock$samples=rbind(mzrock$samples, sample[-delRows, c('rt','intensity','srmid','sample')]);
		} else{
			mzrock$samples=rbind(mzrock$samples, sample[, c('rt','intensity','srmid','sample')]);
		}
		cat("Import", sampleName, "Total number of scans=", nrow(mzrock$samples), "\n");
	}
	
	#renumber rows such that they are sequential
	rownames(mzrock$samples) = 1:nrow(mzrock$samples);
	mzrock$sampleNames = unique(paste(mzrock$samples$sample));

	#set sample colors.. rainbow
	mzrock$sampleColors = rainbow(length(mzrock$sampleNames),start=0.2,end=0.9,gamma=1.5);

	#blank samples are always red
	mzrock$blankSamples = grep("blan", mzrock$sampleNames,ignore.case=TRUE);
	mzrock$sampleColors[ mzrock$blankSamples ] = "red"; 

	#align samples
	alignSamples();

	#construct EICs for each SRM
	mzrock$peaks=calculatePeaks(plot=FALSE);

	#renumber peaks such that they are sequential
	rownames(mzrock$peaks) = 1:nrow(mzrock$peaks);

	#create group names by appending "Group" in front of the groupRank
	mzrock$peaks$group = paste(mzrock$peaks$srmid, "Group", mzrock$peaks$groupRank);
}


calculatePeaks = function(srmIds=NA,plot=FALSE,maxGroupNum=5, group0cutoff=.5) { 
	cat("calculating Peaks\n");

	#get list of srmids 
	if (is.na(srmIds)) { srmIds = unique(mzrock$samples$srmid); }
	else { srmIds = as.factor(srmIds); }

	allpeaks=data.frame();  #data frame with peak location information

	#for each srm
	for(sid in levels(srmIds)) {
		cat(sid,"\n");
		samples = subset(mzrock$samples, srmid == sid ); 	#get all observations for this SRM
		if(nrow(samples)==0) next;

		eics = tapply(1:nrow(samples),samples$sample,list);		#group observations by sample name. ie make EICs
		expectedRt = getUserSpecifiedRt(sid) * 60;				#expected retention of the peak for the compound 


		srmPeaks=data.frame();		#data frame with peaks for a selected SRM
		for(eic in eics ) {										
			if ( length(eic) < 1 ) next;
			peaks=enumeratePeaks(samples[ eic, ],plot=plot);		
			peaks=mergeOverlapingPeaks(peaks);
			if(plot==TRUE) points(peaks$rt,peaks$intensity,col=2,cex=3,pch=5); 
			if(nrow(peaks)>0) srmPeaks=rbind(srmPeaks,peaks);
		}

		if (nrow(srmPeaks)>0) {
			blankSamples = samples[ grep("blank",samples$sample, ignore.case=TRUE), ];	
			
			#group peaks
			srmPeaks=groupPeaks(srmPeaks,plot=plot);
			#srmPeaks = groupPeaksNew(srmPeaks,samples);
			
			#fill in missing peaks
			srmPeaks=fillMissingPeaks(srmPeaks);
			srmPeaks$blankBaseLine=1;
			if (nrow(blankSamples)>0) {
				srmPeaks$blankBaseLine = getBlankBaseLine(srmPeaks,blankSamples);
			}

			#rank groups by intensity
			groupSums = tapply(srmPeaks$fracArea,srmPeaks$group, sum,na.rm=TRUE);
			groupRank=rank(groupSums); 
			groupRank=max(groupRank)+1-groupRank;	#reverse ranking .. bigger groups first
			srmPeaks = srmPeaks[ groupRank[paste(srmPeaks$group)] < maxGroupNum, ];	#remove groups that are too low in rank
			srmPeaks$groupRank = srmPeaks$group= groupRank[paste(srmPeaks$group)];	#set groupNumbers

			#make group with closest rt to expected rt group #0
			if (!is.na(expectedRt)) {	

				groupRt =   tapply(srmPeaks$rt,srmPeaks$group,mean,na.rm=TRUE);
				groupSums = tapply(srmPeaks$fracArea,srmPeaks$group, sum,na.rm=TRUE);
				maxGrp = NA;
				print(groupRt);

				for (g in rownames(groupSums) ) {
					g = paste(g);
					if ( abs(groupRt[g]-expectedRt) > group0cutoff*60 ) next;

					if (!is.na(maxGrp)) {
							if( groupSums[g] > groupSums[maxGrp] ) maxg = g;
					} else {
							maxGrp = g;
					}
				}

					#cat(expectedRt, maxGrp, "\n");

				if ( !is.na(maxGrp) ) {
					srmPeaks$groupRank[ srmPeaks$group == maxGrp ] = 0;
				}
			}
			allpeaks=rbind(allpeaks,srmPeaks);
			#cat(sid, nrow(srmPeaks), nrow(allpeaks), "\n");
			if(plot==TRUE) { points(srmPeaks$rt/60,srmPeaks$intensity,pch=srmPeaks$group,cex=3); readline("Press enter"); }
		}
	}
	return(allpeaks);
}

getBlankBaseLine = function(peaks,blankSamples) {
	baseLine=c();
	for( i in 1:nrow(peaks) ) {
		baseLine[i] = median(blankSamples$intensity[blankSamples$rt >= (peaks$rtmin[i]-20) & blankSamples$rt <= (peaks$rtmax[i]+20)],na.rm=TRUE);
	}
	baseLine[ is.na(baseLine) ] = median(blankSamples$intensity,na.rm=TRUE);
	return(baseLine);
}

fillMissingPeaks=function(peaks) {
		T=table(peaks$sample,peaks$group);	
		if(!any(T==0)) return(peaks);

		groupNames = colnames(T);
		sampleNames = rownames(T);

		srmid = peaks$srmid[1];
		csample=getCompound(srmid);

		missingPeaks=data.frame();
		for( sid in sampleNames) {
			if(!any(T[sid,]==0)) next;
			sample = subset(csample, sample == sid);

			peakPositions=c();
			peakGroupIds=c();
			for( gid in groupNames) {
				if(T[sid,gid]==1) next;
				gpeaks = subset(peaks, group==gid); 
				grpRtMin = mean(gpeaks$rtmin)-10;
				grpRtMax = mean(gpeaks$rtmax)+10;
				pos=0;
				for(i in 1:nrow(sample)) {
					if (sample$rt[i]<grpRtMin) next;	if (sample$rt[i]>grpRtMax) break;
					if (pos==0)  { pos=i; next; }
					if (sample$intensity[pos] < sample$intensity[i] )  pos=i; 
				}
				if (pos != 0 ) {
					peakPositions = c(peakPositions,pos);
					peakGroupIds  = c(peakGroupIds,gid);
				}
				#cat(grpRtMin,grpRtMax,sample$rt[pos],"\n");
			}

			mpeaks =  peakDetailedInfo(sample,peakPositions,noNoiseThresh=0);
			#print(mpeaks);
			#points(mpeaks$rt/60,mpeaks$intensity,col=2,cex=3);
			mpeaks$group = as.numeric(peakGroupIds);
			missingPeaks = rbind(missingPeaks,mpeaks);
		}
		return(rbind(peaks,missingPeaks));
}

groupPeaks = function(cpeaks,plot=FALSE) {

	cat("groupingPeaks\n");
	cpeaks = groupPeaksSimple(cpeaks);
	allpeaksfinale = reduceGroups(cpeaks);

	if ( plot==TRUE) {
		ngroups = max(rpeaks$group);
		plotSamples(sid1,showPeaks=FALSE);
		points(rpeaks$rt/60,rpeaks$intensity,cex=2,pch=19,col=rainbow(ngroups)[rpeaks$group]);
		print(table(rpeaks$sample,rpeaks$group));
	}
	return(allpeaksfinale);
}

groupPeaksSimple = function(pks,rtwin=60) {
	npeaks = nrow(pks); 
	if(npeaks < 2 ) { pks$group=1; return(pks); }
	groups=rep(0,npeaks); 
	groups=cutree(hclust(dist(data.frame(pks$rt,groups)),method="average"),h=rtwin);
	pks$group = groups;
	return(pks);
}

groupPeaksNew = function(allpeaks, EIC) { 
	npeaks = nrow(allpeaks); 
	if(npeaks < 2 ) { allpeaks$group=1; return(allpeaks); }
	allpeaks$group = 1;

	spline = smooth.spline(EIC$rt, EIC$intensity)
	groupBounds = findBounds(spline)
	#baseLine	= findBaseLine(EIC);
	plot(spline,type="l"); points(spline$x[ groupBounds$pos ], spline$y[ groupBounds$pos], type="p"); lines(baseLine, col=2);
	groupIds = rep(0, nrow(allpeaks));
	overlapScore = rep(0, nrow(allpeaks));

	if ( nrow(groupBounds) == 0 ) return(allpeaks); 

	for ( i in 1:nrow(groupBounds)) {
		rtmin = spline$x[ groupBounds$minpos[i] ];
		rtmax = spline$x[ groupBounds$maxpos[i] ];
		for ( j in 1:nrow(allpeaks) ) { 
			overlap = checkOverlap(rtmin,rtmax,allpeaks$rtmin[j], allpeaks$rtmax[j]);
			if ( overlap > overlapScore[j] ) {
					groupIds[j] = i;
					overlapScore[j] = overlap;
			}
		}
	}

	allpeaks$group = groupIds;
	return(reduceGroups(allpeaks));
}

reduceGroups = function(cpeaks) { 
	groupNames =  paste(unique(cpeaks$group));
	sampleNames = paste(unique(cpeaks$sample));

	allpeaksfinale=data.frame(); 
	rpeaks=data.frame();	
	for (gid in groupNames) {
		for(sid in sampleNames) {
			peaks = subset(cpeaks, group==gid & sample==sid); 
			if (ncol(peaks)>0) {
				maxpeak=which.max(peaks$intensity);
				rpeaks=rbind(rpeaks,peaks[maxpeak,]);
			}
		}
	}

	allpeaksfinale=rbind(allpeaksfinale,rpeaks);
	return(allpeaksfinale);
}


alignSamples = function() {
	cat("Aligning samples\n")

	minIntensity = 1000;	#align only peaks with intensity at least this high
	minNObs = 20;			#align samples only if they share this many peaksj
	minCor	= 0.9;			#align samples only if groups have high correlations
	gpeaks =  mzrock$samples$intensity > minIntensity;
	N = sum(gpeaks,na.rm=TRUE);

	rtmatrix=tapply( as.numeric(rownames(mzrock$samples))[gpeaks], 
			list(mzrock$samples$sample[gpeaks],mzrock$samples$srmid[gpeaks]), 
			function(x) { s=mzrock$samples[x,]; if(max(s$intensity,na.rm=TRUE)>minIntensity) { s$rt[which.max(s$intensity)];} else { return(NA); }}
	);

	blankSamples = grep("blank", rownames(rtmatrix), ignore.case=TRUE);	
	if ( length(blankSamples) > 0 ) rtmatrix  = rtmatrix[ -blankSamples, ];	#remove blanks from the alignment

	N=nrow(rtmatrix);
	grpMedianRt=apply(rtmatrix,2, median,na.rm=TRUE);

	for(i in 1:N) { 
		sampleName = rownames(rtmatrix)[i];
		nobs = sum(!is.na(grpMedianRt + rtmatrix[i,])); 
		if(nobs < 20) next;
		if(cor(grpMedianRt,rtmatrix[i,],use="complete.obs") < minCor) next;
		fit=lm(grpMedianRt~rtmatrix[i,]); a=coef(fit)[2]; b=coef(fit)[1]; 
		subset=mzrock$samples$sample == sampleName;
		if ( is.na(a) || is.na(b) ) next;
		mzrock$samples$rt[subset] = mzrock$samples$rt[subset]*a+b;
		cat(sampleName,nobs,a,b,"\n");
	}
}

getQuality <- function(peaks=mzrock$peaks) {
	groupIDs = unique(peaks$group);
	allq=data.frame();
	peakIds = c();
	for(groupID in groupIDs) {
			#cat("getQuality: ", groupID, "\n");
			p = subset(peaks, group==groupID);
			peakIds = c(peakIds,rownames(p));
			q=p[, c("fracArea","baseCorrectedIntensity", "baseCorrectedArea", "noNoiseObs" ) ];

			q$signalNoiseRatio = p$signalNoiseRatio
			q$baseLineSlope = abs(p$baseLineSlope)
			q$window = p$rtmax - p$rtmin;
			q$heightWidthRatio = p$intensity/(q$window+1);
			
			#features related to blank
			q$isBlank=rep(0,nrow(p));
			q$isBlank[ grep("blan", p$sample, ignore.case=TRUE) ] = 1;
			q$signalBaseRatio = p$areaTop/p$blankBaseLine

			#blankPeaks = p[ grep("blan", p$sample, ignore.case=TRUE ), ]
			#blankMax=max(blankPeaks$areaTop);
			#blankMax[ is.na(blankMax) ] = 32;
			#q$blankRatio = p$areaTop/blankMax;
			#q$overlapWithBlank = maxOverlapWithBlank(p,blankPeaks);

			#features related to group
			groupWidth=max(p$rtmax)-min(p$rtmin);
			q$windowRatio = (p$rtmax - p$rtmin)/(groupWidth+1);
			q$rtdiff = abs(p$rt-mean(p$rt))/sd(p$rt); 
			q$rtdiff[ is.na(q$rtdiff) ] = -1;

			peakOverlapCount=rep(0,nrow(p));
			for(i in 1:nrow(p) ) { 
				for(j in i:nrow(p) ) {
					if (i==j) next;
					a=p$rtmin[i]; b=p$rtmax[i];
					c=p$rtmin[j]; d=p$rtmax[j];
					nooverlap = d<a  & d<b | b<c & b<d; 
					if(!nooverlap) {
							peakOverlapCount[i] = peakOverlapCount[i]+1;
							peakOverlapCount[j] = peakOverlapCount[j]+1;
					}
				}
			}
			q$overlapCounter = peakOverlapCount;
			allq = rbind(allq,q);
	}
	allq = allq[ order(as.numeric(peakIds)), ];
	return(allq);
}

groupStats=function(cpeaks) {
	grpMedianRt=tapply(cpeaks$rt,cpeaks$group,median);
	grpMaxInt=tapply(cpeaks$intensity,cpeaks$group,max);
	grpGoodness=tapply(cpeaks$posterior.g,cpeaks$group,sum); 
	grpNGoodSamples=tapply(cpeaks$posterior.g,cpeaks$group,function(x){ sum(x>0.5)}); 
	grpNSamples=tapply(cpeaks$sample,cpeaks$group, function(x) { length(unique(x)); } );
	grpBlankMax =max( cpeaks$intensity[ grep("blan", cpeaks$sample, ignore.case=TRUE) ] );
	srmid=tapply(paste(cpeaks$srmid),cpeaks$group, function(x) { x[1] });
	groupid=tapply(paste(cpeaks$group),cpeaks$group, function(x) { x[1] });
	groupRank=tapply(cpeaks$groupRank,cpeaks$group, function(x) { x[1] });
	df=data.frame(srmid,grpMedianRt,grpMaxInt,grpNSamples,grpNGoodSamples,grpGoodness,grpBlankMax,groupRank);
	rownames(df) = groupid;
	return(df);
}

checkOverlap <-  function(a,b,c,d) { #check for overlap between two line segments# return fraction overlap
	if( (a<c && b<c) || (d<a && d<b) )    return(0.0);	#no overalp
	if( (c<=a && b<=d) || (a<=c && d<=b) ) return(1.0);	#100% overlap
	if  (c<=a) return(1/abs((c-b)/(d-a)));
	if  (a<=c) return(abs((c-b)/(d-a)));
}

maxOverlapWithBlank = function(peaks,blankPeaks) {

	overlapCount = rep(0, nrow(peaks));
	if(nrow(blankPeaks)<=0) return(overlapCount);
	if(nrow(peaks)<=0) return(overlapCount);

	for ( j in 1:nrow(blankPeaks) )  {
		for ( i in 1:nrow(peaks) ) {
	  		overlap=checkOverlap(peaks$rtmin[i],peaks$rtmax[i],blankPeaks$rtmin[j],blankPeaks$rtmax[j]);
			overlapCount[i] = max(overlapCount[i],overlap);
		}
	}
	return(overlapCount);
}


detectSlopeChanges <- function(x) { 
	len=length(x);
	if (len < 3) { return( c() ); }
	slope=diff(x); 
	noiseLevel = median(x);

	x[x<0]=0;
	peaks=c();
	j=0; 
	for(i in 1:(len-2)) {
		if(slope[i] > 0 & slope[i+1] < 0 & x[i+1]>noiseLevel) { j = j+1; peaks[j]=i+1; }
	}
	return(peaks);
}

mergeOverlapingPeaks <- function(peaks,maxRtDiff=60){
	if(nrow(peaks)<=1) return(peaks); 
	i=0;
	while(1) {
		i=i+1;
		rtdiff = abs(peaks[i,]$rt - peaks[i+1,]$rt)
		if (is.null(rtdiff) ) next;

		overlap = checkOverlap(peaks$rtmin[i],peaks$rtmax[i],peaks$rtmin[i+1], peaks$rtmax[i+1]);
		if(rtdiff<maxRtDiff | (overlap>0 & rtdiff < maxRtDiff)) {
			if (peaks[i,]$intensity < peaks[i+1,]$intensity) { removePos=i; } else{ removePos=i+1; }
			peaks=peaks[-removePos,];
			i=i-1;
		 }
		if (i+1 >= nrow(peaks)) break;
	}
	return(peaks);
}

getPeaks <- function(srmID=NA,sampleName=NA,groupID=NA) { 
	srmID = paste(srmID);
	if (!is.na(groupID) & !is.na(sampleName)) {
		return(subset(mzrock$peaks,group == groupID & sample==sampleName));
	} else if (!is.na(sampleName)&!is.na(srmID)) { 
		return(subset(mzrock$peaks,srmid==srmID & sample==sampleName));
	} else if (!is.na(groupID)) {
		return(subset(mzrock$peaks,group==groupID));
	} else if (!is.na(srmID)) {
		return(subset(mzrock$peaks,srmid==srmID));
	} else if (!is.na(sampleName)) {
		return(subset(mzrock$peaks,sample==sampleName));
	} else {
		return(mzrock$peaks);
	}
}

findBaseLine = function(EIC) { 
	q=quantile(EIC$intensity,prob=0.7,na.rm=TRUE); 
	filteredPoint = EIC$intensity<q;
	spline = smooth.spline(EIC$rt[filteredPoint], EIC$intensity[filteredPoint]);
	return(spline);
}
	


findBounds = function(spline) { 
	y=spline$y;
	peakPositions = detectSlopeChanges(spline$y);
	N = length(peakPositions);
	len = length(spline$y);

	peaks = data.frame();
	for (i in 1:length(peakPositions)) { 
		pos = peakPositions[i];
		if ( is.null(pos) ) next;
		ii = pos - 1;
		jj = pos + 1;
		if ( ii < 0   ) next;
		if ( jj > len ) next;

		minposi = ii;	#left bound
		minposj = jj;	#right bound
			
		while(ii > 0){ 
			if (y[ii]<=y[minposi]) minposi=ii;
			if (y[ii]>y[pos]) break;
			if (y[ii]<y[pos]*0.1) break;
			ii=ii-1;
		}

		#peak area calculation.. right bound
		while(jj < len ){ 
			if (y[jj]<=y[minposj]) minposj=jj;
			if (y[jj]>y[pos]) break;
			if (y[jj]<y[pos]*0.1) break;
			jj=jj+1;  
		}
#		cat( pos, minposi, minposj, "\n");
		peaks = rbind( peaks, list('pos'=pos, 'minpos'=minposi, 'maxpos'=minposj));
	}
	return(peaks);
}


peakDetailedInfo  = function(sample,peakPosition,noNoiseThresh=1,snRatio=2) {

	peaks  = data.frame();
	if (length(peakPosition)==0) return(peaks);
	len=nrow(sample); if (len < 4 ) return(peaks); 

	N = length(peakPosition);

	peaks = data.frame( 
					'pos'=numeric(N),
					'minpos'=numeric(N),
					'maxpos'=numeric(N),
					'rt'=numeric(N),
					'intensity'=numeric(N),
					'areaTop'=numeric(N),
					'area'=numeric(N),
					'rtmin'=numeric(N),
					'rtmax'=numeric(N),
					'baseCorrectedArea'=numeric(N),
					'baseCorrectedIntensity'=numeric(N),
					'noNoiseObs'=numeric(N),
					'fracArea'=numeric(N),
					'signalNoiseRatio'=numeric(N),
					'noiseLevel'=numeric(N),
					'baseLineSlope'=numeric(N)
	);


	q50=quantile(sample$intensity,prob=0.5,na.rm=TRUE); 
	noiseLevel=median(sample$intensity[ sample$intensity < q50 ]);
	areaTotal = sum(sample$intensity);

	for (i in 1:length(peakPosition)) { 
			pos = peakPosition[i];
			if ( is.null(pos) ) next;
	
			#get position with maximum intensity 
			a=pos-5; b=pos+5; if (b>len) b=len; if(a<=0)a=1; 
			
			#smoother sometimes shifts location of the peak, look in the window around smoothed line for a real peak
			maxpos = which.max(sample$intensity[a:b]); maxpos=a+(maxpos-1);
			if(maxpos>len) maxpos=len; if (maxpos<=0) maxpos=1; 
			#cat(pos,maxpos,sample$intensity[pos],sample$intensity[maxpos],"\n");

			#peak area calculation.. left bound
			minposi=maxpos-1; if (minposi<=0) minposi=1;
			ii = maxpos; 
			#cat(maxpos,sample$intensity[maxpos],noiseLevel*snRatio,"\n");
			while(ii > 0){ 
				if (sample$intensity[ii]<=sample$intensity[minposi]) minposi=ii;
				if (sample$intensity[ii]>sample$intensity[maxpos]) break;
				if (sample$intensity[ii]<=noiseLevel*snRatio) break; 
				ii=ii-1;
			}
			if (minposi<=0) minposi=1;
	
			#peak area calculation.. right bound
			minposj=maxpos+1; if (minposj>=len) minposj=len;
			jj = maxpos;
			while(jj < len ){ 
				if (sample$intensity[jj]<=sample$intensity[minposj]) minposj=jj;
				if (sample$intensity[jj]>sample$intensity[maxpos]) break;
				if (sample$intensity[jj]<=noiseLevel*snRatio) break; 
				jj=jj+1;  
			}
			if (minposj >= len) { minposj=len; }

			ii=minposi; 
			jj=minposj;

			a=maxpos-1; b=maxpos+1; if(a<0)a=1; if (b>len) b=len;

			#subract area und
			baseCorrected = baselineCorrection(sample$intensity[ii:jj]);

			#if( abs((sample$intensity[jj]-sample$intensity[ii])/(jj-ii))>0.1 ) {
			#	plot(seq(ii-2,jj+2),sample$intensity[(ii-2):(jj+2)],type="l"); 
			#	points(maxpos, sample$intensity[maxpos],col=3,cex=3);
			#	lines(c(ii,jj),sample$intensity[c(ii,jj)],type="l",col=2); 
			#}

			peaks[i, ] = c( maxpos,
						ii,
						jj,
						sample$rt[maxpos],
						sample$intensity[maxpos],
						mean(sample$intensity[a:b]),	
						sum(sample$intensity[ii:jj]),
						sample$rt[ii],
						sample$rt[jj],
						baseCorrected$area,
						baseCorrected$intensity,
						baseCorrected$noNoiseObs,
						baseCorrected$area/areaTotal,
						baseCorrected$signalNoiseRatio,
						noiseLevel,
						(sample$intensity[jj]-sample$intensity[ii])/(jj-ii)
			);
	}
	if ( nrow(peaks)>0) {
		 peaks$srmid= sample$srmid[1]
		 peaks$sample= sample$sample[1]
	}
	return(peaks);
}

baselineCorrection = function(y) {
	ii=1; jj=length(y); 
	x=seq(ii-1,jj-1); 
	slope=(y[jj]-y[ii])/(jj-ii); 
	yp=slope*x+y[ii]; 	#y[ii] = y-intersept
	diff=y-yp; 

	#area
	area = sum(diff[diff>0]);
	
	#intensity
	peakpos = which.max(y); 
	intensityBase=slope*peakpos+y[ii];
	intensity = y[peakpos]-intensityBase;
	
	#noNoiseObs
	noNoiseObs = sum(diff>0);
	
	#signalNoiseRatio
	signalNoiseRatio = intensity/(abs(intensityBase)+1)

	#cat(area,intensity,"\n");
	return(list('area'=area, 'intensity'=intensity, 'noNoiseObs'=noNoiseObs, 'signalNoiseRatio'=signalNoiseRatio));
}


peak.thresh=function(wc, value)
{
    wc.shrink <- wc
    for (i in names(wc)) {
	 if (is.na(i)) next; 
         wci <- wc[[i]]; 
         wc.shrink[[i]] <- wci * (abs(wci) > value);
    }
    wc.shrink
}


enumeratePeaks<- function(sample,plot=FALSE) {
	
	peaks = data.frame();
	len=nrow(sample); if (len < 4 ) return (peaks); 

	d=modwt(sample$intensity,wf="mb4",n.levels=4);	#wavelet transform
    	factor <- quantile(abs(d[["d1"]]),seq(0,1,0.01))
	maxpos = which.max(sample$intensity);

	i=95;
	while(i>0) {  #create smoothed intensity vector, such that  peak with maximum intensity is maintained
		dtmp=peak.thresh(d,factor[i])
		spline=imodwt(dtmp);				
		#plot(sample$intensity,type="l");lines(spline,col=2);
		peakPosition = detectSlopeChanges(spline); 
		npeaks=length(peakPosition);
		#cat(i, npeaks,factor[i], "\n")
		if (maxpos == which.max(spline)) break;
		if (npeaks > 20) break;
		i=i-1;
	}
	
	peaks=peakDetailedInfo(sample,peakPosition);
		
	if (plot) { 
		plot(sample$rt,sample$intensity,type="l"); 
		points(sample$rt,spline,type="l",col=2); 

		q50=quantile(spline,prob=0.5,na.rm=TRUE); 
		noiseLevel=median(spline[ spline < q50 ]);

		#simple average smoother, guassian in shape
		lines(c(0,10000000),c(noiseLevel,noiseLevel),col=3,lwd=3); 

		for(i in 1:nrow(peaks)) {
			peak=peaks[i,];
			points(peak$rt,peak$intensity,type='p',col=i,cex=1.5,pch=19);
		}
	}
	return(peaks);
}

getCompound <- function(srmId,sampleName=NA) {
	if  ( !is.na(sampleName) )  {
		return( subset(mzrock$samples,srmid == srmId & sample == sampleName));
	} else {
		return (subset( mzrock$samples, srmid == srmId));
	}
}

validateSamples <-function() {
	groupIDs = unique(mzrock$peaks$group);
	#groupSums=tapply(mzrock$peaks$intensity, mzrock$peaks$group, sum, na.rm=TRUE);
	#groupIDs = names(groupSums[order(groupSums,decreasing=TRUE)]);
	mzrock$peaks$valid = NA;


	for(groupID in groupIDs) { 
		cat(groupID,"\n");
		if( sum(!is.na(mzrock$peaks$valid))>1000 ) { 
			validated = which(!is.na(mzrock$peaks$valid))
			validated= data.frame( peakid=validated, valid=mzrock$peaks$valid[validated]);
			newmodel=classifySamples(Q,validated); 
		}
		peaks = getPeaks(groupID=groupID); 
		peaksQ = Q[ rownames(peaks), ]
		peakIds = rownames(peaks);
		srmId = paste(peaks$srmid[1]);
		compoundSamples = getCompound(srmId);
		rtmin = min(peaks$rt)-360;
		rtmax = max(peaks$rt)+360;
		blanks=peaks[ grep("blan", peaks$sample, ignore.case=TRUE), ]; 

		if ( exists("newmodel")) {
			peakGoodness = predict(newmodel$model,peaksQ,type="prob");
			print(peakGoodness);
			dotSize=3*peakGoodness[,2];
			#dotSize= predict(newmodel$model,peaksQ)$posterior.g * 3;
		} else {
			dotSize=2;
		}

		 valid   = rep(NA, nrow(peaks));
		 for ( j in 1:nrow(peaks)) {
			peak = peaks[j, ];
			q=     peaksQ[j, ];
			par(mfrow=c(1,1));
			plotSamples(srmId,groupID=groupID,showPeaks=TRUE,rtmin=rtmin,rtmax=rtmax,plotTitle=paste(groupID),showBaseLine=TRUE);
			points(peaks$rt/60,peaks$intensity,pch=19,col=3,cex=dotSize);
			points(blanks$rt/60,blanks$intensity,cex=2,pch=19,col=2);
			points(peak$rt/60,peak$intensity,cex=3,pch=3,col=2,lwd=3);
			legend(rtmax/60-2,max(peaks$intensity)*1.2,sprintf("%-15.10s %3.2f",names(q),q[1,]),cex=.9);

			validAnswer=FALSE;
			while ( ! validAnswer ) {
				ANSWER = readline("\nIgnore[i], Good[g] or Bad[g] Peak? [g,b,G,B,i,I,Q]>"); # prompt	
				ANSWER=gsub('\n','',ANSWER);
				test = grep('[GBgbiIQ]',ANSWER);
				if ( length(test) == 1 ) { validAnswer=TRUE; }
			}

			if ( ANSWER == 'Q' ) return();

			if ( ANSWER == 'B' || ANSWER == 'G' || ANSWER == 'I' ) {
				for (k in j:nrow(peaks)) valid[k] = tolower(ANSWER);
				break;
			} else {
				valid[j] = tolower(ANSWER);
			}
			print(valid);
		 }
		valid[ peaksQ$isBlank == 1 ] = 'b'
		mzrock$peaks$valid[ as.numeric(peakIds)] = valid;
	}
	par(mfrow=c(1,1));
}

classifySamples<-function(Q,validatedList,testSet=NA,plotMisClassified=FALSE) {

	V=validatedList[ validatedList$valid == 'g' | validatedList$valid == 'b', ];
	T= Q[ V$peakid, ]
	T$valid<-as.factor(paste(V$valid));


	#model=NaiveBayes(valid ~ . , data=T);	
	#model=lda(valid ~ . , data=T);
	model=randomForest(valid ~ ., data=T,ntree=500,replace=FALSE,proximity=FALSE);
	#model=rpart(valid ~ ., data=T,method="poisson");

	if (!is.na(testSet)) {
		testSet=testSet[ testSet$valid == 'g' | testSet$valid == 'b', ];
		W = Q[ testSet$peakid, ];
		W$valid<-as.factor(paste(testSet$valid));
		#pred = predict(model,W)$class;	#lda
		pred = predict(model,W);	#RandomForest
		#pred =ifelse(pred[,'g']>0.5,'g','b');	#RPART

		tp = sum(W$valid == 'g' & pred == 'g');
		fn = sum(W$valid == 'g' & pred == 'b');
		fp = sum(W$valid == 'b' & pred == 'g');
		tn = sum(W$valid == 'b' & pred == 'b');
		print(table(W$valid,pred));
		
		if ( plotMisClassified == TRUE ) {
			par(ask=TRUE);
			for (i in 1:nrow(testSet) ) {
				if ( paste(testSet$valid[i]) != paste(pred[i]) ) {
					peak = mzrock$peaks[ testSet$peakid[i], ];
					print(peak);
					plotSamples(paste(peak$srmid),groupID=paste(peak$group), plotTitle=paste(testSet$peakid[i], testSet$valid[i], pred[i])) ;
					points(peak$rt/60,peak$intensity,col=2,cex=2,pch=19);
				}
			}
			par(ask=FALSE);
		}
		return(list(model=model,tp=tp,fp=fp,tn=tn,fn=fn));
	}
	return(list(model=model));
}

wavespline <- function(y,thresh=NA,n.levels=4,wf="mb4") { 
	#x=c(y,rep(0,2^ceiling(log2(length(y)))-length(y))); #increase vector size to power of 2 length

	if (is.na(thresh)) { thresh=max(y)*1.2; if (thresh>5000) thresh=5000; }
	d=modwt(y,wf=wf,n.levels=n.levels);	#wavelet transform
	d=manual.thresh(d,value=thresh);

	spline=imodwt(d);						#inverse wavelet to denoised signal
	#spline=spline[1:length(y)];			#back to original vector size
	return(spline);
}


plotCompounds <- function() {
	srmIds = as.vector(unique(mzrock$samples$srmid));
	par(ask=TRUE);
	for (i in 1:length(srmIds) ) plotSamples(srmIds[i]);
	par(ask=FALSE);
}

plotStackedPlots<-function(srmId,compoundSamples) {
	sampleNames = unique(compoundSamples$sample);
	N=length(sampleNames);

	if (N == 1)		   { layout=c(1,1); }
	else if ( N <= 2 ) { layout=c(1,2); }
	else if ( N <= 4 ) { layout=c(2,2); }
	else if ( N <= 6 ) { layout=c(2,3); }
	else if ( N <= 9 ) { layout=c(3,3); }
	else if ( N <= 12 ) { layout=c(3,4); }
	else if ( N <= 16 ) { layout=c(4,4); }
	else { layout=c(5,5); }

	par(mfrow=layout,cex.axis=0.7, cex.main=0.9,mar=c(5,2,2,2));
	for (i in 1:length(sampleNames)) { 
		plotSamples(srmId, compoundSamples=compoundSamples,
						   sampleName=sampleNames[i], 
						   xlab=sampleNames[i],
						   showBestGuessRt=FALSE, showExpectedRt=TRUE); 
	}
	par(mfrow=c(1,1));
}


plotCategoryGroup = function(srmids,X) { 
	ncompounds = length(srmids);
	par(mfrow=c(5,5),mar=c(2,2,2,2)); 
	for(srmid in srmids){ 
		values=X[rownames(X)==paste(srmid), ];
		compoundName = guessCompoundName(paste(srmid));
		if ( length(values) > 0 ) barplot(as.double(values),main=compoundName, col=mzrock$sampleColors); 
	}
	par(mfrow=c(1,1));
}

#for selected groups write out pdf report
writePredictionPDF <- function(prefix,keepGroupsList) {

	#set up sample colors
	mzrock$sampleNames = unique(paste(mzrock$samples$sample));
	mzrock$sampleColors = rainbow(length(mzrock$sampleNames),start=0.2,end=0.9);
	mzrock$sampleColors[ mzrock$blankSamples ] = "red"; 

	filename = paste(prefix, "predictionReport.pdf", sep="");

 	pdf(filename,width=11,height=8);  # Bryson's fault
	for (gid in keepGroupsList) {
			gpeaks = getPeaks(groupID=gid);
			srmid = paste(gpeaks$srmid[1]);
			compoundName =  guessCompoundName(srmid);
			compoundSamples = getCompound(srmid);
			sampleNames = unique(paste(compoundSamples$sample));

			layout(matrix(c(1,2,1,1), 2, 2, byrow = TRUE)); 

			#centered plot
			rtmin=min(gpeaks$rt)-200; rtmax=max(gpeaks$rt)+350;
			plotSamples(srmid,compoundSamples=compoundSamples,
					showPeaks=TRUE,
					showLegend=TRUE,showBestGuessRt=FALSE,
					plotTitle=paste(compoundName, gid),
					rtmin=rtmin,rtmax=rtmax,groupID=gid
			);
			
			#upper right conner plot
			colors = mzrock$sampleColors[ mzrock$sampleNames %in% sampleNames ];
			bars = rep(0,length(sampleNames)); names(bars)=sampleNames;
			z=tapply(gpeaks$intensity-gpeaks$noiseLevel,gpeaks$sample,max);
			bars[rownames(z)] = z; ymax=max(z,na.rm=TRUE)*1.2; if (ymax<32) ymax=32; 
			barplot(bars,col=colors,xaxt="n",las=1,ylim=c(0,ymax),cex.axis=0.6);

			#black bars with probabilites
			bars = rep(0,length(sampleNames)); names(bars)=sampleNames;
			z=tapply((1-gpeaks$posterior.g),gpeaks$sample,max)*z;
			bars[rownames(z)] = z;
			barplot(z,width=0.3,space=3,col="black",add=TRUE,axes=FALSE,axisnames=FALSE);
			layout(matrix(c(1,1,1,1), 2, 2, byrow = TRUE)); 
	}
	dev.off();
}

#for selected groups write out quality and groups vs samples report
#prefix = directory name where spreadsheets will be written
#Q = quality matrix computed via getQuality
#keepGroupsList = subset list of groups for which report will be written
writeSpreadSheets <- function(prefix,Q,keepGroupsList) {

	#keep only peaks in keep group
	Q = subset(Q, group %in% keepGroupsList);

	#create piviot table
	X=tapply(Q$areaTop, list(Q$group,Q$sample), max);
	Xcompound = tapply(Q$compound, Q$group, function(x){paste(x[1])});
	X = cbind(Xcompound,X);	
	gstats = groupStats(mzrock$peaks);
	X = cbind(gstats[rownames(X), c("groupRank","grpNGoodSamples", "grpMedianRt")], X);		
	
	filename1 = paste(prefix, "highestPeaks.csv", sep="");
	filename2 = paste(prefix, "qualityControl.csv",sep="");

	write.csv(X, filename1,row.names=FALSE);
	write.csv(Q[, c("group", "compound", "sample", "posterior.g", "rt", "areaTop", 
							"baseCorrectedIntensity", "baseCorrectedArea", "noNoiseObs", 
							"isBlank", "groupRank", "rtmin", "rtmax") ],
				file=filename2,
				row.names=FALSE);
}



writeReport <- function(model) { 

	#calculate peak features and various quality features
	Q = getQuality(mzrock$peaks);
	Q$areaTop = mzrock$peaks$areaTop;
	Q$rt = mzrock$peaks$rt;
	Q$rtmin = mzrock$peaks$rtmin;
	Q$rtmax = mzrock$peaks$rtmax;
	Q$srmid = mzrock$peaks$srmid;
	Q$sample = mzrock$peaks$sample;
	Q$group = mzrock$peaks$group;
	Q$groupRank = mzrock$peaks$groupRank;

	compoundNames=c(); 
	for(i in 1:nrow(Q)) {
		compoundNames[i] = guessCompoundName(paste(Q$srmid[i]));
	}
	Q$compound=compoundNames;

	#asign probability score to peaks
	P=predict(model,Q,type="prob");	#predict goodness of peaks
	Q$posterior.g = mzrock$peaks$posterior.g = P[,  'g'];
#	Q$posterior.g = mzrock$peaks$posterior.g = P$posterior[, 'g'];

	#create list of peaks groups we want to keep
	keepGroupsList=c(); 
	for( sid in unique(paste(mzrock$peaks$srmid))) {
		gstats = groupStats(getPeaks(srmID=sid));

		#we want to make sure we keep at least one group for each compound
		#if there is user specified retention time, group we must keep is 0
		#else it is group 1

		keptGroupCount=0;
		for ( grp in rownames(gstats) ) {
			group = gstats[ grp, ];
			#always keep group 0 and any group where there is at least one good peak

			if ( group$groupRank == 0 || group$grpNGoodSamples > 1 ) {
				keepGroupsList=c(keepGroupsList,grp); 
				keptGroupCount=keptGroupCount+1;
			}
		}

		#if we didn't find any groups.. still take at least one 
		if ( keptGroupCount == 0  & nrow(gstats) > 0 ) {  
		     keepGroupsList = c( keepGroupsList, rownames(gstats)[1]);
		     keptGroupCount=keptGroupCount+1;
		}
		#cat(sid, keptGroupCount, "\n");
	}

	#crete directory for reports
	dir.create("reports",showWarnings=FALSE);

	#prefix name of the working directory to the name of output file
	dirname= gsub(".*/", "", getwd(), perl=TRUE);  # Bryson's fault
	prefix = paste("reports/", dirname,"-", sep="");

	#write out reports
	writeSpreadSheets(prefix, Q, keepGroupsList );
	writePredictionPDF(prefix , keepGroupsList );

	#strict report for groups where at least third of samples contain good peaks
	frcGoodSamples=tapply(mzrock$peaks$posterior.g,mzrock$peaks$group,function(x){ sum(x>0.5)/length(x); }); 
	keepGroupsList =  names(frcGoodSamples[ frcGoodSamples > 0.3 ]);
	writePredictionPDF( paste(prefix, "strict-", sep=""), keepGroupsList );
	
	#create categories
	#categories=tapply(COMPOUNDS$category,COMPOUNDS$srmId, function(x){paste(x[1])});
	
	#add category names to highestpeak report
	#X=cbind(X,categories[ rownames(X) ])

	#pdf("reports/categories.pdf",width=11,height=8);
	#for(class in rownames(table(categories))) { 
	#		srmIds=categories[categories==class ]; 
	#		plotCategoryGroup(rownames(srmIds),X); 
	#};
	#dev.off();

	
	return("done");
}

guessCompound<-function(srmid) { 
	#library(gsubfn);
	#strapply("+94.100@cid11.00 [46.700-47.700]", "([-+])(.+)@cid(.+) \\[(.+)-(.+)\\]", c, backref = -5)

	polarity=strsplit(srmid,split="")[[1]][1];		#first char is polarity
	x=strsplit(gsub('\\[|@cid|\\]'," ",srmid,perl=TRUE)," ")[[1]];	 #first number is precursor mass
	y=strsplit(x[length(x)],"-");		#last set of numbers are  product range   X.000-Y.000 

	precursorIntensity = abs(as.numeric(x[1]));			#precursor mass is the first number
	if ( length(x) > 1 ) { cid = as.numeric(x[2]); }	#collision energy should be the second number

	if ( length(y)>0 & length(y[1]) > 0) { 
		productMzMin=as.numeric(y[[1]][1]);
		productMzMax=as.numeric(y[[1]][2]);
	} else {
		cat("Warning..", srmid, "is missing product mass range");
		return(data.frame()); 
	}

	hits = paste(polarity) == paste(COMPOUNDS$polarity);
	dist =abs(COMPOUNDS$precursorIntensity-precursorIntensity);
	dist[ is.na(dist) ] = 1000;
	hits[ is.na(hits) ] = FALSE;
	hits = hits & dist <= 0.5;

	if (!is.na(productMzMax) & !is.na(productMzMin) ) {
		 hits = hits &
				COMPOUNDS$productMz < productMzMax & 
		 		COMPOUNDS$productMz > productMzMin 
	}


	if ( any(hits) && sum(hits) == 1 ) {
			return(COMPOUNDS[hits,]);
	} else if (any(hits) && sum(hits) > 1) { 
		mindist = which.min( dist[hits] );
		return(COMPOUNDS[hits,][mindist,]);
	} else {
		return(data.frame());
	}
}

guessCompoundName<-function(srmid) { 
	compound = guessCompound(srmid);
	if (nrow(compound)==1) return(paste(compound$compound));
	return(paste("X",srmid));
}

getUserSpecifiedRt<-function(srmid) { 
	compound = guessCompound(srmid);
	if (nrow(compound)==1) return(compound$expectedRt);
	return(NA);
}



plotSamples<-function(srmId,compoundSamples=NULL,showPeaks=TRUE,showPoints=FALSE,showBestGuessRt=TRUE,showExpectedRt=TRUE,showLegend=FALSE,showGroups=FALSE, plotTitle=NA,sampleName=NA,groupID=NA,rtmin=NA,rtmax=NA,maxIntensity=NA,xlab=NULL,showBaseLine=FALSE) {

			if(is.null(compoundSamples)) compoundSamples = getCompound(srmId);
			sampleNames = unique(compoundSamples$sample);
			minRt = rtmin; if(is.na(minRt) | minRt < 0  | minRt == -Inf) minRt=min(compoundSamples$rt);
			maxRt = rtmax; if(is.na(maxRt) | maxRt < minRt | maxRt == Inf) maxRt=max(compoundSamples$rt);

			#max Y coordinates
			if(is.na(maxIntensity)) maxIntensity = max(compoundSamples$intensity)*1.3;

			expectedRt=getUserSpecifiedRt(srmId);
			plotStarted=FALSE;

			if (is.na(plotTitle)) { plotTitle=paste(srmId) } 
			else if (!is.na(sampleName)) plotTitle=paste(srmId, sampleName);
			if (!is.na(groupID)) { p = getPeaks(groupID=groupID); if(nrow(p)>0)  maxIntensity = max(p$intensity)*1.2; }
			if (is.null(xlab)) xlab = "Time(min)";

			#offset along Y axis
			offsetIntensityConst = maxIntensity*0.1/max(20,length(sampleNames),na.rm=TRUE);
			if (!is.na(sampleName)) offsetIntensityConst=0;

			for(i in length(sampleNames):1) { 
				sampleName2=sampleNames[i];	
				if (!is.na(sampleName) && paste(sampleName) != sampleName2) next;

				s=subset(compoundSamples, sample == sampleName2);
				#d=modwpt(s$intensity,wf="mb4",n.levels=2);	#wavelet transform
				#spline = d$w2.0;
				spline=s$intensity;
				#spline = convolve(s$intensity,c(1,6,1)/8,type="filter");	#smooth spline
				#spline = c(0,spline,0); #padding

				#shift plots so that they don't overlap
				offsetRt=(i-1)*0.05; 
				offsetIntensity=(i-1)*offsetIntensityConst;
				color = mzrock$sampleColors[ mzrock$sampleNames %in% sampleName2 ];
				if ( !is.na(sampleName) && sampleName == sampleName2 ) { lwd=2; } else{ lwd=1; }
				if (length(grep("blank",sampleName2,ignore.case=TRUE))>0 ) { lwd=3; }
				
				if (!plotStarted) {
					plot((s$rt-offsetRt)/60,spline+offsetIntensity,col=color,type='l',
						xlim=c(minRt/60,maxRt/60),ylim=c(0,maxIntensity+maxIntensity*0.2),
						xlab=xlab, ylab="Intensity", main=plotTitle, lwd=lwd, 
						las=1
					);
					plotStarted=TRUE;

					if ( showExpectedRt & !is.na(expectedRt))  {#draw line at expected line
						 points(expectedRt,maxIntensity+1000,col=2,type="h",lwd=2, lty=2);		
					}

				} else {
					lines((s$rt-offsetRt)/60,spline+offsetIntensity,col=color, lwd=lwd );
				}

				if( showPoints) {
					ss=subset(s,intensity>0);
					points((ss$rt-offsetRt)/60,ss$intensity+offsetIntensity,col=color,pch=19,cex=0.5);
				}

				if ( showPeaks == TRUE ) {
					p = getPeaks(srmID=srmId,groupID=groupID,sampleName=sampleName2);
					if (nrow(p)>0) {
					if(is.null(p$posterior.g)) { cex=0.2; } else { cex = 4 * p$posterior.g; }
					points((p$rt-offsetRt)/60,p$intensity+offsetIntensity,cex=cex,col=color,pch=19);
						
					if (showBaseLine) { 
						for (k in 1:nrow(p)) { pi = p[k,];
						  lines((c(s$rt[pi$minpos],s$rt[pi$maxpos])+offsetRt)/60,
						  (c(s$intensity[pi$minpos],s$intensity[pi$maxpos]))+offsetIntensity,col=k,lwd=3);
						  points((pi$rt-offsetRt)/60,pi$intensity+offsetIntensity,cex=cex,col=k,pch=19);
					
						}
					}
				}
			}
		  }

		if (showGroups) { 
			peaks = getPeaks(srmID=srmId,groupID=groupID); i=1;
			groupIDs = unique(peaks$group);
			for (gid in groupIDs) {
				gpeaks = subset(peaks,group==gid);
				rtmin=min(gpeaks$rt); rtmax=max(gpeaks$rt);
				lines(c(rtmin/60,rtmax/60),c(0,0),lwd=3,col=i);
				i = i+1;
			}
		}

		if ( showLegend ) {
			colors = mzrock$sampleColors[ mzrock$sampleNames %in% sampleNames ];
			legend(minRt/60+0.1,maxIntensity+maxIntensity*0.2,cex=1.1,sampleNames,fill=colors,bty='n',ncol=1);
		}
}



#run functions
COMPOUNDS=loadCompoundNames();

