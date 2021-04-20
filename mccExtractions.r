#This script is needed as part of the Fig4 script to plot phylogeography maps using Seraphim

mccTreeExtraction = function(mcc_tre, mostRecentSamplingDatum)
	{
		mcc_tab = matrix(nrow=dim(mcc_tre$edge)[1], ncol=11)
		colnames(mcc_tab) = c("node1","node2","length","startLon","startLat","endLon","endLat","startNodeL","endNodeL","startYear","endYear")
		mcc_tab[,c("node1","node2")] = mcc_tre$edge
		mcc_tab[,c("length")] = mcc_tre$edge.length
		for (i in 1:length(mcc_tre$annotations))
			{
				annotations = mcc_tre$annotations[[i]]
				mcc_tab[i,c("endLon","endLat")] = cbind(annotations$location2, annotations$location1) # to be edited
				mcc_tab[i,"endNodeL"] = annotations$height
			}
		for (i in 1:length(mcc_tre$annotations))
			{
				index = which(mcc_tab[,"node2"] == mcc_tab[i,"node1"])
				if (length(index) > 0)
					{
						mcc_tab[i,c("startLon","startLat")] = mcc_tab[index,c("endLon","endLat")]
						mcc_tab[i,"startNodeL"] = mcc_tab[index,"endNodeL"]
					}	else		{
						annotations = mcc_tre$root.annotation
						mcc_tab[i,c("startLon","startLat")] = cbind(annotations$location2, annotations$location1) # to be edited
						mcc_tab[i,"startNodeL"] = annotations$height
					}
				mcc_tab[i,"startYear"] = mostRecentSamplingDatum - mcc_tab[i,"startNodeL"]
				mcc_tab[i,"endYear"] = mostRecentSamplingDatum - mcc_tab[i,"endNodeL"]
			}
		mcc_tab = mcc_tab[order(mcc_tab[,"startYear"],decreasing=F),]
		mcc_tab1 = mcc_tab[1,]; mcc_tab2 = mcc_tab[2:dim(mcc_tab)[1],]
		mcc_tab2 = mcc_tab2[order(mcc_tab2[,"endYear"],decreasing=F),]
		mcc_tab = rbind(mcc_tab1, mcc_tab2); return(mcc_tab)
	}
