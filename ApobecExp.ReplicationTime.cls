/// 
/// Class for analysis of distribution of APOBEC mutations 
/// relative to replication time
/// 
Class ApobecExp.ReplicationTime [ Abstract ]
{

/// Function returns value of replication timing
/// for particular cell line and chromosome position
ClassMethod GetRTval(cellline As %String, chr As %String, pos As %Integer) As %Float
{
	 if '$data(^RT(cellline,chr,pos), list)
	 {
		set startpos = $order(^RT(cellline,chr,pos),-1,list) 
		if (startpos = "") {
			return ""	}
		else	
		{
			if (pos > $listget(list,1))	{
				return ""	}
			else {
				return $listget(list,2)	}
		}
	 } 
	 else {
		 return $listget(list,2) }
}

/// Function returns replication timing bin number
/// for particular cell line and chromosome position
ClassMethod GetRTbin(cellline As %String, chr As %String, pos As %Integer) As %Float
{
	 if '$data(^RT(cellline,chr,pos), list)
	 {
		set startpos = $order(^RT(cellline,chr,pos),-1,list) 
		if (startpos = "") {
			return 0	}
		else	
		{
			if (pos > $listget(list,1))	{
				return 0	}
			else {
				return $listget(list,3)	}
		}
	 } 
	 else {
		 return $listget(list,3) }
}

/// Function returns interval of replication timing measurment 
/// for particular cell line and chromosome position
ClassMethod GetRTInterval(cellline As %String, chr As %String, pos As %Integer) As %Library.List
{
	if '$data(^RT(cellline,chr,pos), list)
	{
		set startpos = $order(^RT(cellline,chr,pos),-1,list)
		if (startpos = "") {
			return "" }
		else
		{
			if (pos > $listget(list,1))	{
				return ""	}
			else 
			{
				set endpos = $listget(^RT(cellline,chr,startpos),1)
				return $listbuild(startpos,endpos)	
			}
		}	
	}
	else 
	{
		set startpos = pos
		set endpos = $listget(^RT(cellline,chr,pos),1)
		return $listbuild(startpos,endpos)
	}
}

/// Function returns number of TCW motifs (potential APOBEC targets)
/// for the particular chromosome and position range
ClassMethod GetTCWinRange(chr As %String, startpos As %Integer, endpos As %Integer) As %Integer
{
	set cnt = 0
	if ($data(^TCW(chr,startpos)) '= 0) {
		set cnt = cnt + 1 }
	set pos = $order(^TCW(chr,startpos))
	while ((pos '= "") && (pos <= endpos))
	{
		set cnt = cnt + 1
		set pos = $order(^TCW(chr,pos))
	}
	
	return cnt
}

/// Function divide replication timing values
/// into seven bins for all cell lines
ClassMethod CreateRTBins()
{
	kill ^RTbins
	for i = 1:1:^cellLinesRT
	{
		set cl = ^cellLinesRT(i)
		write cl, !
		
		kill ^RTval
		set cnt = 0
		set totalTCW = 0
		set tmp1 = 0
		set tmp2 = 0
		
		set chr = $order(^RT(cl,""))
		while(chr '= "")
		{
			write chr, !
			
			set pos = $order(^RT(cl,chr,""))
			while (pos '= "")
			{
				//write pos, !
				set startpos = pos
				set endpos = $listget(^RT(cl,chr,pos),1)
				set val = $listget(^RT(cl,chr,pos),2)
				set TCWcnt = ..GetTCWinRange(chr,startpos,endpos)
				if ($data(^RTval(val)) = 0) {
					set ^RTval(val) = TCWcnt 
					set tmp1 = tmp1 + 1}
				else {
					set ^RTval(val) = ^RTval(val) + TCWcnt 
					set tmp2 = tmp2 + 1}
				set cnt = cnt + 1
				set totalTCW = totalTCW + TCWcnt
				set pos = $order(^RT(cl,chr,pos))
			}

			set chr = $order(^RT(cl,chr))
		} 
		
		set ^RTval = $lb(cnt,totalTCW)
		
		set binsize = totalTCW/7
		write totalTCW, $c(9), binsize, $c(9), tmp1, $c(9), tmp2, !
		
		set val = $order(^RTval(""))
		set curTCW = 0
		set binnum = 1
		set j = 0
		while (val '= "")
		{
			set curTCW = curTCW + ^RTval(val)
			if (curTCW > (binsize*binnum))
			{
				set ^RTbins(cl,val) = binnum
				set binnum = binnum + 1
			}
			set pval = val
			set val = $order(^RTval(val))
			set j = j + 1
		}
		set ^RTbins(cl,pval) = binnum
		write j, $c(9), curTCW, $c(9), !
		
	}
}

/// Function writes replication timing bins 
/// into replication timimng global
ClassMethod SetRTBins()
{
	for i = 1:1:^cellLinesRT
	{
		set cl = ^cellLinesRT(i)
		write cl, !

		set chr = $order(^RT(cl,""))
		while(chr '= "")
		{
			write chr, !
			
			set pos = $order(^RT(cl,chr,""))
			while (pos '= "")
			{
				set val = $listget(^RT(cl,chr,pos),2)
				if ($data(^RTbins(cl,val),binnum) = 0)
				{
					set val = $order(^RTbins(cl,val))
					set binnum = ^RTbins(cl,val)
				}
				set ^RT(cl,chr,pos) = $listupdate(^RT(cl,chr,pos),3,binnum)
								
				set pos = $order(^RT(cl,chr,pos))
			}
			set chr = $order(^RT(cl,chr))
		}
	}
}

/// Function calculates and writes number of TCW motifs
/// for replication timing intervals
ClassMethod GetRTNumberOfTCW()
{
	set TOTALTCW = ^TCW
	
	kill ^RTnumTCW
	for i = 1:1:^cellLinesRT
	{
		set cl = ^cellLinesRT(i)
		write cl, !
		
		for j = 1:1:7
		{
			set ^RTnumTCW(cl,j) = 0
		}
		
		set chr = $order(^RT(cl,""))
		while(chr '= "")
		{
			write chr, !
			
			set pos = $order(^RT(cl,chr,""))
			while (pos '= "")
			{
				set startpos = pos
				set endpos = $listget(^RT(cl,chr,pos),1)
				set binnum = $listget(^RT(cl,chr,pos),3)
				set TCWcnt = ..GetTCWinRange(chr,startpos,endpos)
				
				set ^RTnumTCW(cl,binnum) = ^RTnumTCW(cl,binnum) + TCWcnt

				set pos = $order(^RT(cl,chr,pos))
			}
			set chr = $order(^RT(cl,chr))
		}	
		
		set sum = 0
		for j = 1:1:7
		{
			set sum = sum + ^RTnumTCW(cl,j)
		}	
		set ^RTnumTCW(cl,0) = TOTALTCW - sum
	}
}

/// Function calculates and writes number of bases
/// for replication timing intervals
ClassMethod GetRTNumberOfAll()
{
	set TOTALDNALENGTH = 3095677412
	kill ^RTnumAll
	for i = 1:1:^cellLinesRT
	{
		set cl = ^cellLinesRT(i)
		write cl, !
		
		for j = 1:1:7
		{
			set ^RTnumAll(cl,j) = 0
		}
		
		set chr = $order(^RT(cl,""))
		while(chr '= "")
		{
			write chr, !
			
			set pos = $order(^RT(cl,chr,""))
			while (pos '= "")
			{
				set startpos = pos
				set endpos = $listget(^RT(cl,chr,pos),1)
				set binnum = $listget(^RT(cl,chr,pos),3)
				set len = endpos - startpos
	
				set ^RTnumAll(cl,binnum) = ^RTnumAll(cl,binnum) + len

				set pos = $order(^RT(cl,chr,pos))
			}
			set chr = $order(^RT(cl,chr))
		}	
		
		set sum = 0
		for j = 1:1:7
		{
			set sum = sum + ^RTnumAll(cl,j)
		}	
		set ^RTnumAll(cl,0) = TOTALDNALENGTH - sum

	}
}

/// Function calculates and saves number of APOBEC mutations
/// for each cancer type, sample and replication timing bin
ClassMethod MutationRTAnalysis()
{
	 kill ^RTresultsAPOBEC
	 kill ^RTresultsOther
	 set ctype = $order(^mutation(""))
	 while (ctype '= "")
	 {
		 write ctype, !
		 
		 set cl = ^cancerTypeRTCellLines(ctype)
		 
		 set sample = $order(^mutation(ctype,""))
		 while (sample '= "")
		 {
		  	write sample, !
		  	
		  	for i = 0:1:7 
		  	{
			  	set ^RTresultsAPOBEC(ctype,sample,i) = 0 
			  	set ^RTresultsOther(ctype,sample,i) = 0
			}
		  	
		  	set chr = $order(^mutation(ctype,sample,""))
		  	while (chr '= "")
		  	{
			  	//write chr, !
			  	
			  	set pos = $order(^mutation(ctype,sample,chr,""))
			  	while (pos '= "")
			  	{
				  	//write pos, !
				  	
					set binnum = ..GetRTbin(cl,chr,pos) 	
				  	set isAPOBEC = $listget(^mutation(ctype,sample,chr,pos),3)
				  	if (isAPOBEC = 1) {
					   set ^RTresultsAPOBEC(ctype,sample,binnum) = ^RTresultsAPOBEC(ctype,sample,binnum) + 1 }
					else {
					   set ^RTresultsOther(ctype,sample,binnum) = ^RTresultsOther(ctype,sample,binnum) + 1 }
				  	
				  	set pos = $order(^mutation(ctype,sample,chr,pos))
			  	} 
			  	
			  	set chr = $order(^mutation(ctype,sample,chr))	  	
		  	}
		  	
		  	set sample = $order(^mutation(ctype,sample))	 
		 }
		 
		 set ctype = $order(^mutation(ctype))
	 }
}

}
