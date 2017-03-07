/// 
/// Class for the analysis of distribution of APOBEC mutations 
/// relative to the gene expression
/// 
Class ApobecExp.Expression [ Abstract ]
{

/// Set the maximum end of genes among intersecting genes
/// for using in interval tree structure
ClassMethod SetGeneIntervals()
{
	set IndexName = "IndexGene"
	set chr = $order(^ApobecExp.GeneI(IndexName,""))
	while (chr '= "")
	{
		set maxendpos = 0
		set startpos = $order(^ApobecExp.GeneI(IndexName,chr,""))
		while (startpos '= "")
		{
			set endpos = $order(^ApobecExp.GeneI(IndexName,chr,startpos,""))
			while (endpos '= "")
			{
				set ID = $order(^ApobecExp.GeneI(IndexName,chr,startpos,endpos,""))
				if (maxendpos >= startpos) {
					set ^ApobecExp.GeneI(IndexName,chr,startpos,endpos,ID) = maxendpos }
				else {
					set ^ApobecExp.GeneI(IndexName,chr,startpos,endpos,ID) = 0 }	
				if (maxendpos < endpos) {
					set maxendpos = endpos } 
				set endpos = $order(^ApobecExp.GeneI(IndexName,chr,startpos,endpos))
			}
			set startpos = $order(^ApobecExp.GeneI(IndexName,chr,startpos))
		}
		set chr = $order(^ApobecExp.GeneI(IndexName,chr))
	}
}

/// Function returns list of genes covering
/// particular position on particular chromosome
ClassMethod GetGeneList(chr As %String, pos As %Integer) As %Library.List
{
	set ret = ""
	set IndexName = "IndexGene"
	set chrI = " "_$zconvert(chr,"U")
	set startpos = pos
	if ($data(^ApobecExp.GeneI(IndexName,chrI,pos)) '= 0)	{
		set startpos = pos	}
	else {	
		set startpos = $order(^ApobecExp.GeneI(IndexName,chrI,pos),-1) }
	while (startpos '= "")
	{ 
		set endpos = $order(^ApobecExp.GeneI(IndexName,chrI,startpos,""))
		set maxend = 0
		while (endpos '= "")
		{
			set ID = $order(^ApobecExp.GeneI(IndexName,chrI,startpos,endpos,""))
			if (pos <= endpos)	{
				set $list(ret,*+1) = ID }
			set endinterval = ^ApobecExp.GeneI(IndexName,chrI,startpos,endpos,ID)
			if (endinterval > maxend) {
				set maxend = endinterval }
			set endpos = $order(^ApobecExp.GeneI(IndexName,chrI,startpos,endpos))
		}	
		set endinterval = maxend
		if ((endinterval = 0) || (endinterval < pos)) quit
		set startpos = $order(^ApobecExp.GeneI(IndexName,chrI,startpos),-1)
	}
	return ret
}

/// Function returns maximum expression value among genes
/// covering particular position on particular chromosome
/// return: -1 if position in gene but no expression data,
///         -2 if position not in gene
ClassMethod GetMaxExpression(cancer As %String, sample As %String, chr As %String, pos As %Integer)
{
	set genelist = ..GetGeneList(chr,pos)

	if (genelist '= "") {
		set n = $listlength(genelist) }
	else {
		return (-2) }	
	set maxexp = -1
	for j=1:1:n 
	{
		set GeneID = $listget(^ApobecExp.GeneD($listget(genelist,j)),2)
		if ($data(^exp(cancer,sample,GeneID)) '= 0) 
		{
			if (maxexp < ^exp(cancer,sample,GeneID)) {
				set maxexp = ^exp(cancer,sample,GeneID) }
		}
	}

	return (maxexp)
}

/// Function returns expression bin number
/// based on particular expression value
ClassMethod GetExpressionBin(ExpressionValue As %Float) As %Integer
{
	if (ExpressionValue < 0) {
		return ExpressionValue }
	set leftval = $order(^expBins(ExpressionValue),-1)
	set rightval = $order(^expBins(leftval,""))
	if ((ExpressionValue > leftval) && (ExpressionValue <= rightval)) {
		return (^expBins(leftval,rightval)) }
	else {
		return "" }
}

/// Function calculates number of APOBEC and other mutations
/// relative to expression bins
ClassMethod GetExpMutation()
{
	kill ^expResultsAPOBEC
	kill ^expResultsOther
	for i = 1:1:^cancerType
	{
		set cnt = 0
		set cancer = ^cancerType(i)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{
			for j = -2:1:^expBins 
			{
				set ^expResultsAPOBEC(cancer,sample,j) = 0 
				set ^expResultsOther(cancer,sample,j) = 0
			}
			set chr = $order(^mutation(cancer,sample,""))
			while (chr '= "")
			{
				set pos = $order(^mutation(cancer,sample,chr,""))
				while (pos '= "")
				{
				  	set isAPOBEC = $listget(^mutation(cancer,sample,chr,pos),3)

					set maxexp = ..GetMaxExpression(cancer,sample,chr,pos)
					
					set expbin = ..GetExpressionBin(maxexp) 

					if (isAPOBEC = 1) {		
						set ^expResultsAPOBEC(cancer,sample,expbin) = ^expResultsAPOBEC(cancer,sample,expbin) + 1 }
					else {
						set ^expResultsOther(cancer,sample,expbin) = ^expResultsOther(cancer,sample,expbin) + 1 }
							
					set cnt = cnt + 1	
						
					set pos = $order(^mutation(cancer,sample,chr,pos))
				}
				set chr = $order(^mutation(cancer,sample,chr))
			}
			set sample = $order(^mutation(cancer,sample))
		}
		write cancer,"=",cnt,!
	}
}

/// Function calculates and saves number of TCW motifs 
/// (a potential APOBEC target) in genes 
ClassMethod SetTCWinGenes()
{
	kill ^TCWinGenes
	set ^TCWinGenes = 0
	set IndexName = "IndexGene"
	set chr = $order(^ApobecExp.GeneI(IndexName,""))
	while (chr '= "")
	{
		#;write chr,!
		set startpos = $order(^ApobecExp.GeneI(IndexName,chr,""))
		while (startpos '= "")
		{
			#;write startpos,!
			set endpos = $order(^ApobecExp.GeneI(IndexName,chr,startpos,""))
			while (endpos '= "")
			{    
				write chr," ",startpos," ",endpos,!
                set chrlow = $extract(chr,2,*)
                set chrlow = $zconvert(chrlow,"L")   
                set chrlow = ^chrname(chrlow)

				if ($data(^TCW(chrlow,startpos)) '= 0) {
					set pos = startpos }
				else {	
					set pos = $order(^TCW(chrlow,startpos)) }
				while ((pos '= "") && (pos <= endpos))
				{	
					set ^TCWinGenes(chrlow,pos) = 	""
					set ^TCWinGenes = ^TCWinGenes + 1
					set pos = $order(^TCW(chrlow,pos))
				}

			 	set endpos = $order(^ApobecExp.GeneI(IndexName,chr,startpos,endpos))
			}   
			    
			set startpos = $order(^ApobecExp.GeneI(IndexName,chr,startpos))
		}
		set chr = $order(^ApobecExp.GeneI(IndexName,chr))
	}
}

/// Function calculates number of TCW motifs
/// in expression bins
ClassMethod GetExpTCW()
{
	//kill ^expTCW
	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{
			if ($data(^expTCW(cancer,sample,0)) '= 0)
			{
				write cancer," ",sample," Done",!
				set sample = $order(^mutation(cancer,sample))
				continue
			}
			for j = -2:1:7 {
				set ^expTCW(cancer,sample,j) = 0 }
			set chr = $order(^TCWinGenes(""))
			while (chr '= "")
			{
				set cnt = 0
				set pos = $order(^TCWinGenes(chr,""))
				while (pos '= "")
				{
					set maxexp = ..GetMaxExpression(cancer,sample,chr,pos)
					set expbin = ..GetExpressionBin(maxexp) 
					set ^expTCW(cancer,sample,expbin) = ^expTCW(cancer,sample,expbin) + 1
					set cnt = cnt + 1
									
					set pos = $order(^TCWinGenes(chr,pos))
				}
				write chr," ",cnt,!
				set chr = $order(^TCWinGenes(chr))
			}

			write cancer," ",sample,!
			set sample = $order(^mutation(cancer,sample))
		}
	}
}

/// Function calculates number of bases
/// in expression bins
ClassMethod GetExpAll()
{
	//kill ^expAll

	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{
			if ($data(^expAll(cancer,sample,0)) '= 0)
			{
				write cancer," ",sample," Done",!
				set sample = $order(^mutation(cancer,sample))
				continue
			}			
			for j = -2:1:7 {
				set ^expAll(cancer,sample,j) = 0 }	
				write cancer," ",sample,!
				
		    write "Interval part",!				
			set chr = $order(^geneNonoverlap(""))
			while (chr '= "")
			{
				set startpos = $order(^geneNonoverlap(chr,""))
				while (startpos '= "")
				{
					set endpos = $order(^geneNonoverlap(chr,startpos,""))
					while (endpos '= "")
					{
						set pos = startpos
						set maxexp = ..GetMaxExpression(cancer,sample,chr,pos)
						set expbin = ..GetExpressionBin(maxexp) 
						set ^expAll(cancer,sample,expbin) = ^expAll(cancer,sample,expbin) + (endpos - startpos + 1)
				
						set endpos = $order(^geneNonoverlap(chr,startpos,endpos))
					}
					set startpos = $order(^geneNonoverlap(chr,startpos))
				}
				set chr = $order(^geneNonoverlap(chr))
			}
			
			write "Overlapped part",!
			set chr = $order(^geneOverlap(""))
			while (chr '= "")
			{
				set pos = $order(^geneOverlap(chr,""))
				while (pos '= "")
				{
					set maxexp = ..GetMaxExpression(cancer,sample,chr,pos)
					set expbin = ..GetExpressionBin(maxexp) 
					set ^expAll(cancer,sample,expbin) = ^expAll(cancer,sample,expbin) + 1
					set pos = $order(^geneOverlap(chr,pos))
				}
				set chr = $order(^geneOverlap(chr))
			}			
			
			write cancer," ",sample,!					
			set sample = $order(^mutation(cancer,sample))
		}
	}
}

/// Function divide DNA positions covering by genes 
/// into intervals covered by only single gene
/// and intervals covered by multiple genes.
/// Results are saved in globals.
ClassMethod DivideGenesIntoIntervals()
{
	kill ^geneNonoverlap
	set ^geneNonoverlap = 0
	kill ^geneOverlap
	set ^geneOverlap = 0
	set IndexName = "IndexGene"
	set chr = $order(^ApobecExp.GeneI(IndexName,""))
	while (chr '= "")
	{
		set chrlow = $zconvert($extract(chr,2,*),"L")
		set chrlow = ^chrname(chrlow)
		set startpos = $order(^ApobecExp.GeneI(IndexName,chr,""))
		while (startpos '= "")
		{
			set endpos = $order(^ApobecExp.GeneI(IndexName,chr,startpos,""))
			while (endpos '= "")
			{
				set prevn = 0
				for i = startpos:1:endpos
				{
					set genelist = ..GetGeneList(chrlow,i)

					if (genelist '= "") {
						set n = $listlength(genelist) }
					else {
						set n = 0 }		
						
					if (n = 0)
					{
						write "error: no gene, chr=",chrlow,", pos=",i,!
						return
					}		
					elseif (n > 1)
					{
						set ^geneOverlap(chrlow,i) = ""
						set ^geneOverlap = ^geneOverlap + 1
						if (prevn = 1)
						{
							set ^geneNonoverlap(chrlow,istart,iend) = ""
							set ^geneNonoverlap = ^geneNonoverlap + 1
						}
						set prevn = 2
					}		
					elseif(n = 1)
					{
						if (prevn = 1)
						{
							set iend = iend + 1
							if (iend '= i)
							{
								write "error: iend '= i, chr=",chr,", startpos=",startpos,", endpos=",endpos,!
								return 0
							}
							if (i = endpos)
							{
								set ^geneNonoverlap(chrlow,istart,iend) = ""
								set ^geneNonoverlap = ^geneNonoverlap + 1
								set prevn = 0
							}
						}
						else
						{
							set istart = i
							set iend = i
						}
						set prevn = 1
					}
				}

				set endpos = $order(^ApobecExp.GeneI(IndexName,chr,startpos,endpos))
			}
			set startpos = $order(^ApobecExp.GeneI(IndexName,chr,startpos))
		}
		write chr,!
		set chr = $order(^ApobecExp.GeneI(IndexName,chr))
	}
}

}