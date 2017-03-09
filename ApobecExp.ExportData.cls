/// 
/// Class for exporting results of the analysis
/// 
Class ApobecExp.ExportData [ Abstract ]
{

/// Prepare results of the APOBEC mutations density
/// relative replication timing for exporting
ClassMethod PlotRT()
{
	kill ^plotRTAPOBEC
	kill ^plotRTOther
	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set cellline = ^cancerTypeRTCellLines(cancer)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{	
			for j = 1:1:7 {
				set ^plotRTAPOBEC(cancer,sample,j) = ^RTresultsAPOBEC(cancer,sample,j)/^RTnumTCW(cellline,j) }
			set total = 0
			for j = 1:1:7 {
				set total = total + ^plotRTAPOBEC(cancer,sample,j) }
			for j = 1:1:7 {
				set ^plotRTAPOBEC(cancer,sample,j) = ^plotRTAPOBEC(cancer,sample,j) / total }	
			for j = 1:1:7 {
				set ^plotRTOther(cancer,sample,j) = ^RTresultsOther(cancer,sample,j)/^RTnumAll(cellline,j) }
			set total = 0
			for j = 1:1:7 {
				set total = total + ^plotRTOther(cancer,sample,j) }
			for j = 1:1:7 {
				set ^plotRTOther(cancer,sample,j) = ^plotRTOther(cancer,sample,j) / total }	
			write cancer," ",sample,!
			set sample = $order(^mutation(cancer,sample))
		}
	}
}

/// Prepare results of the APOBEC mutations density
/// relative gene expression for exporting
ClassMethod PlotExp()
{
	kill ^plotExpAPOBEC
	kill ^plotExpOther
	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set cellline = ^cancerTypeRTCellLines(cancer)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{	
			for j = 0:1:7 {
				set ^plotExpAPOBEC(cancer,sample,j) = ^expResultsAPOBEC(cancer,sample,j)/^expTCW(cancer,sample,j) }
			set total = 0
			for j = 0:1:7 {
				set total = total + ^plotExpAPOBEC(cancer,sample,j) }
			for j = 0:1:7 {
				set ^plotExpAPOBEC(cancer,sample,j) = ^plotExpAPOBEC(cancer,sample,j) / total }	
			for j = 0:1:7 {
				set ^plotExpOther(cancer,sample,j) = ^expResultsOther(cancer,sample,j)/^expAll(cancer,sample,j) }
			set total = 0
			for j = 0:1:7 {
				set total = total + ^plotExpOther(cancer,sample,j) }
			for j = 0:1:7 {
				set ^plotExpOther(cancer,sample,j) = ^plotExpOther(cancer,sample,j) / total }	
			write cancer," ",sample,!
			set sample = $order(^mutation(cancer,sample))
		}
	}
}

/// Prepare results of the APOBEC mutations density
/// relative replication timing & gene expression for exporting
ClassMethod PlotExpRT()
{
	kill ^plotExpRTAPOBEC
	kill ^plotExpRTOther
	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set cellline = ^cancerTypeRTCellLines(cancer)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{	
			write cancer," ",sample,!
		
			for j = 1:1:7
			{
				for k = 0:1:7
				{
					write "j=",j,", k=",k,!
					write "APOBEC=",^resultsAPOBEC(cancer,sample,j,k),!
					write "TCW=",^expRTTCW(cancer,sample,j,k),!
					
					if (^resultsAPOBEC(cancer,sample,j,k) = 0) {
						set ^plotExpRTAPOBEC(cancer,sample,j,k) = 0 }
					else {
						set ^plotExpRTAPOBEC(cancer,sample,j,k) = ^resultsAPOBEC(cancer,sample,j,k)/^expRTTCW(cancer,sample,j,k) }
				}
			}
			set total = 0
			for j = 1:1:7 
			{
				for k = 0:1:7
				{
					set total = total + ^plotExpRTAPOBEC(cancer,sample,j,k) 
				}
			}
			for j = 1:1:7 
			{
				for k = 0:1:7
				{
					set ^plotExpRTAPOBEC(cancer,sample,j,k) = ^plotExpRTAPOBEC(cancer,sample,j,k) / total
				}
			}
			for j = 1:1:7
			{
				for k = 0:1:7
				{
					write "j=",j,", k=",k,!
					write "Other=",^resultsOther(cancer,sample,j,k),!
					write "All=",^expRTAll(cancer,sample,j,k),!
					
					if (^resultsOther(cancer,sample,j,k) = 0) {
						set ^plotExpRTOther(cancer,sample,j,k) = 0 }
					else {
						set ^plotExpRTOther(cancer,sample,j,k) = ^resultsOther(cancer,sample,j,k)/^expRTAll(cancer,sample,j,k)  }
				}
			}
			set total = 0
			for j = 1:1:7 
			{
				for k = 0:1:7
				{
					set total = total + ^plotExpRTOther(cancer,sample,j,k) 
				}
			}
			for j = 1:1:7 
			{
				for k = 0:1:7
				{
					set ^plotExpRTOther(cancer,sample,j,k) = ^plotExpRTOther(cancer,sample,j,k) / total
				}
			}

			set sample = $order(^mutation(cancer,sample))
		}
}
}

/// Export results of replication timing &
/// gene epression analysis
ClassMethod ExportPlotExpRT(dirname As %String = "/home/mkazanov/ApobecExpResults/")
{
	do $system.Process.SetZEOF(1)
	
	set filename = dirname_"expRTAPOBEC.txt"
	open filename:"WNS"
	use filename	
	
	write "CANCER_TYPE"_$c(9)_"SAMPLE"_$c(9)_"RT_BIN"_$c(9)_"EXP_BIN"_$c(9)_"MUT_DENSITY",!
	
	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{	
			for j = 1:1:7 
			{
				for k = 0:1:7
				{
					write cancer_$c(9)_sample_$c(9)_j_$c(9)_k_$c(9)_^plotExpRTAPOBEC(cancer,sample,j,k),!
				}
			}
			set sample = $order(^mutation(cancer,sample))
		}
	}

    close filename
    
    set filename = dirname_"expRTOther.txt"
	open filename:"WNS"
	use filename	
	
	write "CANCER_TYPE"_$c(9)_"SAMPLE"_$c(9)_"RT_BIN"_$c(9)_"EXP_BIN"_$c(9)_"MUT_DENSITY",!
	
	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{	
			for j = 1:1:7 
			{
				for k = 0:1:7
				{
					write cancer_$c(9)_sample_$c(9)_j_$c(9)_k_$c(9)_^plotExpRTOther(cancer,sample,j,k),!
				}
			}
			set sample = $order(^mutation(cancer,sample))
		}
	}

    close filename
}

/// Export results of replication 
/// timing analysis
ClassMethod ExportPlotRT(dirname As %String = "/home/mkazanov/ApobecExpResults/")
{
	
	do $system.Process.SetZEOF(1)
	
	set filename = dirname_"RTAPOBEC.txt"
	open filename:"WNS"
	use filename
	
	write "CANCER_TYPE"_$c(9)_"SAMPLE"_$c(9)_"RT_BIN"_$c(9)_"MUT_DENSITY",!

	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{	
			for j = 1:1:7 
			{
				write cancer_$c(9)_sample_$c(9)_j_$c(9)_^plotRTAPOBEC(cancer,sample,j),!
			}
			set sample = $order(^mutation(cancer,sample))
		}
	}

    close filename
    
    set filename = dirname_"RTOther.txt"
	open filename:"WNS"
	use filename
	
	write "CANCER_TYPE"_$c(9)_"SAMPLE"_$c(9)_"RT_BIN"_$c(9)_"MUT_DENSITY",!

	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{	
			for j = 1:1:7 
			{
				write cancer_$c(9)_sample_$c(9)_j_$c(9)_^plotRTOther(cancer,sample,j),!
			}
			set sample = $order(^mutation(cancer,sample))
		}
	}

    close filename
}

/// Export results of  gene
/// expression analysis
ClassMethod ExportPlotExp(dirname As %String = "/home/mkazanov/ApobecExpResults/")
{
	
	do $system.Process.SetZEOF(1)
	
	set filename = dirname_"expAPOBEC.txt"
	open filename:"WNS"
	use filename
	
	write "CANCER_TYPE"_$c(9)_"SAMPLE"_$c(9)_"EXP_BIN"_$c(9)_"MUT_DENSITY",!

	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{	
			for j = 0:1:7 
			{
				write cancer_$c(9)_sample_$c(9)_j_$c(9)_^plotExpAPOBEC(cancer,sample,j),!
			}
			set sample = $order(^mutation(cancer,sample))
		}
	}

    close filename
    
    set filename = dirname_"expOther.txt"
	open filename:"WNS"
	use filename
	
	write "CANCER_TYPE"_$c(9)_"SAMPLE"_$c(9)_"EXP_BIN"_$c(9)_"MUT_DENSITY",!

	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set sample = $order(^mutation(cancer,""))
		while (sample '= "")
		{	
			for j = 0:1:7 
			{
				write cancer_$c(9)_sample_$c(9)_j_$c(9)_^plotExpOther(cancer,sample,j),!
			}
			set sample = $order(^mutation(cancer,sample))
		}
	}

    close filename
}

/// Export APOBEC enrichment data
/// 
ClassMethod ExportEnrichment(dirname As %String = "/home/mkazanov/ApobecExpResults/")
{
	set filename = dirname_"sample_enrichent.txt"
	open filename:"WNS"
	use filename
	
	write "CANCER_TYPE"_$c(9)_"SAMPLE"_$c(9)_"APOBEC_ENRICHMENT",!

	for i = 1:1:^cancerType
	{
		set cancer = ^cancerType(i)
		set sample = $order(^sample(cancer,""))
		while (sample '= "")
		{	
			set enrich = $listget(^sample(cancer,sample),1)
			write cancer_$c(9)_sample_$c(9)_enrich,!
			set sample = $order(^mutation(cancer,sample))
		}
	}

	close filename
}

}