/// 
/// Class for the reference data
/// 
Class ApobecExp.DataDictionary [ Abstract ]
{

/// Chromosome names indexed by integer
ClassMethod SetChromosomeNames()
{
	kill ^chr
	set ^chr(1) = "chr1"
	set ^chr(2) = "chr2"
	set ^chr(3) = "chr3"
	set ^chr(4) = "chr4"
	set ^chr(5) = "chr5"
	set ^chr(6) = "chr6"
	set ^chr(7) = "chr7"
	set ^chr(8) = "chr8"
	set ^chr(9) = "chr9"
	set ^chr(10) = "chr10"
	set ^chr(11) = "chr11"
	set ^chr(12) = "chr12"
	set ^chr(13) = "chr13"
	set ^chr(14) = "chr14"
	set ^chr(15) = "chr15"
	set ^chr(16) = "chr16"
	set ^chr(17) = "chr17"
	set ^chr(18) = "chr18"
	set ^chr(19) = "chr19"
	set ^chr(20) = "chr20"
	set ^chr(21) = "chr21"
	set ^chr(22) = "chr22"
	set ^chr(23) = "chrX"
	set ^chr(24) = "chrY"
	set ^chr = 24
}

/// Chromosome names indexed by lowercase chromosome names
ClassMethod SetChromosomeLowerName()
{
	kill ^chrname
	set ^chrname("chr1") = "chr1"
	set ^chrname("chr2") = "chr2"
	set ^chrname("chr3") = "chr3"
	set ^chrname("chr4") = "chr4"
	set ^chrname("chr5") = "chr5"
	set ^chrname("chr6") = "chr6"
	set ^chrname("chr7") = "chr7"
	set ^chrname("chr8") = "chr8"
	set ^chrname("chr9") = "chr9"
	set ^chrname("chr10") = "chr10"
	set ^chrname("chr11") = "chr11"
	set ^chrname("chr12") = "chr12"
	set ^chrname("chr13") = "chr13"
	set ^chrname("chr14") = "chr14"
	set ^chrname("chr15") = "chr15"
	set ^chrname("chr16") = "chr16"
	set ^chrname("chr17") = "chr17"
	set ^chrname("chr18") = "chr18"
	set ^chrname("chr19") = "chr19"
	set ^chrname("chr20") = "chr20"
	set ^chrname("chr21") = "chr21"
	set ^chrname("chr22") = "chr22"
	set ^chrname("chrx") = "chrX"
	set ^chrname("chry") = "chrY"
	set ^chr = 24
}

/// Standard abbreviation of the cancer types
ClassMethod SetCancerTypes()
{
	kill ^cancerType
	set ^cancerType(1) = "BLCA"
	set ^cancerType(2) = "BRCA"
	set ^cancerType(3) = "HNSC"
	set ^cancerType(4) = "LUAD"
	set ^cancerType(5) = "LUSC"
	set ^cancerType = 5
}

/// Cell lines 
ClassMethod SetCellLinesReplication()
{
	kill ^cellLinesRT
	set ^cellLinesRT(1) = "IMR90"
	set ^cellLinesRT(2) = "MCF7"
	set ^cellLinesRT(3) = "NHEK"
	set ^cellLinesRT = 3
}

/// Correspondense of the cancer types and cell lines
ClassMethod SetCancer2RTCellLines()
{
	 set ^cancerTypeRTCellLines("BLCA") = "NHEK"
	 set ^cancerTypeRTCellLines("BRCA") = "MCF7"
	 set ^cancerTypeRTCellLines("HNSC") = "NHEK"
	 set ^cancerTypeRTCellLines("LUAD") = "IMR90"
	 set ^cancerTypeRTCellLines("LUSC") = "IMR90"
}

/// Expression bins bounds in case of 15 bins breakdown
ClassMethod SetExpressionBins14()
{
	set ^expBins14(-9999999,0) = 0
	set ^expBins14(0,10) = 1
	set ^expBins14(10,25) = 2
	set ^expBins14(25,50) = 3
	set ^expBins14(50,100) = 4
	set ^expBins14(100,200) = 5
	set ^expBins14(200,300) = 6
	set ^expBins14(300,400) = 7
	set ^expBins14(400,550) = 8
	set ^expBins14(550,700) = 9	
	set ^expBins14(700,1000) = 10
	set ^expBins14(1000,1400) = 11
	set ^expBins14(1400,2000) = 12
	set ^expBins14(2000,3500) = 13
	set ^expBins14(3500,99999999999) = 14
	set ^expBins14 = 14
}

/// Expression bins bounds in case of 8 bins breakdown
ClassMethod SetExpressionBins()
{
	set ^expBins(-9999999,0) = 0
	set ^expBins(0,25) = 1
	set ^expBins(25,100) = 2
	set ^expBins(100,300) = 3
	set ^expBins(300,550) = 4
	set ^expBins(550,1000) = 5	
	set ^expBins(1000,2000) = 6
	set ^expBins(2000,99999999999) = 7
	set ^expBins = 7
}

}
