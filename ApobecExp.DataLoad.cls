/// 
/// Class for loading all of the required data
/// 
Class ApobecExp.DataLoad [ Abstract ]
{

/// Human genome (DNA sequence) loading 
/// Available from NCBI:
/// ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq/
/// Output global: ^humanDNA
ClassMethod LoadHumanGenome(dirname As %String = "/home/mkazanov/HUMAN")
{
	 kill ^humanDNA
	 do DISABLE^%SYS.NOJRN
	 do $system.Process.SetZEOF(1)
	 for i = 1:1:^chr
	 {
		 set filename = dirname_"/hs_ref_GRCh37.p5_"_^chr(i)_".fa"
		 write filename, !
		 open filename
		 use filename
		 
		 read rec //skip header
		 
		 set j = 1
		 for
		 {
			 read rec
			 quit:$zeof=-1
			 for k=1:1:$Length(rec)
			 {
				 set ch = $Extract(rec,k)
				 set ^humanDNA(^chr(i),j) = ch
				 set j = j + 1
			 }
		 }
		 set ^humanDNA(^chr(i)) = j - 1
		 close filename
	 }
	 do ENABLE^%SYS.NOJRN
}

/// This function creates global with the positions on th DNA
/// of the TCA or TCT (W is T or A) motifs
ClassMethod CreateTCWGlobal()
{
	 kill ^TCW
	 for i = 1:1:24
	 {
		 set k = 1
		 
		 set chpp = ^humanDNA(^chr(i),1)
		 set chp = ^humanDNA(^chr(i),2) 
		 
		 for j = 3:1:^humanDNA(^chr(i))
		 {
			 set ch = ^humanDNA(^chr(i),j)
			 set word = chpp_chp_ch
			 if ((word = "TCA") || (word = "TCT") || (word = "TGA") || (word = "AGA"))
			 {
				 if ((word = "TCA") || (word = "TCT"))
				 {
				 	set ^TCW(^chr(i),(j - 1)) = 1
				 }
				 elseif ((word = "TGA") || (word = "AGA"))
				 {
					 set ^TCW(^chr(i),(j - 1)) = -1
				 }
				 set k = k + 1
			 } 
			 set chpp = chp
			 set chp = ch
		 }
		 set ^TCW(^chr(i)) = (k - 1)
		 write i, !
	 }
}

/// Cancer mutations from 
/// Fredriksson et al (PMID: 25383969).
ClassMethod LoadMutations(dirname As %String = "/home/mkazanov/mutations")
{
	kill ^mutation
	kill ^sample
	do DISABLE^%SYS.NOJRN
	do $system.Process.SetZEOF(1)
	for i = 1:1:^cancerType
	{
		set filename = dirname_"/"_^cancerType(i)_"_short.txt"
		
		write filename, !
		open filename
		use filename
		 
		read rec //skip header

        set mutcnt = 0
        set ^sample(^cancerType(i)) = 0
		for
		{
		     read rec
			 quit:$zeof=-1
			 
			 set tsb = $piece(rec,$c(9),1)
 			 set chr = $piece(rec,$c(9),2)
 			 set pos = $number($piece(rec,$c(9),3))
 			 set ntfrom = $piece(rec,$c(9),4)
 			 set ntto = $piece(rec,$c(9),5)
 			 set isAPOBEC = $number($piece(rec,$c(9),6))
 			 set sampleEnrichment = $number($piece(rec,$c(9),7))
 			 set sampleMutLoad = $number($piece(rec,$c(9),8))
 			 
 			 if ($data(^sample(^cancerType(i),tsb)) = 0)
 			 {
	 			 set ^sample(^cancerType(i),tsb) = $lb(sampleEnrichment,sampleMutLoad)
 			     set ^sample(^cancerType(i)) = ^sample(^cancerType(i)) + 1
 			 }
 			 if ($data(^mutation(^cancerType(i),tsb)) = 0)
 			 {
	 			 set ^mutation(^cancerType(i),tsb) = 0
 			 }
 			 set ^mutation(^cancerType(i),tsb,chr,pos) = $lb(ntfrom,ntto,isAPOBEC)
			 set ^mutation(^cancerType(i),tsb) = ^mutation(^cancerType(i),tsb) + 1
		     set mutcnt = mutcnt + 1
		}
		set ^mutation(^cancerType(i)) = mutcnt
		close filename
	}
	do ENABLE^%SYS.NOJRN
}

/// Replication timing data for IMR90, MCF7, NHEK cell lines from the ENCODE project
/// https://www.encodeproject.org/search/?type=Experiment&assay_title=Repli-seq&assay_slims=Replication+timing
ClassMethod LoadReplicationTiming(dirname As %String = "/home/mkazanov/ReplicationTiming")
{
	set filenames("IMR90") = "wgEncodeUwRepliSeqImr90WaveSignalRep1.mybed"
	set filenames("MCF7") = "wgEncodeUwRepliSeqMcf7WaveSignalRep1.mybed"
	set filenames("NHEK") = "wgEncodeUwRepliSeqNhekWaveSignalRep1.mybed"
	
	do DISABLE^%SYS.NOJRN
	do $system.Process.SetZEOF(1)
	
	set bin = 0
	
	for i = 1:1:^cellLinesRT
	{
		set filename = dirname_"/"_filenames(^cellLinesRT(i))
		
		write filename, !
		open filename
		use filename
		
		for
		{
		     read rec
			 quit:$zeof=-1
			 set chr = $piece(rec,$c(9),1)
			 set startpos = $number($piece(rec,$c(9),2))
			 set endpos = $number($piece(rec,$c(9),3))
			 set rtval = $number($piece(rec,$c(9),4))
			 set ^RT(^cellLinesRT(i),chr,startpos) = $lb(endpos,rtval,bin)
		}
		close filename
	}
	
	do ENABLE^%SYS.NOJRN
}

/// Gene expression data is from 
/// the Broad GDAC Firehose
/// https://gdac.broadinstitute.org/
ClassMethod LoadExpression(dirname As %String = "/home/mkazanov/Expression")
{
	kill ^exp
	do DISABLE^%SYS.NOJRN
	do $system.Process.SetZEOF(1)
	
	for i = 1:1:^cancerType
	{
		set filename = dirname_"/"_"unpivot_expression_"_^cancerType(i)_".txt"
		
		write filename, !
		open filename
		use filename
		 
		for
		{
			read rec
			quit:$zeof=-1
			
			set genenameid = $piece(rec,$c(9),1)
			set geneid = $piece(genenameid,"|",2)
			set genename = $piece(genenameid,"|",1)
			set sample = $extract($piece(rec,$c(9),2),1,16)
			set expval = $number($piece(rec,$c(9),3))
			
			set ^exp(^cancerType(i),sample,geneid) = expval
		}
		close filename
	}
	
	do ENABLE^%SYS.NOJRN
}

/// Gene positions info is
/// from GENCODE 
/// https://www.gencodegenes.org/releases/current.html
ClassMethod LoadGenes(dirname As %String = "/home/mkazanov/Genes")
{
	kill ^gene
	do DISABLE^%SYS.NOJRN
	do $system.Process.SetZEOF(1)

    set filename = dirname_"/"_"export_gene4cache.txt"
		
	write filename, !
	open filename
	use filename

	for
	{
		read rec
		quit:$zeof=-1

        set chr = $piece(rec,$c(9),1)
        set startpos = $number($piece(rec,$c(9),2))
        set endpos = $number($piece(rec,$c(9),3))
        set geneid = $number($piece(rec,$c(9),6))
        set genename = $piece(rec,$c(9),7)
        set strand = $piece(rec,$c(9),4)
        set info = $piece(rec,$c(9),5)
        
        set ^gene(chr,startpos,endpos) = $lb(geneid,genename,strand,info)
	}
	close filename
	
	do ENABLE^%SYS.NOJRN
}

/// Gene positions info is
/// from GENCODE 
/// https://www.gencodegenes.org/releases/current.html
ClassMethod LoadGenesNew(dirname As %String = "/home/mkazanov/Genes")
{
	do DISABLE^%SYS.NOJRN
    
    Set stream=##class(%Stream.FileCharacter).%New()
	set filename = dirname_"/"_"export_gene4cache.txt"
	Set sc=stream.LinkToFile(filename)	
		
	write filename, !
		 
	Set rec=stream.ReadLine()
	
	While 'stream.AtEnd 
	{
			
		 set gene = ##class(ApobecExp.Gene).%New()	 
		 set gene.chr = $piece(rec,$c(9),1)
		 set gene.startpos = $number($piece(rec,$c(9),2))
		 set gene.endpos = $number($piece(rec,$c(9),3))
		 set gene.strand = $piece(rec,$c(9),4)
		 set gene.info = $piece(rec,$c(9),5)
		 set gene.GeneId = $number($piece(rec,$c(9),6))
		 set gene.GeneName = $piece(rec,$c(9),7)
		 set r = gene.%Save()
		 
		 if $$$ISERR(r)
		 {
			 write $System.Status.GetErrorText(r),!
		 }

		Set rec=stream.ReadLine()
 			 
	}
}

}
