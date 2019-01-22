Class DnaTest.DataLoad [ Abstract ]
{

Parameter CHUNKSIZE = 40000;

ClassMethod LoadTest(dirname As %String = "/home/mkazanov/DNATEST")
{
		kill ^dnatest13

		do DISABLE^%SYS.NOJRN
		do $system.Process.SetZEOF(1)
		set filename = dirname_"/dna_random_positions_1M_3.txt"
		
		write filename, !
		open filename
		use filename

		set i = 1
		for
		{
			read rec
			quit:$zeof=-1

			set chr = $number($piece(rec,$c(9),1))
			set pos = $number($piece(rec,$c(9),2))
			set ^dnatest13(i) = $lb(chr,pos)
			set i = i + 1
		}
		set ^dnatest13 = i - 1
		
		close filename
		do ENABLE^%SYS.NOJRN
}

ClassMethod LoadHumanGenome(dirname As %String = "/home/mkazanov/HUMAN")
{
	 kill ^humanDNA
	 do DISABLE^%SYS.NOJRN
	 do $system.Process.SetZEOF(1)
	 for i = 1:1:^chr
	 {
		 set chr = ^chr(i)
		 set filename = dirname_"/hs_ref_GRCh37.p5_chr"_chr_".fa"
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
				 //set chunk = ((j - 1) \ 5) + 1
				 set chunk = ((j - 1) \ ..#CHUNKSIZE) + 1		
				 //set chunkpos = ((j - 1) # 5) + 1				 		 
				 set chunkpos = ((j - 1) # ..#CHUNKSIZE) + 1
				 if (chunkpos = 1){
				    set ^humanDNA(chr,chunk) = ch }
				 else {
					set ^humanDNA(chr,chunk) = ^humanDNA(chr,chunk)_ch }   
				 set j = j + 1
			 }
		 }
		 set ^humanDNA(chr) = j - 1
		 close filename
	 }
	 do ENABLE^%SYS.NOJRN
}

}
