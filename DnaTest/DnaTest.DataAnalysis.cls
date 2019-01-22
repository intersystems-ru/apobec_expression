Class DnaTest.DataAnalysis [ Abstract ]
{

Parameter CHUNKSIZE = 40000;

ClassMethod GetDNA(chr As %String, pos As %Integer) As %String
{
   //set chunk = ((pos - 1) \ 5) + 1
   set chunk = ((pos - 1) \ ..#CHUNKSIZE) + 1
   
   //set chunkPos = ((pos - 1) # 5) + 1
   set chunkPos = ((pos - 1) # ..#CHUNKSIZE) + 1
   
   set str = ""
	try {		     
		set str = $extract(^humanDNA(chr,chunk),chunkPos,chunkPos)      
	}
	catch e
	{
		write "chr=",chr,",startpos=",startpos,",endpos=",endpos,!
		write $zerror
		return
	}   	    	  
	   
   return str
}

ClassMethod DNATest()
{
	#dim sysDate As %TimeStamp = $zdatetime($ztimestamp,3,1,3)
	set i = $order(^dnatest13(""))
	while (i '= "")
	{
		set ch = $listget(^dnatest13(i),1)
		set chr = ^chr(ch)
		set pos = $listget(^dnatest13(i),2)
		set base = ..GetDNA(chr,pos)
		//set base = ##class(DNA.ViaBitMap).substring(chr,pos,pos)
		if (i # 1000 = 0) {
			write i, ! }
		set i = $order(^dnatest13(i))
	}
	#dim duration As %Integer = $System.SQL.DATEDIFF("ms",sysDate,$zdatetime($ztimestamp,3,1,3))
	write duration
}

}
