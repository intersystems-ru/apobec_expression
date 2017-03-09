# apobec_expression

## Background

Somatic mutations is a commonly acceptable cause of human cancer. Known sources of mutations can have both exogenous and endogenous nature. Examples of the endogenous mutagens are ultraviolet light, ionizing radiation and tobacco smoke, whereas known endogenous mutagens are DNA repair deficiency, free radical species and deamination, which is mainly caused, as was recently discovered, by the activity of components of the human immune system - APOBEC (apolipoprotein B mRNA editing enzyme, catalytic polypeptide-like) family enzymes. Recent advances in sequencing of cancer genomes allows to track activity of particular mutagens based on so-called mutational signature, which describe a nucleotide context of mutations being unique for many mutational processes. This make possible to investigate the heterogeneity of the rate of mutations induced by various mutagen along the genome and recognize the possible like of mutation rate with the various epigenetic features such as transcription, replication timing and chromatin organization. The understanding of mutational heterogeneity along the genome is important for implication in statistical methods aiming the identification of the exhaustive catalog of cancer-driving genes. In this project we suggest exploration of the influence of epigenetic features on mutational heterogeneity along the genome for APOBEC enzymes and other endogenous mutagens based on new data from sequenced cancer genomes.

## Data acquisition

All executed methods are from DataLoad class.

1. Load human genome:
 LoadHumanGenome()
2. Load cancer mutations:
 LoadMutations()
3. Load replication timing data:
 LoadReplicationTiming()
4. Load gene expression data:
 LoadExpression()
5. Load gene information:
 LoadGenes()
 
Data dictionaries are from DataDictionary class.

1. Chromosome names dicrionary:
 SetChromosomeNames()
 SetChromosomeLowerName()
2. Cancer types dictionary:
 SetCancerTypes()
3. Cell lines dictionary:
 SetCellLinesReplication()
4. Cancer types - cell lines correspondence:
 SetCancer2RTCellLines()
5. Gene expression bins:
 SetExpressionBins()
 
## APOBEC mutations and replication timing

All methods are from class ReplicationTime

1. Create replication timing bins:
 CreateRTBins()
2. Set replication timing bins to replication timing intervals:
 SetRTBins()
3. Calculate number of TCW motifs in replication timing intervals:
 GetRTNumberOfTCW()
4. Calculate number of bases in replication timing intervals:
 GetRTNumberOfAll()
5. Calculate distribution of APOBEC mutations over replication timing:
 MutationRTAnalysis()
 
## APOBEC mutations and gene expression

All methods are from Expression class.

1. Set the maximum end of genes among intersecting genes
 SetGeneIntervals()
2. Divide DNA positions covering by genes into intervals covered by only single gene 
   and intervals covered by multiple genes:
 DivideGenesIntoIntervals()
3. Calculate number of TCW motifs in expression bins:
 GetExpTCW()
4. Calculate number of bases in expression bins:
 GetExpAll()
5. Calculate number of APOBEC and other mutations relative to expression bins: 
 GetExpMutation()
 
## APOBEC mutations and interaction of replication timing and gene expression

All methods are from ExpressionRT class.

1. Calculate number of TCW motifs in expression & replication timing bins:
 GetExpRTTCW()
2. Calculate number of bases in expression & replication timing bins: 
 GetExpRTAll()
3. Calculate number of APOBEC and other mutations relative to expression & replication timing bins:
 GetRTExpMutation()
 
## Export results

All methods are from ExportData class.

1. Prepare results of the APOBEC mutations density relative to replication timing for exporting:
 PlotRT()
2. Prepare results of the APOBEC mutations density relative to gene expression for exporting:
 PlotExp()
3. Prepare results of the APOBEC mutations density relative to replication timing & gene expression for exporting:
 PlotExpRT()
4. Export results of replication timing analysis:
 ExportPlotRT()
5. Export results of gene expression analysis 
 ExportPlotExp()
6. Export results of replication timing & gene epression analysis:
 ExportPlotExpRT()
 
 

 