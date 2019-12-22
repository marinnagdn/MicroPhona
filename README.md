# MicroPhona

AUTHORS:
- Aimeric Dabin
- Marinna Gaudin
- Jolann Pommellec

University of Nantes, France
Master 2 Bioinformatics

### CONTEXT ###
The study of bacteria interactions in a given environment is a challenge. Co-occurrence networks are generally used to provide a graphic visualization of potential relationships between bacteria, but they are only based on correlations. However a correlation doesn't explain why the interaction might occur. For this reason, we developed a tool, that we called MicroPhona, that enables the revelations of potential bacteria conversations by inferring co-occurrence graphs with metabolic data.

### SOURCE METABOLIC DATA ###
- EMBL database :
D. Machado, S. Andrejev, M. Tramontano, et K. R. Patil, « Fast automated reconstruction of genome-scale metabolic models for microbial species and communities », Nucleic Acids Res., vol. 46,n o 15, p. 7542-7553, 06 2018.

- Recon3D database
« Recon3D enables a three-dimensional view of gene variation in human metabolism | Nature Biotechnology ». [En ligne]. Disponible sur: https://www.nature.com/articles/nbt.4072. 

### INPUTS ###

- .gml files generated with a tool such as FlashWeave (see below for more informations)

	*** We consider the nodes writing with the phylogenetic way
	Ex:	Root;k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium;s__animali

	/!\ KEEP ONLY THE BACTERIA THAT HAVE THEIR GENUS AND SPECIES. ONLY GENUS WILL BREAK THE PROGRAM.
	
- Parameters file :
	
		# 1 for YES, 0 for NO
	
		groups : prefix of your gml files
	
		suffix : suffix of your outputs files
	
		gml : new gml annotated files ?
	
		graph : graph png files ?
	
		oriented : oriented graphs ? 
	
		# oriented:1 -> ONLY IRREVERSIBLE REACTIONS
		# oriented:0 -> Resersible + Irreversible Reactions
	
		human : human as a node ?
		threshold : minimum threshold for absolute correlations
		<optional> bacteria : number of the bacteria you want to highlight. see index_bacteria.txt to find the corresponding of your bacteria of interest ]
		<optional> metabolites : KEGG ID of your metabolites of interest. Go on KEGG Databse to get them. ]

### OUTPUTS ###

- New gml annotated files with metabolic data
- Png graphs for a betetr vizualisation


### EXAMPLE OF USE ###

command line : 

		python3 -W ignore microphona.py parameters.txt

---> See toy_example directory to see the results


### SOURCE OF EXAMPLE ####

“Weight gain in anorexia nervosa does not ameliorate the faecal microbiota, branched chain fatty acid profiles, and gastrointestinal complaints” Isabelle Mack, Ulrich Cuntz, Claudia Grämer, Sabrina Niedermaier, Charlotte Pohl, Andreas Schwiertz, Kurt Zimmermann, Stephan Zipfel, Paul Enck & John Penders Scientific Reports volume 6, Article number: 26752 (2016)

Get DATA (MgniFy Database) :
https://www.ebi.ac.uk/metagenomics/studies/MGYS00001279
-> where we downloaded the tsv files

--> See tsv_files directory to see the origin abundance tables

### GET GML FROM TSV ###

*** FlashWeave : https://github.com/meringlab/FlashWeave.jl/blob/master/README.md

J. Tackmann, J. F. Matias Rodrigues, et C. von Mering, « Rapid Inference of Direct Interactions in Large-Scale Ecological Networks from Heterogeneous Microbial Sequencing Data », Cell Syst, vol. 9, n o 3, p. 286-296.e8, sept. 2019
	
	1) Install Julia
	2) Install FlashWeave
	3) On Julia :
				using FlashWeave
				data = "tsv_files/NW.tsv"
				data_network = learn_network(data, sensitive=true, heterogeneous=false, transposed = true)
				save_network("toy_examples/NW.gml", data_network)


