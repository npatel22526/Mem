# Mem
### MEM : Mendelian Error Method to rapidly detect Deletion in whole exome and whole genome Trios 
#### Default Usage: perl MEM_windowAnalysis.pl -c caseSample.list -b backgroundSamplt.list -w 100000 -s 100000

##### Arguments:
	-a, --analysis			Analysis type [both or cases]
	-b, --bground			List of Control Files, ME files must be in bed format [Require if -a is both]
	-c, --cases			List of Cases Files [Require if run_from is not provided]
	-g, --genome			Genome file containing start and end of autosomal chromosome in b37 version [required, provided along with tool]
	-n, --name			Name to attched to intermediate File
	-r, --run_from			Round1 analysis for intermediate steps, available options are Round1, Merge1, Round2, Merge2
	-s, --slide 			Size of slide [Required][Default 100000 ]
	-w, --window			Size of window [Required][Default 2000000 ]
	-t, --tool			Full path to bedtools [required] [Default : bedtools]
	-M, --rescue_cases_ME 		Min ME count to Rescue the window  after Step1[Default 5 ]
	-X, --filter_bground_ME 	Max ME count in bground to filter window after Step1 [Default 2 ]	
	-Y, --filter_cases_ME 		Min ME count required cases to consider as cluster  [Default 1 ]
	-Z, --filter_bground_sample 	Max Sample count in bground to filter window after Step1 [Default 3 ]	
	-h, --help			Usage summary
	-v, --verbose			Detailed Usage Information
