Scripts associated with aptamer processing. Note that this program requires installation of the Vienna Package. Please see http://www.tbi.univie.ac.at/~ivo/RNA/ for installation and citation of that package.

If you use this software, please cite:

Thiel WH, Bair T, Wyatt Thiel K, Dassie JP, Rockey WM, Howell CA, Liu XY, Dupuy AJ, Huang L, Owczarzy R, Behlke MA, McNamara JO, Giangrande PH. Nucleotide bias observed with a short SELEX RNA aptamer library. Nucleic Acid Ther. 2011 Aug; 21 (4) :253-63. PubMed PMID:21793789; PubMed Central PMCID: PMC3198618.


The main script in this project is : predict_structures.py

An example of usage would be python predict_structures.py -v 2 Final_Rd12.fa the defaults are usually pretty good choices. You may want to invesigate the --prefix and --suffix options depending on how your experiment is structured, these default to empty but you may want to put in sequences if you have constant regions before and after the variable regions of your sequences. You can also choose to use mfold instead of Vienna for structural prediction, both return similar but not identical results, both need to be installed and executable by this script to be used.


The output of predict_structures can be used by create_graph.py

Running with just the default options should produce good inital results, try using the seed option to make a graph that is a little less hairball in appearance.
