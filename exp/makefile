all: score_distribuition.png calibration.png

score_distribuition.png:
	python figure_calibration.py

example-random/Random_Signatures.tsv:
		python generate_random.py

example-random/genesetDP_results.tsv: example-random/Random_Signatures.tsv
	python experiment.py

calibration.png: example-random/genesetDP_results.tsv example-random/Random_Signatures.tsv
	python combine_plot.py
