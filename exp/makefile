all: score_distribuition.png calibration.png

score_distribuition.png:
	python figure_distribuition.py

example-random/Random_Signatures.tsv:
		python generate_random_picks.py 25

example-random/genesetDP_results.tsv: example-random/Random_Signatures.tsv
	python experiment.py 25

calibration.png: example-random/genesetDP_results.tsv example-random/Random_Signatures.tsv
	python figure_calibration.py
