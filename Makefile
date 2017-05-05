all: do-acido-expand

# build acido cDBG
acido/cdbg.gxt: data/acido.fa.gz
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 1e9 data/acido.fa.gz

### 15genome stuff

15genome-clean:
	-rm -r 15genome/

# build cDBG
15genome/cdbg.gxt:
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 4e9 data/15genome.fa.gz -o 15genome

# build catlas
15genome/catlas.csv: 15genome/cdbg.gxt
	python -m spacegraphcats.catlas 15genome 3

# build minhashes
15genome/minhashes.db: 15genome/catlas.csv
	python -m search.make_catlas_minhashes -k 31 --scaled=5000 15genome

# run search!
15genome-search: 15genome/minhashes.db
	python -m search.frontier_search data/15genome.5.fa.sig 15genome 0.1

acido/long-contigs.sig: acido/cdbg.gxt
	extract-long-sequences.py acido/contigs.txt -l 2000 > acido/long-contigs.fa
	sourmash compute -k 31 --scaled 5000 acido/long-contigs.fa -o acido/long-contigs.sig

acido/acido.fa.gz.sig: data/acido.fa.gz
	sourmash compute -k 31 --scaled 5000 data/acido.fa.gz -o acido/acido.fa.gz.sig

output/frontier.02.sig: acido/long-contigs.sig 15genome/minhashes.db
	-mkdir output
	python -m search.frontier_search acido/long-contigs.sig 15genome 0.0 -o output/frontier.00.sig
	python -m search.frontier_search acido/long-contigs.sig 15genome 0.1 -o output/frontier.01.sig
	python -m search.frontier_search acido/long-contigs.sig 15genome 0.2 -o output/frontier.02.sig

do-acido-expand: 15genome/minhashes.db output/frontier.02.sig acido/acido.fa.gz.sig
	sourmash search acido/acido.fa.gz.sig output/frontier*.sig  acido/long-contigs.sig --containment -n 100
	sourmash search acido/acido.fa.gz.sig output/frontier*.sig  acido/long-contigs.sig -n 100
