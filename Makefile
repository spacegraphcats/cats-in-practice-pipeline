all: do-acido-expand

clean:
	-rm -r 15genome
	-rm -r acido
	-rm -r output

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

output/frontier.p25.sig: acido/long-contigs.sig 15genome/minhashes.db
	-mkdir output
	python -m search.frontier_search acido/long-contigs.sig 15genome 0.0 -o output/frontier.p00.sig --purgatory
	python -m search.frontier_search acido/long-contigs.sig 15genome 0.05 -o output/frontier.p05.sig --purgatory
	python -m search.frontier_search acido/long-contigs.sig 15genome 0.1 -o output/frontier.p10.sig --purgatory
	python -m search.frontier_search acido/long-contigs.sig 15genome 0.15 -o output/frontier.p15.sig --purgatory
	python -m search.frontier_search acido/long-contigs.sig 15genome 0.2 -o output/frontier.p20.sig --purgatory
	python -m search.frontier_search acido/long-contigs.sig 15genome 0.25 -o output/frontier.p25.sig --purgatory

do-acido-expand: 15genome/minhashes.db output/frontier.p25.sig acido/acido.fa.gz.sig
	sourmash search acido/acido.fa.gz.sig output/frontier*.sig  acido/long-contigs.sig --containment -n 100
	sourmash search acido/acido.fa.gz.sig output/frontier*.sig  acido/long-contigs.sig -n 100
	scripts/frontier-containment-similarity.py acido/acido.fa.gz.sig output 0 0.05 .1 .15 .2 .25 -o output/acido-expand-frontier.csv

acido-expand-frontier.csv: do-acido-expand

####

#
# shewanella.mappedreads.fa is a collection of reads from podar data
# that maps to the Shewanella OS228 genome via bwa aln.  "Real" data,
# with known answer.
#

# prepared reads -- this is here only for record keeping & never
# needs to be done again.
XXXshewanella.abundtrim.gz:
	trim-low-abund.py --normalize 12 -V -Z 10 -M 2e9 -C 3 -k 21 shewanella.mappedreads.fa -o shewanella.abundtrim.gz --gzip

# download the prepared reads (27 MB) from OSF
shewanella.abundtrim.gz:
	curl -L 'https://osf.io/7az9p/?action=download' > shewanella.abundtrim.gz

# build cDBG
shew-reads/cdbg.gxt: shewanella.abundtrim.gz
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 4e9 shewanella.abundtrim.gz -o shew-reads

# build catlas
shew-reads/catlas.csv: shew-reads/cdbg.gxt
	python -m spacegraphcats.catlas shew-reads 1

# build minhashes
shew-reads/minhashes.db: shew-reads/catlas.csv shew-reads/contigs.fa.gz
	python -m search.make_catlas_minhashes -k 31 --scaled=1000 shew-reads

# download the shewanella genome from OSF
shew-reads/shewanella.fa.gz:
	mkdir -p shew-reads
	curl -L 'https://osf.io/fx4ew/?action=download' > shew-reads/shewanella.fa.gz

# compute shewanella genome signature
shew-reads/shewanella.fa.gz.sig: shew-reads/shewanella.fa.gz
	sourmash compute -k 31 --scaled=1000 shew-reads/shewanella.fa.gz -o shew-reads/shewanella.fa.gz.sig

# download the Shewanella baltica OS185 (strain variant) genome, too.
shew-reads/shewanella-OS185.fa.gz:
	mkdir -p shew-reads
	curl -L 'https://osf.io/9756q/?action=download' > shew-reads/shewanella-OS185.fa.gz

# compute OS185 genome signature
shew-reads/shewanella-OS185.fa.gz.sig:
	sourmash compute -k 31 --scaled=1000 shew-reads/shewanella-OS185.fa.gz -o shew-reads/shewanella-OS185.fa.gz.sig

# do various searches
output/shew-OS185-f10.sig: shew-reads/catlas.csv \
	shew-reads/shewanella-OS185.fa.gz.sig \
	shew-reads/shewanella.fa.gz.sig
	python -m search.frontier_search shew-reads/shewanella.fa.gz.sig shew-reads 0.0 --purgatory -o output/shew-f00.sig

	python -m search.frontier_search shew-reads/shewanella-OS185.fa.gz.sig shew-reads 0.0 --purgatory -o output/shew-OS185-f00.sig
	python -m search.frontier_search shew-reads/shewanella-OS185.fa.gz.sig shew-reads 0.1 --purgatory -o output/shew-OS185-f01.sig
	python -m search.frontier_search shew-reads/shewanella-OS185.fa.gz.sig shew-reads 0.5 --purgatory -o output/shew-OS185-f05.sig
	python -m search.frontier_search shew-reads/shewanella-OS185.fa.gz.sig shew-reads 1.0 --purgatory -o output/shew-OS185-f10.sig

shew-strain-search: output/shew-OS185-f10.sig
	sourmash search shew-reads/shewanella.fa.gz.sig output/shew-f00.sig output/shew-OS185-f*.sig -n 100 --containment
	sourmash search shew-reads/shewanella.fa.gz.sig output/shew-f00.sig output/shew-OS185-f*.sig -n 100
