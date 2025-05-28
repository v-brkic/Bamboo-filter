# Bioinformatics 1
University of Zagreb
Faculty of Electrical Engineering and Computing
www.fer.unizg.hr/predmet/bio1

# Bamboo Filter

```bash
# inside the project directory
make          # produces an executable named `bamboo_filter`
./bamboo_filter --genome=genomes/ecoli.fa \
  --kmer=20 --numKmers=10000 \          #starts the implementation with default values

./bamboo_filter --genome=genomes/ecoli.fa \
  --kmer=20 --numKmers=10000 \
  --capacity=4096 --bucketSize=4 \
  --loadFactor=0.8 --maxIter=2000           #starts the implementation with expected values
```
