# ProNorm

Normalise 2 sets of peptide abundance data, with any number of biological replicates, by using an algorithm similar
to the Progenesis QI normalise-to-all-proteins option. The outputs are the individual pair-wise fold changes and the
average fold change between the two sets across replicates. This JavaScript version allows custom data filtering
where Progenesis does not, allowing it to work with samples that contain consistently large numbers of missing data
points between the 2 sets.

## Install

Requires Node.js to be installed, and then in the root directory of this repo:

```
npm install
```

## Usage

```
./normalise.js --in ./data.csv --set1=2,3,4,5 --set2=6,7,8,9 --gene-col=10 --group > output.csv
```

The `--in` option is required and is the raw data in CSV form with the first row containing column headings.

Use `--set1=` and `--set2=` to provide a list of column numbers containing set 1 and set 2 respectively. The first
column is column 1.

Specify the column containing the gene or protein name with `--gene-col`.

Finally it is useful to have the algorithm group the peptides into proteins for the output so use `--group`. Note that
it does not filter by those with only one peptide so you may need to exclude this option and filter manually.

## Algorithm

1. Calculate log fold change for all run combinations
2. Using the median and Median Absolute Deviation (MAD), iteratively filter the outliers in each
3. Find mean or median which should be 0 and therefore corresponds to a log-space shift to be removed
4. Convert shift to abundance-space scalar
5. Scale the appropriate run in a given pair of runs using the above scalar
6. Group peptides into proteins (optional)
7. Recalculate log fold changes using scaled runs

