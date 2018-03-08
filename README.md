This repository includes scripts to analyze and integrate Illumina calls on the HGSV Illumina trio samples. There are two major steps in the process:
#  Extract SVs from each individual algorithm and re-organize them in bed format
### Script: ./scripts/step1.standardize_vcfs_to_bed.py
### usage:
```
step1.standardize_vcfs_to_bed.py

positional arguments:
  input_path    directory of input vcf
  output_file   name of output bed file
  reference     reference genome used for bam alignment
  caller_names  file containing names of algorithms to be standardized
  contig_names  file containing names of contigs to be standardized

optional arguments:
  -h, --help    show this help message and exit

```

#  Integrate SVs from each algorithm into a non-redundant SV discovery set
### Script: ./scripts/step2.integrate_ILL_svs.py 
### usage:
```
step2.integrate_ILL_svs.py

positional arguments:
  input_path  directory of input vcf
  reference   reference genome used for bam alignment
  blacklist   blacklist regions to be excluded in the integration

optional arguments:
  -h, --help  show this help message and exit

```

