# variant_filtering_comparison
Code for comparing the variants found in a case using different filtering rules.

This script requires a TSV containing cases which were previously reported.

To run:
```
python get_variants_old_new_filtering.py \
--cases [TSV of reported cases and variants found] \
--dx_old [path to DNAnexus folder(s) of VCFs filtered with old filtering rules] \
--dx_new [path to DNAnexus folder(s) of VCFs filtered with new filtering rules] \
--old_path [full local path to folder where the VCFs filtered with old filtering rules will be downloaded]
--new_path [full local path to folder where the VCFs filtered with new filtering rules will be downloaded]
--out [name of output .xlsx file]
```