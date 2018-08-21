This is a collection of short scripts for fits operation, coordinate transfer, etc. They are mainly developed using Python 3.6, while some are written in bash. Download these scripts to your executable folder, such as bin, for using them in command line.

## Examples

### automatically save the ADS Bibtex entry for a paper to your .bib file

```bash
dbib.py 2018MNRAS.479.4041S
```
Note: 

replace 2018MNRAS.479.4041S with your own Bibcode

replace the file name in dbib.py to be your .bib file

The family name of the first author will be added to the front of the bibcode in your saved file


### open a fits file using DS9 with your own setting ("-zoom to fit" works well)
```bash
ds8 name1.fits name2.fits
```

### For the 53 functions defined in tools.py, below is an example to use cut_fits
```python
from tools import cut_fits
cut_fits(fits_in=a.fits, fits_ou=b.fits, x=0, y=0, size_x=1, size_y=1)
```

