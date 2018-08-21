This is a collection of short scripts for fits operation, coordinate transfer, etc. They are mainly developed using Python 3.6, while some are written in bash. For using them in command line, you need to download these scripts to your executable folder, such as bin. You may need to revise the scripts to make them work for your purpose.

## Examples

### For the 53 functions defined in tools.py, below is an example to use cut_fits
```python
from tools import cut_fits
cut_fits(fits_in=a.fits, fits_ou=b.fits, x=0, y=0, size_x=1, size_y=1)
```
Note:

This will make a cutout of a fits image using Cutout2D and save it to the given file

fits_in: input filename of fits image

fits_ou: output filename of fits image

x: ra or galactic longitude of the image center in degrees

y: dec or galactic latitude of the image center in degrees

size_x: the width of the image in degrees

size_y: the height of the image in degrees

### Open a fits file using DS9 with your own setting ("-zoom to fit" works well)                                   
```bash
ds8 name1.fits name2.fits
```

### Automatically save the ADS Bibtex entry for a paper to your .bib file
```bash
dbib.py 2018MNRAS.479.4041S
```
Note: 

replace 2018MNRAS.479.4041S with your own Bibcode

replace the file name in dbib.py to be your .bib file

The family name of the first author will be added to the front of the bibcode in your saved file

### Quickly find a paper in ADS using three numbers (year, volume, page)
```bash
ref 2018 479 4041
```
Note:

This will open the ADS webpage of that paper 
