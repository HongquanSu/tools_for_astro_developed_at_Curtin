#!/usr/bin/python

from os import system as s
import sys

def download_ads_bibcode(bibcode):
    link = 'http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=' + bibcode + '\&data_type=BIBTEX\&db_key=AST\&nocookieset=1'
    s('wget -O /Users/hongquansu/Downloads/bib.html ' + link)

    # get first author
    with open('/Users/hongquansu/Downloads/bib.html') as f0:
        for i, line in enumerate(f0):
            if i == 6:
                author_name = line.split('{{')[1].split('}')[0]
                break

    # write bibcode to bibtex file
    with open('/Users/hongquansu/Downloads/bib.html') as f1:
        with open('/Users/hongquansu/Dropbox/references/bibtex.bib', 'aw') as f2:
            for i, line in enumerate(f1):
                if i == 5:
                    line5_split = line.split('{')
                    line5_add_name = line5_split[0] + '{' + author_name + line5_split[1]
                    f2.write(line5_add_name)
                elif i > 5:
                    f2.write(line)


if __name__ == '__main__':
    # download the ADS bibcode and save it to my bibtex.bib file
    bibcode = sys.argv[1]
    download_ads_bibcode(bibcode)
