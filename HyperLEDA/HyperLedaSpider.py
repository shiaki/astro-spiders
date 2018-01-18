#!/usr/bin/python3
#vcc:3AFCDB39F9CB43216C7B5EA689CF252C10734020C56C68DE3EA11776533D85C5
#sig:9107C73E452A44CD4B99B13BEB0442BEDE0C64FC48B13EBC56A193DC142138BA

'''
    HyperLedaSpider: Get galaxy information from HyperLEDA database
'''

import os, sys, re, six, requests
import urllib.parse as ps
from bs4 import BeautifulSoup

from astropy.io import votable

class HyperLedaSpider:

    hyperleda_mirror = 'http://leda.univ-lyon1.fr/'
    req_timeout = 30

    def __init__(self, name):

        # construct api
        req_url = self.hyperleda_mirror + 'ledacat.cgi?o=' + ps.quote(name)
        req = requests.get(req_url, timeout=self.req_timeout)

        # cook
        soup = BeautifulSoup(req_text, "html5lib")

        # find coordinates and alternate names
        header_tab = soup.find_all('a', text=re.compile(r'Celestial position'))
        header_tab = header_tab[0].parent.parent.parent.findChildren('table')
        # first: celestial coord, second: alt names

        # read the table of coordinates
        coord_rows = header_tab[0].find_all('tr')
        coords = [tuple([td_i.text for td_i in \
                tr_j.find_all('td')]) for tr_j in coord_rows]
        if coords[-1][0] == '': # convert 'Precision' key
            coords[-1] = tuple([w.strip() for w in coords[-1][1].split(':')])
        coords = {key_i: val_i for key_i, val_i in coords} # convert to dict

        # make RA/Dec machine-readable # TODO, only J2000?

        # read the table of alternate names
        alt_names = [w.text for w in header_tab[1].find_all('td')]

        # find data table
        data_tab = soup.find_all('a', text=re.compile(r'objtype'))
        data_tab = data_tab[0].parent.parent.parent
        data_rows = data_tab.find_all('tr')
        data_dic = {}
        for row_i in data_rows[1:]: # skip the header
            row_i = [w.text.strip() for w in row_i.find_all('td')]
            par_i, unit_i, descr_i = row_i[0], row_i[2], row_i[3]
            # split value and error.
            if u'±' not in row_i[1]: val_i, err_i = row_i[1], 'NaN'
            else: val_i, err_i = [w.strip() for w in row_i[1].split(u'±')]
            try: val_i, err_i = float(val_i), float(err_i)
            except: pass # convert numeric values to float type
            data_dic[par_i] = (val_i, err_i, unit_i, descr_i)
        print (data_dic)

if __name__ == '__main__': HyperLedaSpider('M87')
