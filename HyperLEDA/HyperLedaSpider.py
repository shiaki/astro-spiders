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
import astropy.coordinates as coord
import astropy.units as u
import numpy as np

class HyperLedaSpider:

    hyperleda_mirror = 'http://leda.univ-lyon1.fr/'
    req_timeout = 30

    def __init__(self, name=None, radec=None, objtype='G'):

        # if name is available, try this first,
        if name is not None:

            # argument type check
            if not isinstance(name, str):
                raise TypeError("'name' is not a string.")

            # query with name
            data_dic = self._init_by_name(name)
            if data_dic != {}: self.data = data_dic

        # if failed with name, use RA/Dec,
        if not hasattr(self, "data") and (radec is not None):

            # argument type check
            if not isinstance(radec, (tuple, list)):
                raise TypeError("'radec' is not a 2-tuple.")

            data_dic = self._init_by_radec(radec, objtype)
            if data_dic != {}: self.data = data_dic

        # complain if failed.
        if not hasattr(self, 'data'):
            raise RuntimeError('')

    def _init_by_radec(self, radec, objtype, sort_by_dist):

        # construct IAU-style sky coord. # TODO: other formats
        crd = coord.SkyCoord(ra=radec[0], dec=radec[1], unit=('deg', 'deg'))
        crd_str = crd.to_string(style='hmsdms').encode('utf-8')
        crd_str = 'J' + crd.ra.to_string(precision=1, sep=''\
                ).encode('utf-8') + crd.dec.to_string(precision=0, sep='', \
                alwayssign=True).encode('utf-8')
        req_url = self.hyperleda_mirror + 'ledacat.cgi?o=' + ps.quote(crd_str)

        # query and parse,
        req = requests.get(req_url, timeout=self.req_timeout)
        soup = BeautifulSoup(req.text, "html5lib")

        # dict to hold the results
        query_info = []

        # get name tables
        name_tabs = soup.find_all('a', text=re.compile(r'search in the field'))
        for i_name_tab, name_tab_i in enumerate(name_tabs):
            name_tab_i = name_tab_i.parent.find_all('a')
            name_i = name_tab_i[0].text
            query_info.append({})
            query_info[i_name_tab] = {'name': name_i}

        # get header tables for each object
        header_tabs = soup.find_all('a', text=re.compile(r'Celestial position'))
        for i_header_tab, header_tab_i in enumerate(header_tabs):
            header_tab_i = header_tab_i.parent.parent.parent.findChildren('table')
            coord_rows = header_tab_i[0].find_all('tr')
            coords = [tuple([td_i.text for td_i in \
                    tr_j.find_all('td')]) for tr_j in coord_rows]
            if coords[-1][0] == '': # convert 'Precision' key
                coords[-1] = tuple([w.strip() for w in coords[-1][1].split(':')])
            coords = {key_i: val_i for key_i, val_i in coords} # convert to dict
            if len(header_tab_i) > 1: # if there is alternative names, make a list
                alt_names = [w.text for w in header_tab_i[1].find_all('td')]
            else: alt_names = [] # otherwise, give an empty list.
            query_info[i_header_tab]['coords'] = coords
            query_info[i_header_tab]['alt_names'] = alt_names

        # get property table
        data_tabs = soup.find_all('td', text=re.compile(r'Parameter'))
        for i_data_tab, data_tab_i in enumerate(data_tabs):
            data_rows = data_tab_i.parent.parent.parent.find_all('tr')
            data_dic_i = {}
            for row_i in data_rows[1:]: # skip the header
                row_i = [w.text.strip() for w in row_i.find_all('td')]
                par_i, unit_i, descr_i = row_i[0], row_i[2], row_i[3]
                # split value and error.
                if u'±' not in row_i[1]: val_i, err_i = row_i[1], 'NaN'
                else: val_i, err_i = [w.strip() for w in row_i[1].split(u'±')]
                try:  val_i, err_i = float(val_i), float(err_i)
                except: pass # convert numeric values to float type
                data_dic_i[par_i] = (val_i, err_i, unit_i, descr_i)
            query_info[i_data_tab]['data'] = data_dic_i

        # take only galaxies, if necessary
        if objtype is not None:
            query_info_tr = []
            for idx, val in query_info.items():
                if val['data']['objtype'][0] != objtype: continue
                query_info_tr.append(val)
            query_info = query_info_tr # overwrite

        # sort with separation, if necessary
        if sort_by_dist is True:
            dist, query_info_sorted = []
            for val in query_info:
                crd_i = coord.SkyCoord(val['coords']['J2000'], unit=('hour', 'deg'))
                dist.append(crd.separation(crd_i).arcmin)
            for sort_idx_i in np.argsort(dist):
                query_info_sorted.append(query_info[sort_idx_i])


    def _init_by_name(self, name):

        # query
        req_url = self.hyperleda_mirror + 'ledacat.cgi?o=' + ps.quote(name)
        req = requests.get(req_url, timeout=self.req_timeout)

        # parse
        soup = BeautifulSoup(req.text, "html5lib")

        # find coordinates and alternate names
        header_tab = soup.find_all('a', text=re.compile(r'Celestial position'))
        if len(header_tab) == 0: return {} # return empty dic if query failed.
        header_tab = header_tab[0].parent.parent.parent.findChildren('table')
        # first: celestial coord, second: alt names

        # read the table of coordinates
        coord_rows = header_tab[0].find_all('tr')
        coords = [tuple([td_i.text for td_i in \
                tr_j.find_all('td')]) for tr_j in coord_rows]
        if coords[-1][0] == '': # convert 'Precision' key
            coords[-1] = tuple([w.strip() for w in coords[-1][1].split(':')])
        coords = {key_i: val_i for key_i, val_i in coords} # convert to dict

        # make RA/Dec machine-readable # TODO

        # read the table of alternate names
        if len(header_tab) > 1:
            alt_names = [w.text for w in header_tab[1].find_all('td')]
        else: alt_names = []

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

        return data_dic

if __name__ == '__main__':
    a = HyperLedaSpider('M87')
    print(a.data)
