#!/usr/bin/python3
#vcc:E8DBD5521DAD0AF4939F3908EE541D1EFD5835B910FBD015BEF4ADD3D7F2D2C4
#sig:4447BE16A1EDE778F21AEF915AFEABA4C52C2EFA02E4E708ABED63CA69258110

'''
    NedSpider: Get information from NED, the NASA/IPAC Extragalactic Database
    ** In construction
    YJ Qin, Jan 12 / Tucson, AZ
'''

import os
import sys
import warnings

import re
import six

import requests
import urllib.parse as ps
from bs4 import BeautifulSoup

import numpy as np
from astropy.io import votable
from astropy.coordinates import SkyCoord
from astropy.utils.exceptions import AstropyWarning

# suppress astropy warnings
warnings.simplefilter('ignore', category=AstropyWarning)

# aux functions
def _split_link(bs_elm):
    text = bs_elm.text.strip().replace('*', '')
    try: href = bs_elm.find_all('a')[0]['href']
    except: href = None
    return (text, href)

def _s2i(text, bad_value=-1):
    try: val = int(text.strip())
    except: val = bad_value
    return val

def _s2f(text, bad_value=np.nan):
    try: val = float(text.strip())
    except: val = bad_value
    return val

# custom exception
class IdentificationError(Exception):
    strerror = "Failed to identify the source."
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class NedSpider:

    '''
        NedSpider: Get galaxy information from NED
        YJ Qin, Jan 12 2018 @ Tucson
    '''

    # NED mirror settings
    _ned_mirror_url = 'https://ned.ipac.caltech.edu/'
    _req_timeout = 30

    # NED cosmology settings
    _hconst, _omegam, _omegav, _corr_z = 73, 0.27, 0.73, 1

    _ned_galaxy_types = ['G', 'GClstr', 'GGroup', 'GPair', 'GTrpl', 'G_Lens']
    # _ned_galaxy_types = ['G']

    def __init__(self, name=None, radec=None, radius=1.,
            ned_mirror_url=None, req_timeout=None,
            match_only_galaxy=True, match_unique=False):

        # if there is a name, try it first.
        if isinstance(name, str): self._init_by_name(name=name, \
                ned_mirror_url=ned_mirror_url, req_timeout=req_timeout, \
                match_only_galaxy=match_only_galaxy, match_unique=match_unique)

        # if name query succeeded, done.
        if hasattr(self, '_candidates'):
            if len(self._candidates) > 0: return

        # otherwise, use RA/Dec
        if radec is not None: self._init_by_radec(radec=radec, radius=radius, \
                ned_mirror_url=ned_mirror_url, req_timeout=req_timeout, \
                match_only_galaxy=match_only_galaxy, match_unique=match_unique)

        # if RA/Dec query succeeded, done.
        if hasattr(self, '_candidates'):
            if len(self._candidates) > 0: return

        # cannot identify by either name or radec,
        raise IdentificationError('Failed to identify.')

    def _init_by_name(self, name, ned_mirror_url=None, req_timeout=None,
            match_only_galaxy=True, single_only=True, match_unique=False):

        '''
            Search a galaxy by its name.
            'name': str, name of the galaxy,
            'ned_mirror_url': mirror of NED service, default:
                    'https://ned.ipac.caltech.edu/'
            'req_timeout': timeout of each query, default: 30 sec.
            'match_unique': if multiple galaxies identified by this name,
                    take the first one (True) or raise an error (False).
        '''

        # check if name is valid
        if not isinstance(name, str):
            raise TypeError('Galaxy name is not a string.')
        if name.lower() in ['', 'none', 'nan', 'n/a', 'na', 'unknown']:
            raise RuntimeError('Invalid galaxy name.')
        if name.lower() in ['mw', 'lmc', 'smc', 'milky way']:
            raise RuntimeError('Galaxy is a local object.')
            # NED does have SMC & LMC

        # overwrite default mirror settings
        if ned_mirror_url is None:
            self.ned_mirror_url = NedSpider._ned_mirror_url
        else: self.ned_mirror_url = ned_mirror_url

        if req_timeout is None:
            self.req_timeout = NedSpider._req_timeout
        else: self.req_timeout = req_timeout

        # construct API
        req_url = self.ned_mirror_url + '/cgi-bin/objsearch?objname=' + \
                ps.quote(name) + '&of=table' + '&extend=NO' + '&img_stamp=NO'
        #print(req_url)
        # TODO: use urllib to construct url.
        # did not use VOTable because it does not include objid.
        req = requests.get(req_url, timeout=self.req_timeout)

        # convert this to votable // obsolete
        '''
        try:
            src_tab = votable.parse(six.BytesIO(req.content), pedantic=False)
            src_tab = src_tab.get_first_table().to_table().as_array()
            # any better solution? votable.Table is not Table or recarray
        except IndexError:
            raise RuntimeError('Failed to retrieve table')
        except: raise
        '''
        # no longer used. VOTable features not well supported.

        # parse the page
        soup = BeautifulSoup(req.text, "html5lib")
        src_tab = soup(text=re.compile(r'SOURCE LIST'))
        if len(src_tab) == 0: # object not identified
            raise IdentificationError('Source not identified.')
        # get the inner table (source list)
        src_tab = src_tab[0].parent.parent.find('table', recursive=False)
        src_rows = src_tab.find_all('tr')[2:] # skip header rows

        # get candidate sources
        self._candidates = []
        for row_i in src_rows: # for every identified object (row),
            cols_i = row_i.find_all('td') # break into table cells
            self._candidates.append(tuple( \
                    [_split_link(col_i) for col_i in cols_i]))
        # convert to list of tuple of (value, link)

        # take only galaxies?
        if match_only_galaxy: self._candidates = list(filter(\
                lambda x: x[4][0] in self._ned_galaxy_types, self._candidates))
        if single_only: self._candidates = list(filter(
                lambda x: x[4][0] == 'G', self._candidates))
        #for cad_i in self._candidates[:25]: print(cad_i[1][0]) # DEBUG


        # check how many objects matched
        N_src = len(self._candidates)
        if N_src < 1: raise IdentificationError( \
                'No source identified by this name.') # warum?
        if match_unique is True and N_src > 1:
            raise IdentificationError('Multiple sources identified.')

        # convert to tab.
        self._tabulate_candidates()

        return

    def _init_by_radec(self, radec, radius=1.,
            ned_mirror_url=None, req_timeout=None, match_only_galaxy=True,
            single_only=True, match_unique=False):

        '''
            Search a galaxy by its RA/Dec
        '''

        # convert RA/Dec to SkyCoord if necessary format.
        if not isinstance(radec, SkyCoord):

            ra, dec, crd = None, None, None # should be a SkyCoord in the end.
            radec_format_complaint = "'radec' should be a two-element " \
                    "tuple/list/array of RA and Dec, or a SkyCoord object."

            try: # try to construct like a two-element tuple/list/array
                # print(radec)
                ra, dec = radec # split into two
                if isinstance(ra, str) and isinstance(dec, str): # two strings
                    if (':' in ra) or ('h' in ra): # HMS/DMS format
                        crd = SkyCoord(ra, dec, unit=('hour', 'deg'))
                    else: # RA/Dec in decimal degrees, in string form
                        crd = SkyCoord(ra, dec, unit=('deg', 'deg'))
                elif isinstance(ra, float) and isinstance(dec, float): # floats
                    crd = SkyCoord(ra, dec, unit=('deg', 'deg'))
                else: raise TypeError(radec_format_complaint)

            # case: not iterable
            except TypeError: raise TypeError(radec_format_complaint)

            # case: iterable but wrong number
            except ValueError: # usually "too many to unpack"
                if isinstance(radec, str): # not a single string?
                    if (':' in radec) or ('h' in 'radec'):
                        crd = SkyCoord(radec, unit=('hour', 'deg'))
                    else: crd = SkyCoord(radec, unit=('deg', 'deg'))
                else: raise TypeError(radec_format_complaint)

            # complain anyway
            except: raise

            # sanity check
            if crd is None: raise RuntimeError('Failed to understand "radec".')
            else: radec = crd # overwrite radec.

        # overwrite default mirror settings
        if ned_mirror_url is None:
            self.ned_mirror_url = NedSpider._ned_mirror_url
        else: self.ned_mirror_url = ned_mirror_url

        if req_timeout is None:
            self.req_timeout = NedSpider._req_timeout
        else: self.req_timeout = req_timeout

        # convert radec in SkyCoord to hms:dms string.
        ra, dec = radec.to_string(style='hmsdms').split(' ')

        # construct query string
        req_url = self.ned_mirror_url + '/cgi-bin/objsearch?' \
                + 'lon=' + ps.quote(ra)  + '&' \
                + 'lat=' + ps.quote(dec) + '&' \
                + 'radius=' + ps.quote(str(radius)) + '&' \
                + 'in_objtypes1=Galaxies&' \
                + 'in_objtypes1=GPairs&' \
                + 'in_objtypes1=GTriples&' \
                + 'in_objtypes1=GGroups&' \
                + 'in_objtypes1=GClusters&' \
                + 'obj_sort=Distance+to+search+center&' \
                + 'of=table&' \
                + 'list_limit=0&' \
                + 'extend=NO&' \
                + 'img_stamp=NO&' \
                + 'search_type=Near+Position+Search'
        # NOTE: NED uses arcmin. Me too. TODO: use better solution

        # get and parse
        req = requests.get(req_url, timeout=self.req_timeout)
        soup = BeautifulSoup(req.text, "html5lib")
        src_tab = soup(text=re.compile(r'SOURCE LIST'))
        if len(src_tab) == 0: # object not identified
            raise IdentificationError('Source not identified.')
        # get the inner table (source list)
        src_tab = src_tab[0].parent.parent.find('table', recursive=False)
        src_rows = src_tab.find_all('tr')[2:] # skip header rows

        # get candidate sources
        self._candidates = []
        for row_i in src_rows: # for every identified object (row),
            cols_i = row_i.find_all('td') # break into table cells
            self._candidates.append(tuple( \
                    [_split_link(col_i) for col_i in cols_i]))
        # convert to list of tuple of (value, link)

        # take only galaxies?
        if match_only_galaxy: self._candidates = list(filter(\
                lambda x: x[4][0] in self._ned_galaxy_types, self._candidates))
        if single_only: self._candidates = list(filter(
                lambda x: x[4][0] == 'G', self._candidates))
        # for cad_i in candidates: print(cad_i)

        # check how many objects matched
        N_src = len(self._candidates)
        if N_src < 1: raise IdentificationError( \
                'No source identified by this name.') # warum?
        if match_unique is True and N_src > 1:
            #for i in self._candidates: print (i[1][0])
            #print('radius', radius)
            raise IdentificationError('Multiple sources identified.')

        # convert to tab.
        self._tabulate_candidates()

        return

    def _tabulate_candidates(self):

        '''
            Convert search results (list of tuples of tuples...) into recarray
            TODO: doc for myself.
        '''

        # Format of self._candidates
        # 0 row id & link, 1 name, 2 RA, 3 dec, 4 type, 5 RV, 6 z, 7 qual,
        # 8 mag, 9 ang sep, 10 refs, 11 notes, 12 phot, 13 posn, 14 vel & z,
        # 15 diam, 16 assoc, 17 imgs, 18 spec, 19 row & link
        self.candidates = np.empty(shape=(len(self._candidates),),
                dtype=[('name', 'S48'), ('RA_h', 'S24'), ('Dec_h', 'S24'),
                       ('RA_d', 'f8'), ('Dec_d', 'f8'), ('type', 'S8'),
                       ('vel', 'f8'), ('z', 'f8'), ('zqual', 'S8'),
                       ('mag', 'f8'), ('dist', 'f8'), ('N_ref', 'i4'),
                       ('N_note', 'i4'), ('N_phot', 'i4'), ('N_posn', 'i4'),
                       ('N_velz', 'i4'), ('N_diam', 'i4'), ('N_assoc', 'i4')])

        # convert and copy columns
        for i_cd, cd_i in enumerate(self._candidates):

            # convert hexadecimal HMS/DMS to degrees
            crd_i = SkyCoord(cd_i[2][0], cd_i[3][0], frame='icrs')

            # copy columns
            self.candidates[i_cd] = (
                cd_i[1][0], cd_i[2][0], cd_i[3][0], crd_i.ra.deg, crd_i.dec.deg,
                cd_i[4][0], _s2f(cd_i[5][0]), _s2f(cd_i[6][0]), cd_i[7][0],
                _s2f(cd_i[8][0]), _s2f(cd_i[9][0]), _s2i(cd_i[10][0]),
                _s2i(cd_i[11][0]), _s2i(cd_i[12][0]), _s2i(cd_i[13][0]),
                _s2i(cd_i[14][0]), _s2i(cd_i[15][0]), _s2i(cd_i[16][0])
            )

        # how many candidates in the list
        self.N_candidates = self.candidates.size

    def _retrieve_page(self, idx=None, name=None,):

        '''
            Retrieve NED detailed search result of a source.
            TODO: doc for myself.
        '''

        # create "cache" of soup objects if necessary
        if not hasattr(self, '_soup'): self._soup = {}

        # if the desired page is already retrieved, return soup directly
        if idx  in self._soup: return self._soup[idx]
        if name in self._soup: return self._soup[name]

        # when calling with neither name nor idx,
        if idx is None and name is None:
            if self.N_candidates > 1: # having multiple sources matched?
                raise RuntimeError('Multiple candidates matched in the list, ' \
                      'Need either idx or name to specify.')
            else: # otherwise, use the only source in the list.
                idx, name = 0, self.candidates[0]['name']

        # did we get source index and name in the previous step?
        if idx is None and name is None:

            if isinstance(idx, int): # if we have a valid index, use it.
                name = self.candidates[idx]['name']
                # do we need to validate the value/range of idx?

            else: # if idx is not explicitly given, match by name.
                id_msk = self.candidates['name'] == name.strip()
                N_fcad = id_msk.astype(int).sum() # number of matched
                if N_fcad > 1: raise IdentificationError(\
                        "More than one candidates matched by the name.")
                if N_fcad < 1: raise IdentificationError(\
                        "No candidate is matched by this name.")
                idx = idx_msk.tolist().index(True) # write idx with the matched
                name = name.strip() # overwrite name with stripped as we did.

        # url for the page of detailed information. parse and modify
        req_url = ps.urlsplit(self.ned_mirror_url + '/' + \
                self._candidates[idx][0][1]) # now in namedtuple
        query_par = dict(ps.parse_qsl(req_url.query)) # parse query str
        query_par['img_stamp'] = 'NO' # saves network traffic
        query_par['of'] = 'table' # return HTML table

        # apply our cosmological seetings
        query_par['hconst'] = str(self._hconst)
        query_par['omegam'] = str(self._omegam)
        query_par['omegav'] = str(self._omegav)
        query_par['corr_z'] = str(self._corr_z)

        # re-assemble url and send request
        req_url = req_url._replace(query=ps.urlencode(query_par)) # better way?
        req = requests.get(ps.urlunsplit(req_url), timeout=self.req_timeout)

        '''
        # else: # take the first object.
        # FIXME TODO TODO
        src_info = src_tab[0]
        self.name = src_info[1].decode('utf-8')
        self.ra, self.dec = src_info[2], src_info[3]
        self.obj_type = src_info[4].decode('utf-8')
        self.vel, self.z = src_info[5], src_info[6]

        # construct query with single source
        req_url = self.ned_mirror_url + \
                'cgi-bin/objsearch?objname=' + ps.quote(self.name) + \
                '&extend=NO&hconst=73&omegam=0.27&omegav=0.73' + \
                '&corr_z=1&of=table&zv_breaker=30000.0' + \
                '&list_limit=5&img_stamp=NO'
        # TODO: use urllib
        '''
        # obsolete. debug.

        # bs it and put into pot
        soup_i = BeautifulSoup(req.text, "html5lib")
        self._soup[idx]  = soup_i
        self._soup[name] = soup_i

        return soup_i

        '''
        # do we need to get its properties?
        if get_prop is not None:

            # convert to tuple, if get_prop is a string
            if isinstance(get_prop, str):
                if get_prop.lower() == 'all':
                    get_prop = ['alias', 'pos', 'diam', 'velz', 'phot', \
                          'imgs', 'spec', 'assoc']
                else: get_prop = (get_prop,)

            if not isinstance(get_prop, (list, tuple)):
                raise TypeError('"get_prop" should be a list, tuple or str.')

            # get properties of the source.
            if 'phot' in get_prop: self.phot = self.get_phot()
            if 'spec' in get_prop: self.spec = self.get_spec()

        return
        '''
        # obsolete, for the change of API. (see NOTE180114A)

    def alias(self, name=None, idx=None, galaxy_only=False, single_only=False,
              expand_aliases=False):

        '''
            Get cross-identified names of this source.
            TODO: doc
        '''

        # retrieve the detailed information page
        #print(self.candidates)
        soup_i = self._retrieve_page(name=name, idx=idx)

        # find the table of cross-identified names
        cid_elm = soup_i(text=re.compile(\
                r'CROSS-IDENTIFICATIONS'))[0].parent.parent
        cid_tab = cid_elm.find('table').find_all('td') # find table cells
        cid_names = [w.text.strip() for w in cid_tab] # get identified names
        cid_names = [re.sub(' +', ' ', w) for w in cid_names]
        # remove double space in source names. e.g. IRAS  03174-1935
        cid_names = list(zip(cid_names[0::2], cid_names[1::2]))
        # group name and type in pairs
        cid_names = list(filter(lambda x: x[0] not in \
                ['', 'Object Names', 'Type'], cid_names))
        # remove table headers and empty entries

        # remove other kind of sources (IrS, UvS, etc) if required
        if galaxy_only: cid_names = list(filter( \
                lambda w: w[1] in self._ned_galaxy_types, cid_names))
        if galaxy_only: cid_names = list(filter(lambda w: w[1] == 'G', cid_names))

        if expand_aliases: # expand NED-style names to other.
            expand_list = []
            for name_i, type_i in cid_names:
                name_spl_i = name_i.split(' ') # break into parts
                if name_spl_i[0] in ['NGC', 'PGC', 'UGC', 'IC', 'MESSIER']:
                    name_spl_i[1] = name_spl_i[1].lstrip('0') # remove leading 0
                    expand_list.append((' '.join(name_spl_i), type_i))
                    if name_spl_i[0] == 'MESSIER': # special case for M galaxies
                        expand_list.append(('M' + name_spl_i[1], type_i))
            cid_names = cid_names + expand_list

        # print('Object cross-identified as:')
        # for i in cid_names: print(i) # done.
        # debug

        # return cross-identified names in list of tuple of (name, type)
        return cid_names

    def photometry(self, name=None, idx=None, as_recarray=False):

        '''
            Get photometry of this object.
            TODO: doc
        '''

        N_phot, phot_link = None, None
        if idx is None and name is None:
            if self.N_candidates == 1:
                N_phot, phot_link = self._candidates[0][12]
            else: raise RuntimeError('Object not uniquely identified' + \
                    ' in the candidate list.')
            # in most cases, this branch should be triggered.

        if N_phot is None and phot_link is None: # failed in the previous step
            if isinstance(idx, int): # if we have a valid index, use it.
                N_phot, phot_link = self._candidates[idx][12]
            else: # try to indentify by name
                for cad_i in self._candidates:
                    if cad_i[1][0] == name: # exact match
                        if N_phot is None: N_phot, phot_link = cad_i[12]
                        else: raise IdentificationError(\
                                'More than one candidate matched by name')
                if N_phot is None:
                    raise IdentificationError( \
                            'No candidate is matched by this name.')

        '''
        # search photometry & SED link in the webpage
        phot_link = self._soup.find_all('a',
                text=re.compile(r'Photometric data point(s)'))
        print(phot_link)
        '''
        # obsolete: change of API. (see NOTE180111A)

        # no photometry measurement, return None
        if phot_link is None: return []

        # modify the link so that we can have a XML table
        # print(phot_link, "this is phot_link", N_phot, "this is N_phot")
        req_url = ps.urlsplit(self.ned_mirror_url + phot_link) # now in namedtuple
        query_par = dict(ps.parse_qsl(req_url.query)) # parse query str
        query_par['img_stamp'] = 'NO'
        query_par['of'] = 'xml_main'

        # apply our cosmological seetings
        query_par['hconst'] = str(self._hconst)
        query_par['omegam'] = str(self._omegam)
        query_par['omegav'] = str(self._omegav)
        query_par['corr_z'] = str(self._corr_z)

        # combine into url and send request
        req_url = req_url._replace(query=ps.urlencode(query_par))
        req = requests.get(ps.urlunsplit(req_url), timeout=self.req_timeout)

        # get the first table and convert to numpy recarray
        phot_table = votable.parse(six.BytesIO(req.content), \
                pedantic=False).get_first_table().to_table( \
                use_names_over_ids=True)
        phot_table.convert_bytestring_to_unicode(python3_only=False)
        # the project is partly in py2

        # if required, return numpy recarray
        if as_recarray: return phot_table.as_array()
        else: return phot_table # otherwise, return astropy table (default)

    def classification(self, name=None, idx=None):

        '''
            Get detailed classification of this source
            TODO: doc
        '''

        ned_objname = None # the 'official' object name, not the user-provided

        # if just single object,
        if idx is None and name is None:
            if self.N_candidates == 1:
                ned_objname = self._candidates[0][1][0]
            else: raise RuntimeError('Target not uniquely identified' \
                    ' in the candidate list.')

        # if we have a valid index, use it.
        if isinstance(idx, int):
            ned_objname = self._candidates[idx][1][0]
        else: # try to indentify by name
            for cad_i in self._candidates:
                if cad_i[1][0] == name: # exact match
                    if N_phot is None:
                        ned_objname = self._candidates[idx][1][0]
                    else: raise IdentificationError(\
                            'More than one candidate matched by name')
            if ned_objname is None:
                raise IdentificationError( \
                    'No candidate is matched by this name.')

        # construct query url
        req_url = self.ned_mirror_url + \
                '/cgi-bin/NEDatt?objname=' + ps.quote(ned_objname)
        req = requests.get(req_url)

        # get classification result table
        soup_i = BeautifulSoup(req.text, "html5lib")
        cid_elm = soup_i.find_all('table',
                attrs={'summary': 'Classification Results'})[0]
        cid_rows = cid_elm.find_all('tr')
        cid_results, sec_i = {}, ''
        for row_i in cid_rows:
            cells_i = row_i.find_all('th')
            if len(cells_i) == 1: # this is a section header.
                sec_i = cells_i[0].text
                cid_results[sec_i] = []
            else: # this is a row.
                cells_i = [i.text for i in row_i.find_all('td')][:-1]
                for i_cell, cell_i in enumerate(cells_i):
                    if cell_i == '...': cells_i[i_cell] = ''
                    else: cells_i[i_cell] = cell_i.replace(\
                            u'\xa0i\xa0\xa0', u'').replace('    ', ' ')
                if sec_i != '': cid_results[sec_i].append(tuple(cells_i))

        return cid_results

    def image(self, idx=None, name=None):

        # lower priority.
        raise NotImplementedError('In construction.')

        '''
            Get the list of archival images
        '''

        if idx is None and name is None:
            if self.N_candidates == 1:
                image_link = self._candidates[0][17][1]
            else: raise RuntimeError('Object not uniquely identified' + \
                    ' in the candidate list.')

        # if we have a valid index, use it.
        image_link = None, None
        if isinstance(idx, int):
            image_link = self._candidates[idx][17][1]
        else: # try to indentify by name
            for cad_i in self._candidates:
                if cad_i[1][0] == name: # exact match
                    if image_link is None: image_link = cad_i[17][1]
                    else: raise IdentificationError(\
                            'More than one candidate matched by name')
            if image_link is None:
                raise IdentificationError('No photometry information is ' \
                        'matched by this name.')

        # modify the link
        req_url = ps.urlsplit(self.ned_mirror_url + image_link) # in namedtuple
        query_par = dict(ps.parse_qsl(req_url.query)) # parse query str

        # apply our cosmological seetings
        query_par['hconst'] = str(self._hconst)
        query_par['omegam'] = str(self._omegam)
        query_par['omegav'] = str(self._omegav)
        query_par['corr_z'] = str(self._corr_z)

        # combine into url and send request
        req_url = req_url._replace(query=ps.urlencode(query_par))
        req = requests.get(ps.urlunsplit(req_url), timeout=self.req_timeout)

        # parse and analyze the table
        soup_i = BeautifulSoup(req.text, 'html5lib')

        # TODO: convert HTML table to some organized tabular format.
        # see: NOTE180111A

#a = NedSpider(name='NGC 1300')
#a = NedSpider(radec=('01h36m41.7s', '+15d47m01s'))
#print(a.alias(galaxy_only=True, expand_aliases=True))

#phot = a.photometry()
#print(phot)
#a.get_phot()
#NedSpider._get_phot(25832)

#print (NedSpider._get_class('')['Galaxy Morphology'])
