"""
The classes defined here are hold information on an individual pulsar
or an interable list of pulsars.
"""

from __future__ import print_function, division

import warnings

from six import string_types, iteritems

from .config import PSR_ALL_PARS, PSR_ALL
from .utils import get_version

class Pulsar(object):
    """
    An object to hold a single pulsar. The class requires a pulsar name is required.

    Args:
        psrname (str): a string containing a pulsar name
        version (str): a string with the ATNF version to use for queries

    Additional keyword arguments are any of the valid queriable pulsar parameters.

    """

    def __init__(self, psrname, version=None, **kwargs):
        self._name = psrname
        self._version = version if not version else get_version()
        self._ephemeris = None

        for key, value in iteritems(kwargs):
            setattr(self, key, value)

    def __repr__(self):
        """
        Define the method return by repr
        """

        if self.have_ephemeris():
            return self.get_ephemeris()
        else:
            return self.name

    def keys(self):
        """
        Return a list of the class attribute names for allowed pulsar parameters
        """
        return [key for key in self.__dict__ if key in PSR_ALL_PARS+[par+'_ERR' for par in PSR_ALL_PARS]]

    def items(self):
        """
        Return a list of the class attribute values
        """
        return [value for key, value in iteritems(self.__dict__) if key in PSR_ALL_PARS+[par+'_ERR' for par in PSR_ALL_PARS]]

    @property
    def name(self):
        """
        Return the pulsar name
        """

        return self._name

    def __getitem__(self, key):
        """
        If the class has a attribute given by the key then return it, otherwise generate a
        query for that key to set it.
        
        Args:
            key (str): an item to get
        """

        ukey = key.upper()
        pulsarname = self.name

        param = None
        if key in self.__dict__:
            return self.__dict__[key]
        elif ukey in self.__dict__:
            return self.__dict__[ukey]
        else:
            if ukey[-4:] == '_ERR': # an error parameter
                tkey = ukey[:-4] # parameter name without error
            else:
                tkey = ukey

            if tkey not in PSR_ALL_PARS:
                raise KeyError('"{}" is not a recognised pulsar parameter'.format(tkey))
            else:
                # generate a query for the key and add it
                try:
                    from .search import QueryATNF
                    q = QueryATNF(params=tkey, psrs=pulsarname, version=self._version, include_errs=True)
                except IOError:
                    raise Exception('Problem querying ATNF catalogue')

            if q.num_pulsars != 1:
                raise Exception('Problem getting parameter "{}"'.format(tkey))

            param = q.get_dict()[ukey][0] # required output parameter
            setattr(self, ukey, param)    # set output parameter value

            # set parameter value if an error value was requested
            if PSR_ALL[tkey]['err']:
                if tkey != ukey: # asking for error, so set actual value
                    setattr(self, tkey, q.get_dict()[tkey][0]) # set parameter value
                else: # asking for value, so set error
                    setattr(self, tkey+'_ERR', q.get_dict()[tkey+'_ERR'][0]) # set error value

        return param

    def __getattr__(self, key):
        """
        If the class has a attribute given by the key then return it, otherwise generate a
        query for that key to set it (use the already defined __getitem__)
        """

        ukey = key.upper()

        # swapped from using hasattr to try...except... (see https://hynek.me/articles/hasattr/)
        try:
            return self.__dict__[key]
        except KeyError:
            try:
                return self.__dict__[ukey]
            except KeyError:
                try:
                    if ukey in PSR_ALL_PARS:
                        return self[ukey]
                except KeyError:
                    raise AttributeError(key)

    def __dir__(self):
        """
        Set this to what is returned for ipython's autocomplete (otherwise the custom
        ``__getattr__`` caused problems!)
        """

        return self.keys()

    def have_ephemeris(self):
        """
        Check whether we already have an ephemeris.
        
        Returns:
            bool: True if an ephemeris has already be downloaded.
        """

        if self._ephemeris:
            return True
        else:
            return False

    def get_ephemeris(self):
        """
        Query the ATNF to get the ephemeris for the given pulsar
        """

        ephem = self._ephemeris

        if not self.have_ephemeris():
            pulsarname = self.name

            try:
                from .search import QueryATNF
                q = QueryATNF(psrs=pulsarname, version=self._version, include_errs=True, get_ephemeris=True)
            except IOError:
                raise Exception('Problem querying ATNF catalogue')

            # set any parameters that can be set from the returned ephemeris
            ephem = q.get_dict()[pulsarname]

            self.set_ephemeris(ephem)

        return ephem

    def set_ephemeris(self, ephem=None):
        """
        Set attributes from the returned ephemeris

        Args:
            ephem (str): the ephemeris string
        """

        if not self._ephemeris and ephem:
            self._ephemeris = ephem # set ephemeris if it doesn't already exist

        assert isinstance(self._ephemeris, string_types), 'Ephemeris must be a string'

        # get ephemeris values
        ephemvals = [ev.split() for ev in ephem.split('\n') if len(ev.split()) > 1]

        for ev in ephemvals:
            if ev[0].upper() in PSR_ALL_PARS and not hasattr(self, ev[0].upper()):
                if PSR_ALL[ev[0].upper()]['format'][0] == 'S': # string type
                    setattr(self, ev[0].upper(), ev[1])
                elif PSR_ALL[ev[0].upper()]['format'][0] == 'i': # int type
                    try:
                        setattr(self, ev[0].upper(), int(ev[1]))
                    except ValueError:
                        warnings.warn('Could not set attribute for parameter "{}"'.format(ev[0].upper()), UserWarning)
                else: # float type
                    try:
                        setattr(self, ev[0].upper(), float(ev[1]))
                    except ValueError:
                        warnings.warn('Could not set attribute for parameter "{}"'.format(ev[0].upper()), UserWarning)

                # get errors if given
                if len(ev) == 3:
                    if PSR_ALL[ev[0].upper()]['err']:
                        try:
                            setattr(self, ev[0].upper()+'_ERR', float(ev[2]))
                        except ValueError:
                            pass

    def __str__(self):
        """
        Define the string method as a call to get_ephemeris and output the ephemeris
        """

        if self.have_ephemeris():
            return self.get_ephemeris()
        else:
            return self.name

    def __eq__(self, other):
        """
        Define '==' rich comparison methods. True if pulsars have the same name.
        """

        if not isinstance(other, Pulsar):
            return False
        else:
            if self.name == other.name:
                return True
            else:
                return False

    def __ne__(self, other):
        """
        Define '!=' rich comparison methods. False if pulsars have the same name.
        """

        assert isinstance(other, Pulsar), "You are not comparing two Pulsar types!"

        if self.name == other.name:
            return False
        else:
            return True

    def __copy__(self):
        """
        Define how the object should be copied with copy
        """

        attrs = {}
        for key, value in zip(self.keys(), self.items()):
            attrs[key] = value
        newpsr = type(self)(self.name, version=self._version, **attrs)
        newpsr.set_ephemeris(ephem=self._ephemeris)

        return newpsr


class Pulsars(object):
    """
    Class to contain multiple :class:`~psrqpy.pulsar.Pulsar` objects.
    """

    def __init__(self):
        self._num_pulsars = 0 # number of pulsars in the object
        self._psrs = {}       # dictionary of Pulsar objects in the object, keyed to the name
        self._got_ephemerides = False # set whether ephemerides have been got for all pulsars
        self._version = None

    def __iter__(self):
        """
        Iterator for the class
        """
        for psr in self._psrs:
            yield psr

    def __getitem__(self, key):
        """
        Define getitem to get a Pulsar object from the _psrs dictionary
        """

        if key in self._psrs.keys():
            return self._psrs[key]
        else:
            return None

    def __len__(self):
        """
        Define len method as the number of pulsars in the object
        """

        return self._num_pulsars

    def add_pulsar(self, psr):
        """
        Add a pulsar into the object.

        Args:
            psr (:class:`~psrqpy.pulsar.Pulsar`, :class:`~psrqpy.pulsar.Pulsars`): a
                :class:`~psrqpy.pulsar.Pulsar` object, or :class:`~psrqpy.pulsar.Pulsars`
                object
        """

        assert isinstance(psr, Pulsar) or isinstance(psr, Pulsars), 'psr is not a Pulsar type'

        if isinstance(psr, Pulsar):
            if psr.name not in self._psrs:
                self._num_pulsars += 1 # add one pulsar
                self._psrs[psr.name] = psr

                # check if the added pulsar already has an ephemeris
                if not psr.have_ephemeris():
                    self._got_ephemerides = False
        else:
            # check for duplicates
            for psrname in psrs:
                if psrname not in self._psrs.keys(): # don't add duplicates
                    self._psrs[psrname] = psrs[psrname]
                    self._num_pulsars += 1

                    # check whether any pulsars already have ephemerides
                    if not psrs[psrname].have_ephemeris() and self._got_ephemerides:
                        self._got_ephemerides = False

    def remove_pulsar(self, psrname):
        """
        Remove a pulsar from the object. Only do one at a time.

        Args:
            psrname (str): a string with the name of a pulsar
        """

        assert isinstance(psrname, string_types), 'psrname is not a string'

        if psrname in self._psrs:
            del self._psrs[psrname]

    def pop(self, psrname):
        """
        Remove a pulsar from the object and return the removed pulsar.

        Args:
            psrname (str): a string with the name of a pulsar
        """
        assert isinstance(psrname, string_types), 'psrname is not a string'

        if psrname in self._psrs:
            return self._psrs.pop(psrname)
        else:
            return None

    def have_ephemerides(self):
        """
        Check whether we have ephemerides for all pulsars
        """

        if self._got_ephemerides:
            return True
        else:
            return False

    def get_ephemerides(self, version=None):
        """
        Query the ATNF to get the ephemerides for all pulsars in the object.

        Args:
            version (str): the ATNF pulsar catalogue version to use.
        """

        if not self.have_ephemerides():
            self._version = version if not version else get_version() # get version of the ATNF catalogue to use

            psrnames = self._psrs.keys() # list of pulsar names

            try:
                from .search import QueryATNF
                q = QueryATNF(psrs=psrnames, version=self._version, include_errs=True, get_ephemeris=True)
            except IOError:
                raise Exception('Problem querying ATNF catalogue')

            for pulsarname in psrnames:
                if not self._psrs[pulsarname].have_ephemeris():
                    # set any parameters that can be set from the returned ephemeris
                    ephem = q.get_dict()[pulsarname]

                    self._psrs[pulsarname].set_ephemeris(ephem)

            self._got_ephemerides = True

        # return list of ephemerides
        return [self._psrs[psr].get_ephemeris() for psr in self._psrs]

    def __str__(self):
        """
        Define string method
        """

        if self.have_ephemerides():
            return '\n'.join(self.get_ephemerides())
        else:
            return '\n'.join([self._psrs[psr].name for psr in self._psrs])
