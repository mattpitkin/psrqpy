"""
The classes defined here are hold information on an individual pulsar
or an interable list of pulsars.
"""


from .config import PSR_ALL_PARS, PSR_ALL


class Pulsar(object):
    """
    An object to hold a single pulsar. The class requires a pulsar name. The
    pulsar parameters are set as attributes of the class.

    Args:
        psrname (str): a string containing a pulsar name
        query (:class:`psrqpy.QueryATNF`): a query

    Additional keyword arguments are any of the valid queriable pulsar
    parameters.

    """

    def __init__(self, psrname, query=None, **kwargs):
        self._name = psrname
        self._query = query

        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        """
        Define the method return by repr
        """

        return self.name

    def keys(self):
        """
        Return a list of the class attribute names for allowed pulsar
        parameters.
        """

        keys = (
            PSR_ALL_PARS
            + [par + "_ERR" for par in PSR_ALL_PARS]
            + [par + "_REF" for par in PSR_ALL_PARS]
        )
        return [key for key in self.__dict__ if key in keys]

    def items(self):
        """
        Return a list of the class attribute values.
        """

        keys = (
            PSR_ALL_PARS
            + [par + "_ERR" for par in PSR_ALL_PARS]
            + [par + "_REF" for par in PSR_ALL_PARS]
        )
        return [value for key, value in self.__dict__.items() if key in keys]

    @property
    def name(self):
        """
        Return the pulsar name
        """

        return self._name

    def __getitem__(self, key):
        """
        If the class has a attribute given by the key then return it, otherwise
        generate a query for that key to set it.

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
            if ukey.endswith("_ERR") or ukey.endswith(
                "_REF"
            ):  # an error or reference parameter
                tkey = ukey[:-4]  # parameter name without error
            else:
                tkey = ukey

            if tkey not in PSR_ALL_PARS:
                raise KeyError(
                    '"{}" is not a recognised pulsar ' "parameter".format(tkey)
                )
            else:
                # generate a query for the key and add it
                if self._query is None:
                    try:
                        from .search import QueryATNF

                        self._query = QueryATNF()
                    except IOError:
                        raise Exception("Problem querying ATNF catalogue")

            psrrow = self._query.get_pulsar(pulsarname)

            if psrrow is None:
                raise Exception('Pulsar "{}" is unknown'.format(pulsarname))

            param = psrrow[ukey][0]  # required output parameter
            setattr(self, ukey, param)  # set output parameter value

            # set parameter value if an error value was requested
            if PSR_ALL[tkey]["err"]:
                if tkey != ukey and ukey.endswith(
                    "_ERR"
                ):  # asking for error, so set actual value
                    setattr(self, tkey, psrrow[tkey][0])
                else:  # asking for value, so set error
                    setattr(self, tkey + "_ERR", psrrow[tkey + "_ERR"][0])

            # set parameter value if a reference values was requested
            if PSR_ALL[tkey]["ref"]:
                if tkey != ukey and ukey.endswith(
                    "_REF"
                ):  # asking for reference, so set actual value
                    setattr(self, tkey, psrrow[tkey][0])
                else:  # asking for value, so set reference
                    setattr(self, tkey + "_REF", psrrow[tkey + "_REF"][0])

        return param

    def __getattr__(self, key):
        """
        If the class has a attribute given by the key then return it, otherwise
        generate a query for that key to set it (use the already defined
        __getitem__)
        """

        ukey = key.upper()

        # swapped from using hasattr to try...except... (see
        # https://hynek.me/articles/hasattr/)
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

    def __getstate__(self):
        """
        Define to allow pickling.
        """

        return self.__dict__

    def __setstate__(self, d):
        """
        Define to allow pickling.
        """

        self.__dict__.update(d)

    def __dir__(self):
        """
        Set this to what is returned for ipython's autocomplete (otherwise the
        custom ``__getattr__`` caused problems!)
        """

        return self.keys()

    def __str__(self):
        """
        Define the string method.
        """

        return self.name

    def __eq__(self, other):
        """
        Define '==' rich comparison methods. True if pulsars have the same
        name.
        """

        if not isinstance(other, Pulsar):
            return False

        return bool(self.name == other.name)

    def __ne__(self, other):
        """
        Define '!=' rich comparison methods. False if pulsars have the same name.
        """

        if not isinstance(other, Pulsar):
            return True

        return bool(self.name != other.name)

    def __copy__(self):
        """
        Define how the object should be copied with copy
        """

        attrs = {}
        for key, value in zip(self.keys(), self.items()):
            attrs[key] = value
        newpsr = type(self)(self.name, **attrs)

        return newpsr


class Pulsars(object):
    """
    Class to contain multiple :class:`~psrqpy.pulsar.Pulsar` objects.
    """

    def __init__(self):
        self._num_pulsars = 0  # number of pulsars in the object
        self._psrs = {}  # dictionary of Pulsar objects in the object, keyed to the name

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

        return None

    def __getstate__(self):
        """
        Define to allow pickling.
        """

        return self.__dict__

    def __setstate__(self, d):
        """
        Define to allow pickling.
        """

        self.__dict__.update(d)

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

        assert isinstance(psr, (Pulsar, Pulsars)), "psr is not a Pulsar type"

        if isinstance(psr, Pulsar):
            if psr.name not in self._psrs:
                self._num_pulsars += 1  # add one pulsar
                self._psrs[psr.name] = psr

        else:
            # check for duplicates
            for psrname in psr:
                if psrname not in self._psrs.keys():  # don't add duplicates
                    self._psrs[psrname] = psr[psrname]
                    self._num_pulsars += 1

    def remove_pulsar(self, psrname):
        """
        Remove a pulsar from the object. Only do one at a time.

        Args:
            psrname (str): a string with the name of a pulsar
        """

        assert isinstance(psrname, str), "psrname is not a string"

        if psrname in self._psrs:
            del self._psrs[psrname]
            self._num_pulsars -= 1

    def pop(self, psrname):
        """
        Remove a pulsar from the object and return the removed pulsar.

        Args:
            psrname (str): a string with the name of a pulsar
        """
        assert isinstance(psrname, str), "psrname is not a string"

        if psrname in self._psrs:
            self._num_pulsars -= 1
            return self._psrs.pop(psrname)

        return None

    def __str__(self):
        """
        Define string method
        """

        return "\n".join([self._psrs[psr].name for psr in self._psrs])
