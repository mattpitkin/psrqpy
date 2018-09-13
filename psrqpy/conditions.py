"""
Parser for logical conditions.
"""

from __future__ import print_function, division

import re
import numpy as np
from .config import *

from six import string_types

# set formatting of warnings to not include line number and code (see
# e.g. https://pymotw.com/3/warnings/#formatting)
def warning_format(message, category, filename, lineno, file=None, line=None):
    return '{}: {}'.format(category.__name__, message)

warnings.formatwarning = warning_format

# string of logical expressions for use in regex parser
LOGEXPRS = (r'(\bAND\b'        # logical AND
            r'|\band\b'        # logical AND
            r'|\&\&'           # logical AND
            r'|\bOR\b'         # logical OR
            r'|\bor\b'         # logical OR
            r'|\|\|'           # logical OR
            r'|!='             # not equal to
            r'|=='             # equal to
            r'|<='             # less than or equal to
            r'|>='             # greater than or equal to
            r'|<'              # less than
            r'|>'              # greater than
            r'|\('             # left opening bracket
            r'|\)'             # right closing bracket
            r'|\bNOT\b'        # logical NOT
            r'|\bnot\b'        # logical NOT
            r'|!'              # logical NOT
            r'|~'              # logical NOT
            r'|\bASSOC\b'      # pulsar association
            r'|\bassoc\b'      # pulsar association
            r'|\bTYPE\b'       # pulsar type
            r'|\btype\b)'      # pulsar type
            r'|\bBINCOMP\b'    # pulsar binary companion type
            r'|\bbincomp\b'    # pulsar binary companion type
            r'|\bEXIST\b'      # pulsar parameter exists in the catalogue
            r'|\bexist\b'      # pulsar parameter exists in the catalogue
            r'|\bERROR\b'      # condition on parameter error
            r'|\berror\b)'     # condition on parameter error


def condition(table, expression, exactMatch=False):
    """
    Apply a logical expression to a table of values.

    Args:
        table (:class:`astropy.table.Table`, :class:`pandas.DataFrame`): a
            table of pulsar data
        expression (str, :class:`~numpy.ndarray): a string containing a set of
            logical conditions with respect to pulsar parameter names (also
            containing `conditions
            <http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal#condition>`_
            allowed when accessing the ATNF Pulsar Catalogue), or a boolean
            array of the same length as the table.
        exactMatch

    Returns:
        table: the table of values conforming to the input condition.
            Depending on the type of input table the returned table will either
            be a :class:`astropy.table.Table` or :class:`pandas.DataFrame`.

    Example:
        Some examples this might be:

        1. finding all pulsars with frequencies greater than 100 Hz

        >>> newtable = condition(psrtable, 'F0 > 100')

        2. finding all pulsars with frequencies greater than 50 Hz and
        period derivatives less than 1e-15 s/s.

        >>> newtable = condition(psrtable, '(F0 > 50) & (P1 < 1e-15)')

        3. finding all pulsars in binary systems

        >>> newtable = condition(psrtable, 'TYPE(BINARY)')

        4. parsing a boolean array equivalent to the first example

        >>> newtable = condition(psrtable, psrtable['F0'] > 100)

    """

    from astropy.table import Table
    from pandas import DataFrame

    # check if expression is just a boolean array
    if isinstance(expression, np.ndarray):
        if expression.dtype != np.bool:
            raise TypeError("Numpy array must be a boolean array")
        elif len(expression) != len(table):
            raise Exception("Boolean array and table must be the same length")
        else:
            return table[expression]
    else:
        if not isinstance(expression, string_types):
            raise TypeError("Expression must be a boolean array or a string")

    # parse the expression string and split into tokens
    reg = re.compile(LOGEXPRS)
    tokens = reg.split(expression)
    tokens = [t.strip() for t in tokens if t.strip() != '']

    if isinstance(table, Table):
        # convert astropy table to pandas DataFrame
        tab = table.to_pandas()
    elif not isinstance(table, DataFrame)
        raise TypeError("Table must be a pandas DataFrame or astropy Table")
    else:
        tab = table

    matchTypes = ['ASSOC', 'TYPE', 'BINCOMP', 'EXIST', 'ERROR']

    # parse through tokens and replace as required
    ntokens = len(tokens)
    newtokens = []
    while i < ntokens:
        if tokens[i] in ['&&', 'AND']:
            # replace synonyms for '&' or 'and'
            newtokens.append('&')
        elif tokens[i] in ['||', 'OR']:
            # replace synonyms for '|' or 'or'
            newtokens.append('|')
        elif tokens[i] in ['!', 'NOT', 'not']:
            # replace synonyms for '~'
            newtokens.append('~')
        elif tokens[i].upper() in matchTypes:
            if ntokens < i+3:
                warnings.warn("A '{}' must be followed by a '(NAME)': ignoring in query".format(tokens[i].upper()), UserWarning)
            elif tokens[i+1] != '(' or tokens[i+3] != ')'
                warnings.warn("A '{}' must be followed by a '(NAME)': ignoring in query".format(tokens[i].upper()), UserWarning)
            else:
                if  tokens[i].upper() == 'ASSOC':
                    if 'ASSOC' not in tab.keys():
                        warnings.warn("'ASSOC' parameter not in table: ignoring in query", UserWarning)
                    elif exactMatch:
                        newtokens.append('(ASSOC == "{}")'.format(tokens[i+2]))
                    else:
                        assoc = np.array([tokens[i+2] in a for a in table['ASSOC']])
                        newtokens.append('(@assoc)')
                elif tokens[i].upper() == 'TYPE':
                    if tokens[i+2].upper() == 'BINARY':
                        if 'BINARY' not in tab.keys():
                            warnings.warn("'BINARY' parameter not in table: ignoring in query", UserWarning)
                        else:
                            newtokens.append('(BINARY != "None")')
                    else:
                        if 'TYPE' not in tab.keys():
                            warnings.warn("'TYPE' parameter not in table: ignoring in query", UserWarning)
                        elif exactMatch:
                            newtokens.append('(TYPE == "{}")'.format(tokens[i+2]))
                        else:
                            ttype = np.array([tokens[i+2] in a for a in table['TYPE']])
                            newtokens.append('(@ttype)')
                elif tokens[i].upper() == 'BINCOMP':
                    if 'BINCOMP' not in tab.keys():
                        warnings.warn("'BINCOMP' parameter not in table: ignoring in query", UserWarning)
                    elif exactMatch:
                        newtokens.append('(BINCOMP == "{}")'.format(tokens[i+2]))
                    else:
                        bincomp = np.array([tokens[i+2] in a for a in table['BINCOMP']])
                        newtokens.append('(@bincomp)')
                elif tokens[i].upper() == 'EXIST':
                    if tokens[i+2] not in tab.keys():
                        warnings.warn("'{}' does not exist for any pulsar".format(tokens[i+2]), UserWarning)
                        # create an empty DataFrame
                        tab = DataFrame(columns=table.keys())
                        break
                    else:
                        newtokens.append('({} != None)'.format(tokens[i+2]))
                elif tokens[i].upper() == 'ERROR':
                    if tokens[i+2]+'_ERR' not in tab.keys():
                        warnings.warn("Error value for '{}' not present: ignoring in query".format(tokens[i+2]), UserWarning)
                    else:
                        newtokens.append('{}_ERR'.format(tokens[i+2]))
            i += 2
        else:
            newtokens.append(tokens[i])

        i += 1

    # evaluate the expression
    try:
        newtab = tab.query(''.join(newtokens))
    except RuntimeError:
        raise RuntimeError("Could not parse the query")

    if isinstance(table, Table):
        # convert back to an astropy table
        tab = Table.from_pandas(tab)

        # re-add any units/types
        for key in table.colnames:
            tab.columns[key].unit = table.columns[key].unit
            tab[key] = tab[key].astype(table[key].dtype)

    return tab
