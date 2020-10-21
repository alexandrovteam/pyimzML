from datetime import datetime
from warnings import warn

from .uo import terms as uo_terms
from .ms import terms as ms_terms
from .ims import terms as ims_terms

all_terms = {}
all_terms.update(uo_terms)
all_terms.update(ms_terms)
all_terms.update(ims_terms)

DTYPE_MAPPING = {
    'xsd:string': str,
    'xsd:anyURI': str,
    'xsd:float': float,
    'xsd:double': float,
    'xsd:decimal': float,
    'xsd:nonNegativeFloat': float,
    'xsd:int': int,
    'xsd:integer': int,
    'xsd:positiveInteger': int,
    'xsd:nonNegativeInteger': int,
    'xsd:boolean': bool,
    'xsd:dateTime': datetime,
}

FLOAT_TYPE_FIX_MAPPING = {
'32-bit float': 'MS:1000521',
'64-bit float': 'MS:1000523',
}


def convert_xml_value(dtype, value):
    try:
        if dtype is not None:
            return DTYPE_MAPPING[dtype](value)
        elif value is None or value is '':
            # Many cv_params are flags and have either a None or empty-string value.
            # Replace their value with True in these cases, so their existance isn't so ambiguous.
            return True
        else:
            return value
    except KeyError:
        return value
    except ValueError:
        return None


def convert_term_name(accession):
    return all_terms.get(accession, (accession, None))[0]


def convert_cv_param(accession, value):
    """
    Looks up a term by accession number, and convert the provided value to the expected type.
    """
    name, dtype = all_terms.get(accession, (accession, None))
    converted_value = convert_xml_value(dtype, value)
    return converted_value


def lookup_and_convert_cv_param(accession, raw_name, value, unit_accession=None):
    """
    Looks up a term by accession number, and returns the term name, its value converted into
    the expected datatype, and the unit name (if a unit accession number is also given).
    """
    name, dtype = all_terms.get(accession, (raw_name or accession, None))
    converted_value = convert_xml_value(dtype, value)
    unit_name = all_terms.get(unit_accession, (unit_accession, None))[0]

    if accession not in all_terms:
        warn('Unrecognized accession in <cvParam>: %s (name: "%s").' % (accession, raw_name))
    elif name != raw_name:
        #
        if accession in FLOAT_TYPE_FIX_MAPPING.values() and name in FLOAT_TYPE_FIX_MAPPING:
            fixed_accession = FLOAT_TYPE_FIX_MAPPING[raw_name]
            warn(
                'Accession %s found with incorrect name "%s" (expected "%s"). '
                'This is a known issue with some imzML conversion software - updating accession '
                'to %s.' % (accession, raw_name, name, fixed_accession)
            )
            accession = fixed_accession
            name = raw_name
        else:
            warn(
                'Accession %s found with incorrect name "%s". Updating name to "%s".'
                % (accession, raw_name, name)
            )

    return accession, name, converted_value, unit_name



