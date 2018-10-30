# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

from mass.exceptions import MassSBMLError as _MassSBMLError

from mass.io.json import (
    write_json_model, read_json_model)
from mass.io.enzymes import (
    create_enzyme_model)

try:
    from mass.io.sbml import (
        write_sbml_model, read_sbml_model)
except _MassSBMLError:
    pass
