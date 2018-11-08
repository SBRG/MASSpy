# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

from mass.exceptions import MassSBMLError as _MassSBMLError

from mass.io.json import (
    load_json_model, save_json_model)

try:
    from mass.io.sbml import (
        write_sbml_model, read_sbml_model)
except _MassSBMLError:
    pass
