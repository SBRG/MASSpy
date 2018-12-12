# -*- coding: utf-8 -*-
"""This module is used to configure global variables.

Global variables include the logger for logging capabilities, and the global
zero tolerance value.
"""

from __future__ import absolute_import

import logging

log = logging.getLogger(__name__)
ZERO_TOLERANCE = 1e-8
