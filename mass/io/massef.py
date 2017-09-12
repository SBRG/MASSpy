# -*- coding: utf-8 -*-

# Compatibility with Python 2.7
from __future__ import absolute_import

# Import necesary packages
import re
import io
from warnings import warn
from six import iteritems

# from mass
from mass.core.massmodel import MassModel
from mass.core.massreaction import MassReaction
from mass.core.massmetabolite import MassMetabolite

# Class begins
## Precompiled regular expressions
_rxn_id_finder = re.compile("(\w+):")
rxn_file_re = re.compile("(\S+)\_reactions")
kf_file_re = re.compile("(\S+)\_kf|(\S+)\_forward_rate_constant")
Keq_file_re = re.compile("(\S+)\_Keq|(\S+)\_equilibrium_constant")
kr_file_re = re.compile("(\S+)\_kr|(\S+)\_reverse_rate_constant")

def create_enzyme_model(rxn_file, kf_file=None, Keq_file=None, kr_file=None):
    """Create an enzyme model from text files

    The name for the file containing reactions  must be (MODEL_ID)_reactions,
    where the MODEL_ID has no spaces. Each reaction begins on a newline in the
    text file, and must be of a format readable by massmodel.string_to_mass.
    For example:

    Filename: 'MODEL_reactions'
        RID1: s[MET1, **kwargs] + s[MET2, **kwargs] <=>  s[MET3, **kwargs]
        RID2: s[MET3, **kwargs] + s[MET4, **kwargs] <=>  s[MET5, **kwargs]

    To import values for  kf, Keq, and kr values, the filename must be of the
    format '(MODEL_ID)_(argument), where each RID : value pair begins on a
    newline in the textfile. For example:

    Filename: 'MODEL_Keqs'
        RID1: 1
        RID2: 200.123

    Parameters
    ----------
    reaction_file : string
        Name of the file containing the reactions. Must adhere to the format of
        mass.massmodel.string_to_mass
    kf_file : string or None
        Name of the file containing the forward rate constants. Filename must
        be of format: '(MODEL_ID)_kf' or '(MODEL_ID)_forward_rate_constant'
    Keq_file : string or None
        Name of the file containing the equilibrium constants. Filename must
        be of format: '(MODEL_ID)_Keq' or '(MODEL_ID)_equilibrium_constant'
    kr_file : string or None
        Name of the file containing the reverse rate constants. Filename must
        be of format: '(MODEL_ID)_kr' or '(MODEL_ID)_reverse_rate_constant'

    See Also
	--------
    massmodel.string_to_mass
        For reactions from strings.
    """
    # Check file input
    if not rxn_file_re.search(rxn_file):
        raise ValueError("Could not find the file of reactions.Please ensure"
                " that the filename is of the format: '(MODEL_ID)_reactions'"
                " where the MODEL_ID has no spaces.")
    else:
        model_id = rxn_file_re.search(rxn_file).group(1)
        new_model = MassModel(model_id)

    # Get reaction strings
    with io.open(rxn_file, 'r') as file:
        reaction_strings = file.readlines()
    # Remove new line characters
    reaction_strings = [rxn.strip('\n') for rxn in reaction_strings]
    # Create metabolites and reactions for the model
    new_model.string_to_mass(reaction_strings)

    # Add other parameters if files given.
    if kf_file is not None:
        new_model = import_forward_rate_constants(new_model, kf_file)
    if Keq_file is not None:
        new_model = import_equilibrium_constants(new_model, Keq_file)
    if kr_file is not None:
        new_model = import_reverse_rate_constants(new_model, kr_file)

    return new_model

def import_forward_rate_constants(model, kf_file):
    """Import forward rate constant values for an enzyme model

    Note: Reactions must already exist in the massmodel.

    Parameters
    ----------
    model : mass.MassModel
        The enzyme massmodel.
    kf_file : string or None
        Name of the file containing the forward rate constants. Filename must
        be of format: '(MODEL_ID)_kf' or '(MODEL_ID)_forward_rate_constant'
    """
    # Check input
    if not kf_file_re.search(kf_file):
        raise ValueError("Could not find the file of forward rate constants."
                    "Please ensure that the filename is of the format: "
                    "'(MODEL_ID)_kf' or '(MODEL_ID)_forward_rate_constants' "
                    "where the MODEL_ID has no spaces.")
    if not isinstance(model, MassModel):
        raise TypeError("model must be a massmodel")
    else:
        model_id = kf_file_re.search(kf_file).group(1)
        if not re.match(model.id, model_id):
            raise ValueError("Filename model ID does not match given model ID")

    with io.open(kf_file, 'r') as file:
        params_and_values = file.readlines()

    params_and_values = [param.strip('\n') for param in params_and_values]
    for param_string in params_and_values:
        try:
            res = _rxn_id_finder.search(param_string)
            rxn_id = res.group(1)
            param_string = param_string[res.end():]
        except:
            ValueError("Could not find an ID for '%s'" % param_string)

        reaction = model.reactions.get_by_id(rxn_id)
        reaction.kf = float(param_string)

    return model

def import_equilibrium_constants(model, Keq_file):
    """Import equilibrium constant values for an enzyme model

    Note: Reactions must already exist in the massmodel.

    Parameters
    ----------
    model : mass.MassModel
        The enzyme massmodel.
    Keq_file : string or None
        Name of the file containing the equilibrium constants. Filename must
        be of format: '(MODEL_ID)_Keq' or '(MODEL_ID)_equilibrium_constant'
    """
    if not Keq_file_re.search(Keq_file):
        raise ValueError("Could not find the file of equilibrium constants."
                    "Please ensure that the filename is of the format: "
                    "'(MODEL_ID)_Keq' or '(MODEL_ID)_equilibrium_constants' "
                    "where the MODEL_ID has no spaces.")

    if not isinstance(model, MassModel):
        raise TypeError("model must be a massmodel")
    else:
        model_id = Keq_file_re.search(Keq_file).group(1)
        if not re.match(model.id, model_id):
            raise ValueError("Filename model ID does not match given model ID")

    with io.open(Keq_file, 'r') as file:
        params_and_values = file.readlines()

    params_and_values = [param.strip('\n') for param in params_and_values]
    for param_string in params_and_values:
        try:
            res = _rxn_id_finder.search(param_string)
            rxn_id = res.group(1)
            param_string = param_string[res.end():]
        except:
            ValueError("Could not find an ID for '%s'" % param_string)

        reaction = model.reactions.get_by_id(rxn_id)
        reaction.Keq = float(param_string)

    return model

def import_reverse_rate_constants(model, kr_file):
    """Import reverse rate constant values for an enzyme model

    Note: Reactions must already exist in the massmodel.

    Parameters
    ----------
    model : mass.MassModel
        The enzyme massmodel.
    kr_file : string or None
        Name of the file containing the reverse rate constants. Filename must
        be of format: '(MODEL_ID)_kr' or '(MODEL_ID)_reverse_rate_constants'
    """
    if not kr_file_re.search(kr_file):
        raise ValueError("Could not find the file of reverse rate constants."
                    "Please ensure that the filename is of the format: "
                    "'(MODEL_ID)_kr' or '(MODEL_ID)_reverse_rate_constants' "
                    "where the MODEL_ID has no spaces.")
    if not isinstance(model, MassModel):
        raise TypeError("model must be a massmodel")
    else:
        model_id = kr_file_re.search(kr_file).group(1)
        if not re.match(model.id, model_id):
            raise ValueError("Filename model ID does not match given model ID")

    with io.open(kr_file, 'r') as file:
        params_and_values = file.readlines()

    params_and_values = [param.strip('\n') for param in params_and_values]
    for param_string in params_and_values:
        try:
            res = _rxn_id_finder.search(param_string)
            rxn_id = res.group(1)
            param_string = param_string[res.end():]
        except:
            ValueError("Could not find an ID for '%s'" % param_string)

        reaction = model.reactions.get_by_id(rxn_id)
        reaction.kr = float(param_string)

    return model
