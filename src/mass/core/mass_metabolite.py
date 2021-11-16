# -*- coding: utf-8 -*-
r"""MassMetabolite is a class for holding information regarding metabolites.

The :class:`MassMetabolite` class inherits and extends the
:class:`~cobra.core.metabolite.Metabolite` class in :mod:`cobra`. It contains
additional information required for simulations and other :mod:`mass`
functions and workflows.

Some key differences between the
:class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>` and the
:class:`mass.MassMetabolite <mass.core.mass_metabolite.MassMetabolite>` are
listed below:

    * Unlike the :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>`
      the ``id`` initialization argument has been replaced with
      ``id_or_specie`` to be clear in that a new :class:`MassMetabolite` can
      be instantiated with the same properties of a previously created
      metabolite by passing the metabolite object to this argument.

    * The :attr:`formula`\ s can have moieties by placing them within square
      brackets (e.g. ``"[ENZYME]"``)

"""
import re
from warnings import warn

from cobra.core.metabolite import Metabolite, element_re, elements_and_molecular_weights
from cobra.util.context import resettable
from cobra.util.util import format_long_string
from six import iteritems

from mass.util.expressions import generate_ode
from mass.util.util import ensure_non_negative_value, get_public_attributes_and_methods


class MassMetabolite(Metabolite):
    """Class for holding information regarding a metabolite.

    Parameters
    ----------
    id_or_specie : str, ~cobra.core.metabolite.Metabolite, MassMetabolite
        A string identifier to associate with the metabolite, or an existing
        metabolite object. If an existing metabolite object is
        provided, a new :class:`MassMetabolite` is instantiated with the same
        properties as the original metabolite.
    name : str
        A human readable name for the metabolite.
    formula : str
        Chemical formula associated with the metabolite (e.g. H2O).
    charge : float
        The charge number associated with the metabolite.
    compartment : str
        Compartment of the metabolite.
    fixed : bool
        Whether the metabolite concentration should remain at a fixed value.
        Default is ``False``.

    """

    def __init__(
        self,
        id_or_specie=None,
        name="",
        formula=None,
        charge=None,
        compartment=None,
        fixed=False,
    ):
        """Initialize the MassMetabolite."""
        super(MassMetabolite, self).__init__(
            id=str(id_or_specie),
            formula=formula,
            name=name,
            charge=charge,
            compartment=compartment,
        )
        if isinstance(id_or_specie, (Metabolite, MassMetabolite)):
            # Instiantiate a new MassMetabolite with state identical to
            # the provided Metabolite/MassMetabolite object.
            self.__dict__.update(id_or_specie.__dict__)
            self._reaction = set()
            self._model = None

        # If is not a MassMetabolite object, initialize additional attributes
        if not isinstance(id_or_specie, MassMetabolite):
            # Whether the concentration of the metabolite is fixed.
            self._fixed = fixed
            # The initial condition of the metabolite.
            self._initial_condition = None

    # Public
    @property
    def elements(self):
        """Get or set ``dict`` of elements in the :attr:`formula`.

        Parameters
        ----------
        elements_dict : dict
            A ``dict`` representing the elements of the chemical formula where
            keys are elements and values are the amount.

        Notes
        -----
        * Enzyme and macromolecule moieties can be recognized by enclosing them
          in brackets (e.g. [ENZYME]) when defining the chemical formula.
          They are treated as one entity and therefore are only counted once.
        * Overrides
          :meth:`~cobra.core.metabolite.Metabolite.elements` of the
          :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>`
          to allow for the use of moieties.

        """
        tmp_formula = self.formula
        if tmp_formula is None:
            return {}
        if "*" in tmp_formula:
            warn("invalid character '*' found in formula '{0}'".format(self.formula))
            tmp_formula = tmp_formula.replace("*", "")
        composition = {}
        # Identify any moieties
        while "[" in tmp_formula and "]" in tmp_formula:
            s = tmp_formula.index("[")
            e = tmp_formula.index("]") + 1
            moiety = tmp_formula[s:e]
            composition[moiety] = 1
            tmp_formula = tmp_formula.replace(tmp_formula[s:e], "")

        # Count elements
        for (element, count) in element_re.findall(tmp_formula):
            if count == "":
                count = 1
            else:
                try:
                    count = float(count)
                    if count == int(count):
                        count = int(count)
                    else:
                        warn(
                            "{0} is not an integer (in formula {1})".format(
                                count, self.formula
                            )
                        )
                except ValueError:
                    warn(
                        "failed to parse {0} (in formula {1})".format(
                            count, self.formula
                        )
                    )
                    return None
            if element in composition:
                composition[element] += count
            else:
                composition[element] = count

        return composition

    @elements.setter
    def elements(self, elements_dict):
        """Set the formula using a ``dict`` of elements."""

        def stringify(element, number):
            return element if number == 1 else element + str(number)

        self.formula = "".join(
            stringify(e, n) for e, n in sorted(iteritems(elements_dict))
        )

    @property
    def initial_condition(self):
        """Get or set the initial condition of the metabolite.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Notes
        -----
        Initial conditions of metabolites cannot be negative.

        Parameters
        ----------
        value : float
            A non-negative number for the concentration of the metabolite.

        Raises
        ------
        ValueError
            Occurs when trying to set a negative value.

        """
        return getattr(self, "_initial_condition")

    @initial_condition.setter
    @resettable
    def initial_condition(self, initial_condition):
        """Set the initial condition of the metabolite."""
        initial_condition = ensure_non_negative_value(initial_condition)
        setattr(self, "_initial_condition", initial_condition)

    @property
    def fixed(self):
        """Get or set whether the metabolite remains constant.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        Parameters
        ----------
        value : bool
            Whether the metabolite should remain constant, meaning that the
            metabolite ODE is 0.

        """
        return getattr(self, "_fixed")

    @fixed.setter
    @resettable
    def fixed(self, value):
        """Set whether the metabolite remains constant.

        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.

        """
        if not isinstance(value, bool):
            raise TypeError("value must be a bool")

        setattr(self, "_fixed", value)

    @property
    def ordinary_differential_equation(self):
        """Return a :mod:`sympy` expression of the metabolite's associated ODE.

        Will return ``None`` if metabolite is not associated with a
        :class:`~.MassReaction`, or ``0``. if the :attr:`fixed` attribute
        is set as ``True``.
        """
        return generate_ode(self)

    @property
    def formula_weight(self):
        """Calculate and return the formula weight of the metabolite.

        Does not consider any moieties enclosed in brackets.

        Notes
        -----
        Overrides :meth:`~cobra.core.metabolite.Metabolite.formula_weight`
        of the :class:`cobra.Metabolite <cobra.core.metabolite.Metabolite>`
        to allow for the use of moieties.

        """
        # Remove moieties
        element_dict = {
            k: v
            for k, v in iteritems(self.elements)
            if not (k.startswith("[") and k.endswith("]"))
        }
        # Calculate formula weight
        try:
            return sum(
                [
                    count * elements_and_molecular_weights[element]
                    for element, count in sorted(iteritems(element_dict))
                ]
            )
        except KeyError as e:
            warn("The element {0} does not appear in the periodic table".format(e))

    @property
    def model(self):
        """Return the :class:`.MassModel` associated with the metabolite."""
        return getattr(self, "_model")

    @property
    def ic(self):
        """Alias for the :attr:`initial_condition`."""
        return self.initial_condition

    @ic.setter
    def ic(self, value):
        """Alias for the :attr:`initial_condition`."""
        self.initial_condition = value

    @property
    def ode(self):
        """Alias for the :attr:`ordinary_differential_equation`."""
        return self.ordinary_differential_equation

    def _remove_compartment_from_id_str(self):
        """Remove the compartment from the ID str of the metabolite.

        Warnings
        --------
        This method is intended for internal use only.

        """
        met_id_str = str(self)
        if self.compartment:
            met_id_str = re.sub("_" + self.compartment + "$", "", met_id_str)

        return met_id_str

    def _repr_html_(self):
        """HTML representation of the overview for the MassMetabolite.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return """
        <table>
            <tr>
                <td><strong>MassMetabolite identifier</strong></td>
                <td>{id}</td>
            </tr><tr>
                <td><strong>Name</strong></td>
                <td>{name}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{address}</td>
            </tr><tr>
                <td><strong>Formula</strong></td>
                <td>{formula}</td>
            </tr><tr>
                <td><strong>Compartment</strong></td>
                <td>{compartment}</td>
            </tr><tr>
                <td><strong>Initial Condition</strong></td>
                <td>{fixed}{ic}</td>
            </tr><tr>
                <td><strong>In {n_reactions} reaction(s)</strong></td>
                <td>{reactions}</td>
            </tr>
        </table>""".format(
            id=self.id,
            name=format_long_string(self.name),
            formula=self.formula,
            address="0x0%x" % id(self),
            compartment=self.compartment,
            fixed="Fixed at " if self.fixed else "",
            ic=self.initial_condition,
            n_reactions=len(self.reactions),
            reactions=format_long_string(", ".join(r.id for r in self.reactions), 200),
        )

    def __dir__(self):
        """Override default dir() implementation to list only public items.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return get_public_attributes_and_methods(self)


__all__ = ("MassMetabolite",)
