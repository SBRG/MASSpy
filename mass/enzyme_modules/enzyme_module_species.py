# -*- coding: utf-8 -*-
r"""
EnzymeModuleSpecies is a class for holding information regarding enzyme module species.

The :class:`EnzymeModuleSpecies` class inherits and extends the
:class:`~.MassMetabolite` class. It is designed to represent various bound
states and conformations of the enzymes represented through the
:class:`~.EnzymeModule` class.

The enzyme specific attributes on the :class:`EnzymeModuleSpecies` are the
following:

    * :attr:`~EnzymeModuleSpecies.enzyme_module_id`
    * :attr:`~EnzymeModuleSpecies.bound_catalytic`
    * :attr:`~EnzymeModuleSpecies.bound_effectors`

The :class:`EnzymeModuleSpecies` contains the attributes
:attr:`~EnzymeModuleSpecies.bound_catalytic` and
:attr:`~EnzymeModuleSpecies.bound_effectors`, designed to hold the
:class:`~.MassMetabolite`\ (s) that could be bound to the catalytic site and
regularory sites. There are no differences between the attributes other than
their name.

Some other important points about the :class:`EnzymeModuleSpecies` include:

    * If the :attr:`name` attribute is not set upon initializing, it is
      automatically generated using the enzyme specific attributes.

    * If the :attr:`formula` or charge attributes are not set upon
      initialization, it is inferred using the formulas and charges set on the
      :class:`~.MassMetabolite`\ (s) found in
      :attr:`~EnzymeModuleSpecies.bound_catalytic` and
      :attr:`~EnzymeModuleSpecies.bound_effectors`. A moiety is also included
      for the formula using the :attr:`~EnzymeModuleSpecies.enzyme_module_id`.

    * The purpose of the generated formula and charge is to ensure reactions
      remained mass and charge balanaced as metabolite species are bound and
      altered by the :class:`~.EnzymeModuleReaction`\ s of the
      :class:`~.EnzymeModule`.

"""  # noqa
import re
from collections import defaultdict
from itertools import chain

from six import integer_types, iteritems

from mass.core.mass_metabolite import MassMetabolite


class EnzymeModuleSpecies(MassMetabolite):
    r"""Class representing an enzyme species of an :class:`~.EnzymeModule`.

    Accepted ``kwargs`` are passed to the initialization method of the base
    class, :class:`.MassMetabolite`.

    Parameters
    ----------
    id_or_specie : str, MassMetabolite, EnzymeModuleSpecies
        A string identifier to associate with the enzyme module species, or an
        existing metabolite object. If an existing metabolite object is
        provided, a new :class:`EnzymeModuleSpecies` is instantiated with the
        same properties as the original metabolite.
    enzyme_module_id : str
        The identifier of the associated :class:`~.EnzymeModule`.
    bound_catalytic : dict
        A ``dict`` representing the ligands bound to the enzyme's active
        site(s), with :class:`~.MassMetabolite`\ s as keys and number bound as
        values. If the final coefficient for a metabolite is 0, it is removed.
    bound_effectors : dict
        A ``dict`` representing the ligands bound to the enzyme's
        regulatory site(s), with :class:`~.MassMetabolite`\ s as keys and
        number bound as values. If the final coefficient for a metabolite
        is 0, it is removed.
    **kwargs
        name :
            ``str`` representing a human readable name for the enzyme module
            species.
        formula :
            ``str`` representing a chemical formula associated with the
            enzyme module species.
        charge :
            ``float`` representing the charge number associated with the
            enzyme module species.
        compartment :
            ``str`` representing the compartment where the enzyme module
            species is located.
        fixed :
            ``bool`` indicating whether the enzyme module species
            concentration should remain at a fixed value. Default is ``False``.

    """

    def __init__(self, id_or_specie=None, enzyme_module_id="",
                 bound_catalytic=None, bound_effectors=None, **kwargs):
        """Initialize the EnzymeModuleSpecies."""
        # Initialize MassMetabolite parent class
        super(EnzymeModuleSpecies, self).__init__(
            id_or_specie=id_or_specie,
            name=kwargs.get("name", ""),
            formula=kwargs.get("formula", None),
            charge=kwargs.get("charge", None),
            compartment=kwargs.get("compartment", None),
            fixed=kwargs.get("fixed", False))
        if isinstance(id_or_specie, EnzymeModuleSpecies):
            # Instiantiate a new EnzymeModuleSpecies with state identical to
            # the provided EnzymeModuleSpecies object.
            self.__dict__.update(id_or_specie.__dict__)
        else:
            # Set the id of the enzyme represented by the EnzymeModuleSpecies
            self.enzyme_module_id = enzyme_module_id

            # Set metabolites bound to active site(s) of the enzyme form
            self._bound_catalytic = {}
            self.bound_catalytic = bound_catalytic

            # Set metabolites bound to regulatory site(s) of the enzyme form
            self._bound_effectors = {}
            self.bound_effectors = bound_effectors

            if not isinstance(id_or_specie, MassMetabolite):
                # Set formula, charge, and compartment attributes if
                # if a MassMetabolite was not used to initialize object.
                for attr in ["formula", "charge"]:
                    val = kwargs.get(attr, None)
                    if val is None:
                        val = self.__class__.__dict__[
                            "generate_species_" + attr](self)
                    setattr(self, attr, val)

    @property
    def bound_catalytic(self):
        r"""Get or set ligands bound to the enzyme's active site(s).

        Notes
        -----
        Assigning a ``dict`` to this property updates the current ``dict`` of
        ligands bound at the active site(s) with the new values.

        Parameters
        ----------
        value : dict
            A ``dict`` where keys are :class:`~.MassMetabolite` and values
            are the number currently bound to the active site(s).
            An empty ``dict`` will reset the bound ligands.

        """
        return getattr(self, "_bound_catalytic")

    @bound_catalytic.setter
    def bound_catalytic(self, value):
        """Set the ``dict`` of ligands bound to the active site(s)."""
        self._set_bound_dict("bound_catalytic", value)

    @property
    def bound_effectors(self):
        r"""Get or set ligands bound to the enzyme's regulatory site(s).

        Notes
        -----
        Assigning a ``dict`` to this property updates the current ``dict`` of
        ligands bound at the regulatory site(s) with the new values.

        Parameters
        ----------
        value : dict
            A ``dict`` where keys are :class:`~.MassMetabolite` and values
            are the number currently bound to the regulatory site(s).
            An empty ``dict`` will reset the bound ligands.

        """
        return getattr(self, "_bound_effectors")

    @bound_effectors.setter
    def bound_effectors(self, value):
        """Set the ``dict`` of ligands bound to the active site(s)."""
        self._set_bound_dict("bound_effectors", value)

    def generate_enzyme_module_species_name(self, update_enzyme=False):
        """Generate name for the enzyme module species based on bound ligands.

        Notes
        -----
        * The :attr:`~EnzymeModuleSpecies.enzyme_module_id`,
          :attr:`~EnzymeModuleSpecies.bound_catalytic`, and
          :attr:`~EnzymeModuleSpecies.bound_effectors` attributes are used in
          generating the name.
        * If the :attr:`~EnzymeModuleSpecies.enzyme_module_id` attributes are
          not set, the string ``'Enzyme'`` will be used in its place.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`name` attribute of the enzyme
            module species in addition to returning the generated name.
            Default is ``False``.

        Returns
        -------
        str
            String representing the name of the :class:`EnzymeModuleSpecies`.

        """
        if self.enzyme_module_id:
            name = self.enzyme_module_id
        else:
            name = "Enzyme"

        # Add the ligands bound to the active site(s)
        catalytic_str = "-".join([met._remove_compartment_from_id_str()
                                  for met in self.bound_catalytic])
        if catalytic_str:
            catalytic_str = "-" + catalytic_str + " complex"
        name += catalytic_str

        # Add the ligands bound to the effector site(s)
        effector_str = _make_bound_attr_str_repr(self.bound_effectors)
        if effector_str:
            name += "; " + effector_str

        if update_enzyme:
            self.name = name

        return name

    def generate_species_formula(self, update_enzyme=False):
        """Generate the chemical formula for the enzyme module species.

        This function is primarily utilized for keeping reactions between
        :class:`EnzymeModuleSpecies` mass and charge balanced.

        Notes
        -----
        The :attr:`~EnzymeModuleSpecies.enzyme_module_id`,
        :attr:`~EnzymeModuleSpecies.bound_catalytic`, and
        :attr:`~EnzymeModuleSpecies.bound_effectors` attributes are used in
        generating the formula.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`formula` attribute of the enzyme
            module species in addition to returning the generated formula.
            Default is ``False``.

        Returns
        -------
        str
            String representing the formula of the
            :class:`EnzymeModuleSpecies`.

        """
        formula = ""
        if self.enzyme_module_id:
            moiety = self.enzyme_module_id.lower()
            if not re.match('^[A-Z]+[a-z]+$', moiety):
                moiety = re.sub("[0-9]", "", moiety)

            formula += "[" + moiety + "]"

        total_elements = defaultdict(list)
        for dictionary in [self.bound_catalytic, self.bound_effectors]:
            if not dictionary:
                continue
            elem_iters = [iteritems({k: v * num_bound
                                     for k, v in iteritems(met.elements)})
                          for met, num_bound in iteritems(dictionary)]

            for k, v in chain(*elem_iters):
                total_elements[k].append(v)

        total_elements = {k: sum(v) for k, v in iteritems(total_elements)}
        if total_elements and formula:
            formula += "-"
        for k, v in iteritems(total_elements):
            formula += k + str(v)

        if update_enzyme:
            self.formula = formula

        return formula

    def generate_species_charge(self, update_enzyme=False):
        """Generate the charge for the enzyme module species.

        This function is primarily utilized for keeping reactions between
        :class:`EnzymeModuleSpecies` mass and charge balanced.

        Notes
        -----
        The attr:`~EnzymeModuleSpecies.bound_catalytic`, and
        :attr:`~EnzymeModuleSpecies.bound_effectors` attributes are used
        in generating the charge.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`charge` attribute of the enzyme
            module species in addition to returning the generated charge.
            Default is ``False``.

        Returns
        -------
        float
            Value representing the charge of the :class:`EnzymeModuleSpecies`.

        """
        charge = 0
        for dictionary in [self.bound_catalytic, self.bound_effectors]:
            if not dictionary:
                continue
            for met, num_bound in iteritems(dictionary):
                charge += met.charge * num_bound

        if update_enzyme:
            self.charge = charge

        return charge

    # Internal
    def _set_bound_dict(self, attribute, value):
        """Set a ``dict`` of metabolites for a site.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if attribute not in {"bound_catalytic", "bound_effectors"}:
            raise ValueError("Must be either {'bound_catalytic', "
                             "'bound_effectors'}")

        if not isinstance(value, dict) and value is not None:
            raise TypeError(attribute + " must be a dict")

        if value:
            bound_dict = {}
            for met, num_bound in iteritems(value):
                if not isinstance(met, MassMetabolite) or \
                   not isinstance(num_bound, integer_types):
                    raise ValueError("value must be a dict where keys are "
                                     "MassMetabolites and values are ints")
                if num_bound != 0:
                    bound_dict[met] = num_bound
            getattr(self, "_" + attribute).update(bound_dict)
        else:
            setattr(self, "_" + attribute, {})

    def _set_id_with_model(self, value):
        """Set the id of the EnzymeModuleSpecies to the associated MassModel.

        Warnings
        --------
        This method is intended for internal use only.

        """
        super(EnzymeModuleSpecies, self)._set_id_with_model(value)
        self.model.enzyme_module_species._generate_index()

    def _repair_bound_obj_pointers(self):
        """Repair object pointer for metabolites in bound dict attributes.

        Requires a model to be associated with the EnzymeModuleSpecies.

        Warnings
        --------
        This method is intended for internal use only.

        """
        if self.model is not None:
            for attr in ["_bound_catalytic", "_bound_effectors"]:
                bound_dict = getattr(self, attr)
                try:
                    bound_dict = {
                        self.model.metabolites.get_by_id(str(met)): num
                        for met, num in iteritems(bound_dict)}
                except KeyError as e:
                    raise KeyError("'{0}' does not exist in model metabolites."
                                   .format(str(e)))
                # Add the metabolites into the module that don't already exist
                bound_dict = setattr(self, attr, bound_dict)

    def _repr_html_(self):
        """HTML representation of the overview for the EnzymeModuleSpecies.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return """
        <table>
            <tr>
                <td><strong>EnzymeModuleSpecies identifier</strong></td>
                <td>{id}</td>
            </tr><tr>
                <td><strong>Name</strong></td>
                <td>{name}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{address}</td>
            </tr><tr>
                <td><strong>Enzyme Module</strong></td>
                <td>{enzyme}</td>
            </tr><tr>
                <td><strong>Compartment</strong></td>
                <td>{compartment}</td>
            </tr><tr>
                <td><strong>Bound Catalytic</strong></td>
                <td>{catalytic}</td>
            </tr><tr>
                <td><strong>Bound Effectors</strong></td>
                <td>{effectors}</td>
            </tr><tr>
                <td><strong>Initial Condition</strong></td>
                <td>{ic}</td>
            </tr><tr>
                <td><strong>In {n_reactions} reaction(s)</strong></td>
                <td>{reactions}</td>
            </tr>
        </table>""".format(
            id=self.id, name=self.name, address='0x0%x' % id(self),
            enzyme=str(self.enzyme_module_id), compartment=self.compartment,
            catalytic=_make_bound_attr_str_repr(self.bound_catalytic),
            effectors=_make_bound_attr_str_repr(self.bound_effectors),
            ic=self._initial_condition, n_reactions=len(self.reactions),
            reactions=', '.join(r.id for r in self.reactions))


def _make_bound_attr_str_repr(attribute_dict):
    """Make the string representation for bound ligands.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return "; ".join([
        "{0} {1}".format(v, k)
        for k, v in iteritems(attribute_dict) if v != 0])


__all__ = ("EnzymeModuleSpecies",)
