# -*- coding: utf-8 -*-
r"""
EnzymeModuleForm is a class for holding information regarding enzyme species.

The :class:`EnzymeModuleForm` class inherits and extends the
:class:`~.MassMetabolite` class. It is designed to represent various bound
states and conformations of the enzymes represented through the
:class:`~.EnzymeModule` class.

The enzyme specific attributes on the :class`EnzymeModuleForm` are the
following:

    * :attr:`~EnzymeModuleForm.enzyme_module_id`
    * :attr:`~EnzymeModuleForm.bound_catalytic`
    * :attr:`~EnzymeModuleForm.bound_effectors`

The :class:`EnzymeModuleForm` contains the attributes
:attr:`~EnzymeModuleForm.bound_catalytic` and
:attr:`~EnzymeModuleForm.bound_effectors`, designed to hold the
:class:`~.MassMetabolite`\ (s) that could be bound to the catalytic site and
regularory sites. There are no differences between the attributes other than
their name.

Some other important points about the :class:`EnzymeModuleForm` include:

    * If the :attr:`name` attribute is not set upon initializing, it is
      automatically generated using the enzyme specific attributes.
    * If the :attr:`formula` or charge attributes are not set upon
      initialization, it is inferred using the formulas and charges set on the
      :class:`~.MassMetabolite`\ (s) found in
      :attr:`~EnzymeModuleForm.bound_catalytic` and
      :attr:`~EnzymeModuleForm.bound_effectors`. A moiety is also included for
      the formula using the :attr:`~EnzymeModuleForm.enzyme_module_id`.
    * The purpose of the generated formula and charge is to ensure reactions
      remained mass and charge balanaced as metabolite species are bound and
      altered by the :class:`~.EnzymeModuleReaction`\ s of the
      :class:`~.EnzymeModule`.

"""
import re
from collections import defaultdict
from itertools import chain

from six import integer_types, iteritems

from mass.core.massmetabolite import MassMetabolite


class EnzymeModuleForm(MassMetabolite):
    r"""Class representing enzymatic species in an :class:`~.EnzymeModule`.

    Parameters
    ----------
    id_or_specie : str, MassMetabolite, EnzymeModuleForm
        A string identifier to associate with the enzymatic species, or an
        existing metabolite object. If an existing metabolite object is
        provided, a new :class:`EnzymeModuleForm` is instantiated with the
        same properties as the original metabolite.
    name : str
        A human readable name for the enzymatic species.
    formula : str
        Chemical formula associated with the enzymatic species.
    charge : float
        The charge number associated with the enzymatic species.
    compartment : str
        The compartment where the enzymatic species is located.
    fixed : bool
        Whether the enzymatic species concentration should remain at a
        fixed value. Default is ``False``.

    Attributes
    ----------
    enzyme_module_id : str
        The identifier of the associated :class:`~.EnzymeModule`.
    bound_catalytic : dict
        A dictionary representing the ligands bound to the enzyme's active
        site(s), with :class:`~.MassMetabolite`\ s as keys and number bound as
        values. If the final coefficient for a metabolite is 0, it is removed.
    bound_effectors : dict
        A dictionary representing the ligands bound to the enzyme's
        regulatory site(s), with :class:`~.MassMetabolite`\ s as keys and
        number bound as values. If the final coefficient for a metabolite is 0,
        it is removed.

    """

    def __init__(self, id_or_specie=None, name="", formula=None, charge=None,
                 compartment=None, fixed=False, enzyme_module_id="",
                 bound_catalytic=None, bound_effectors=None):
        """Initialize the EnzymeModuleForm."""
        # Initialize MassMetabolite parent class
        super(EnzymeModuleForm, self).__init__(
            id_or_specie=id_or_specie, name=name, formula=formula,
            charge=charge, compartment=compartment, fixed=fixed)
        if isinstance(id_or_specie, EnzymeModuleForm):
            # Instiantiate a new EnzymeModuleForm with state identical to
            # the provided EnzymeModuleForm object.
            self.__dict__.update(id_or_specie.__dict__)
        else:
            # Set the id of the enzyme represented by the EnzymeModuleForm
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
                for attr, val in zip(["formula", "charge"], [formula, charge]):
                    if val is None:
                        val = self.__class__.__dict__[
                            "get_species_" + attr](self)
                    setattr(self, attr, val)

    @property
    def bound_catalytic(self):
        r"""Get or set ligands bound to the enzyme's active site(s).

        Notes
        -----
        Assigning a dict to this property updates the current dict of
        ligands bound at the active site(s) with the new values.

        Parameters
        ----------
        value : dict
            Dictionary where keys are :class:`~.MassMetabolite` and values are
            the number currently bound to the active site(s).
            An empty ``dict`` will reset the bound ligands.

        """
        return getattr(self, "_bound_catalytic")

    @bound_catalytic.setter
    def bound_catalytic(self, value):
        """Set the dictionary of ligands bound to the active site(s)."""
        self._set_bound_dict("bound_catalytic", value)

    @property
    def bound_effectors(self):
        r"""Get or set ligands bound to the enzyme's regulatory site(s).

        Notes
        -----
        Assigning a dict to this property updates the current dict of
        ligands bound at the regulatory site(s) with the new values.

        Parameters
        ----------
        value : dict
            Dictionary where keys are :class:`~.MassMetabolite` and values are
            the number currently bound to the regulatory site(s).
            An empty ``dict`` will reset the bound ligands.

        """
        return getattr(self, "_bound_effectors")

    @bound_effectors.setter
    def bound_effectors(self, value):
        """Set the dictionary of ligands bound to the active site(s)."""
        self._set_bound_dict("bound_effectors", value)

    def generate_enzyme_module_form_name(self, update_enzyme=False):
        """Generate a name for the enzymatic species based on bound ligands.

        Notes
        -----
        * The :attr:`~EnzymeModuleForm.enzyme_module_id`,
          :attr:`~EnzymeModuleForm.bound_catalytic`, and
          :attr:`~EnzymeModuleForm.bound_effectors` attributes are used in
          generating the name.
        * If the :attr:`~EnzymeModuleForm.enzyme_module_id` attributes are not
          set, the string ``'Enzyme'`` will be used in its place.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`name` attribute of the enzymatic
            species in addition to returning the generated name.
            Default is ``False``.

        Returns
        -------
        str
            String representing the name of the :class:`EnzymeModuleForm`.

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

    def get_species_formula(self, update_enzyme=False):
        """Generate the chemical formula for the enzymatic species.

        Notes
        -----
        The :attr:`~EnzymeModuleForm.enzyme_module_id`,
        :attr:`~EnzymeModuleForm.bound_catalytic`, and
        :attr:`~EnzymeModuleForm.bound_effectors` attributes are used in
        generating the name.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`formula` attribute of the enzymatic
            species in addition to returning the generated formula.
            Default is ``False``.

        Returns
        -------
        str
            String representing the formula of the :class:`EnzymeModuleForm`.

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

    def get_species_charge(self, update_enzyme=False):
        """Generate the charge for the enzymatic species.

        Notes
        -----
        The attr:`~EnzymeModuleForm.bound_catalytic`, and
        :attr:`~EnzymeModuleForm.bound_effectors` attributes are used in
        generating the name.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`charge` attribute of the enzymatic
            species in addition to returning the generated charge.
            Default is ``False``.

        Returns
        -------
        float
            Value representing the charge of the :class:`EnzymeModuleForm`.

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
        """Set a dictionary of metabolites for a site.

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
        """Set the id of the EnzymeModuleForm to the associated MassModel.

        Warnings
        --------
        This method is intended for internal use only.

        """
        super(EnzymeModuleForm, self)._set_id_with_model(value)
        self.model.enzyme_module_forms._generate_index()

    def _repair_bound_obj_pointers(self):
        """Repair object pointer for metabolites in bound dict attributes.

        Requires a model to be associated with the EnzymeModuleForm.

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
        """HTML representation of the overview for the EnzymeModuleForm.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return """
        <table>
            <tr>
                <td><strong>EnzymeModuleForm identifier</strong></td>
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
        <table>""".format(
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


__all__ = ("EnzymeModuleForm",)
