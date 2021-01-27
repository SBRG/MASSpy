# -*- coding: utf-8 -*-
r"""
EnzymeModuleForm is a class for holding information regarding enzyme module forms.

The :class:`EnzymeModuleForm` class inherits and extends the
:class:`~.MassMetabolite` class. It is designed to represent various bound
states and conformations of the enzymes represented through the
:class:`~.EnzymeModule` class.

The enzyme specific attributes on the :class:`EnzymeModuleForm` are the
following:

    * :attr:`~EnzymeModuleForm.enzyme_module_id`
    * :attr:`~EnzymeModuleForm."bound_metabolites"`

The :class:`EnzymeModuleForm` contains the attribute
:attr:`~EnzymeModuleForm."bound_metabolites", designed to hold the
:class:`~.MassMetabolite`\ (s) that could be bound to the sites on the enzyme.

Some other important points about the :class:`EnzymeModuleForm` include:

    * If the :attr:`name` attribute is not set upon initializing, it is
      automatically generated using the enzyme specific attributes.

    * If the :attr:`formula` or charge attributes are not set upon
      initialization, it is inferred using the formulas and charges set on the
      :class:`~.MassMetabolite`\ (s) found in
      :attr:`~EnzymeModuleForm.bound_metabolites`. A moiety is also included
      for the formula using the :attr:`~EnzymeModuleForm.enzyme_module_id`.

    * The purpose of the generated formula and charge is to ensure reactions
      remained mass and charge balanaced as metabolite species are bound and
      altered by the :class:`~.EnzymeModuleReaction`\ s of the
      :class:`~.EnzymeModule`.

"""
import re
from collections import defaultdict
from itertools import chain

from six import integer_types, iteritems

from mass.core.mass_metabolite import MassMetabolite


class EnzymeModuleForm(MassMetabolite):
    r"""Class representing an enzyme forms of an :class:`~.EnzymeModule`.

    Accepted ``kwargs`` are passed to the initialization method of the base
    class, :class:`.MassMetabolite`.

    Parameters
    ----------
    id_or_specie : str, MassMetabolite, EnzymeModuleForm
        A string identifier to associate with the enzyme module forms, or an
        existing metabolite object. If an existing metabolite object is
        provided, a new :class:`EnzymeModuleForm` is instantiated with the
        same properties as the original metabolite.
    enzyme_module_id : str
        The identifier of the associated :class:`~.EnzymeModule`.
    bound_metabolites : dict
        A ``dict`` representing the ligands bound to the enzyme,
        with :class:`~.MassMetabolite`\ s or their identifiers as
        keys and the number bound as values.
    **kwargs
        name :
            ``str`` representing a human readable name for the enzyme module
            form.
        formula :
            ``str`` representing a chemical formula associated with the
            enzyme module form.
        charge :
            ``float`` representing the charge number associated with the
            enzyme module form.
        compartment :
            ``str`` representing the compartment where the enzyme module
            form is located.
        fixed :
            ``bool`` indicating whether the enzyme module form
            concentration should remain at a fixed value. Default is ``False``.

    """

    def __init__(
        self, id_or_specie=None, enzyme_module_id="", bound_metabolites=None, **kwargs
    ):
        """Initialize the EnzymeModuleForm."""
        # Initialize MassMetabolite parent class
        super(EnzymeModuleForm, self).__init__(
            id_or_specie=id_or_specie,
            name=kwargs.get("name", ""),
            formula=kwargs.get("formula", None),
            charge=kwargs.get("charge", None),
            compartment=kwargs.get("compartment", None),
            fixed=kwargs.get("fixed", False),
        )
        if isinstance(id_or_specie, EnzymeModuleForm):
            # Instiantiate a new EnzymeModuleForm with state identical to
            # the provided EnzymeModuleForm object.
            self.__dict__.update(id_or_specie.__dict__)
        else:
            # Set the id of the enzyme represented by the EnzymeModuleForm
            self.enzyme_module_id = enzyme_module_id

            # Set metabolites bound to site(s) of the enzyme form
            self._bound_metabolites = {}
            self.bound_metabolites = bound_metabolites

            if not isinstance(id_or_specie, MassMetabolite):
                # Set formula, charge, and compartment attributes if
                # if a MassMetabolite was not used to initialize object.
                for attr in ["formula", "charge"]:
                    val = kwargs.get(attr, None)
                    if val is None:
                        val = self.__class__.__dict__["generate_form_" + attr](self)
                    setattr(self, attr, val)

    @property
    def bound_metabolites(self):
        r"""Get or set metabolites bound to the enzyme's site(s).

        Notes
        -----
        Assigning a ``dict`` to this property updates the current ``dict`` of
        ligands bound at the enzyme site(s) with the new values.

        Parameters
        ----------
        value : dict
            A ``dict`` where keys are :class:`~.MassMetabolite` and values
            are the number currently bound to the site(s).
            An empty ``dict`` will reset the bound ligands.

        """
        return getattr(self, "_bound_metabolites")

    @bound_metabolites.setter
    def bound_metabolites(self, value):
        """Set the ``dict`` of ligands bound to the enzyme site(s)."""
        if not isinstance(value, dict) and value is not None:
            raise TypeError("bound_metabolites must be a dict")

        if value:
            bound_metabolites = {}
            for met, num_bound in iteritems(value):
                if not isinstance(met, MassMetabolite) or not isinstance(
                    num_bound, integer_types
                ):
                    raise ValueError(
                        "value must be a dict where keys are "
                        "MassMetabolites and values are ints"
                    )
                if num_bound != 0:
                    bound_metabolites[met] = num_bound
            getattr(self, "_bound_metabolites").update(bound_metabolites)
        else:
            setattr(self, "_bound_metabolites", {})

    def generate_enzyme_module_form_name(self, update_enzyme=False):
        """Generate name for the enzyme module form based on bound ligands.

        Notes
        -----
        * The :attr:`~.EnzymeModuleForm.bound_metabolites` attribute is used
          in generating the name.
        * If the :attr:`~EnzymeModuleForm.enzyme_module_id` attributes are
          not set, the string ``'Enzyme'`` will be used in its place.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`name` attribute of the enzyme
            module form in addition to returning the generated name.
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
        bound_str = "-".join(
            [met._remove_compartment_from_id_str() for met in self.bound_metabolites]
        )
        if bound_str:
            bound_str = "-" + bound_str + " complex"
        name += bound_str

        if update_enzyme:
            self.name = name

        return name

    def generate_form_formula(self, update_enzyme=False):
        """Generate the chemical formula for the enzyme module form.

        This function is primarily utilized for keeping reactions between
        :class:`EnzymeModuleForm` mass and charge balanced.

        Notes
        -----
        The :attr:`~.EnzymeModuleForm.bound_metabolites` attribute is used
        in generating the formula.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`formula` attribute of the enzyme
            module form in addition to returning the generated formula.
            Default is ``False``.

        Returns
        -------
        str
            String representing the formula of the
            :class:`EnzymeModuleForm`.

        """
        formula = ""
        if self.enzyme_module_id:
            moiety = self.enzyme_module_id.upper()
            if not re.match("^[A-Z]+[a-z]+$", moiety):
                moiety = re.sub("[0-9]", "", moiety)

            formula += "[" + moiety + "]"

        total_elements = defaultdict(list)
        if self.bound_metabolites:
            elem_iters = [
                iteritems({k: v * num_bound for k, v in iteritems(met.elements)})
                for met, num_bound in iteritems(self.bound_metabolites)
            ]

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

    def generate_form_charge(self, update_enzyme=False):
        """Generate the charge for the enzyme module form.

        This function is primarily utilized for keeping reactions between
        :class:`EnzymeModuleForm` mass and charge balanced.

        Notes
        -----
        The :attr:`~.EnzymeModuleForm.bound_metabolites` attribute is used
        in generating the charge.

        Parameters
        ----------
        update_enzyme : bool
            If ``True``, update the :attr:`charge` attribute of the enzyme
            module form in addition to returning the generated charge.
            Default is ``False``.

        Returns
        -------
        float
            Value representing the charge of the :class:`EnzymeModuleForm`.

        """
        charge = 0
        if self.bound_metabolites:
            for met, num_bound in iteritems(self.bound_metabolites):
                if met.charge is not None:
                    charge += met.charge * num_bound

        if update_enzyme:
            self.charge = charge

        return charge

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
            try:
                self._bound_metabolites = {
                    self.model.metabolites.get_by_id(str(met)): num
                    for met, num in iteritems(self.bound_metabolites)
                }
            except KeyError as e:
                raise KeyError(
                    "'{0}' does not exist in model metabolites.".format(str(e))
                )

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
                <td><strong>Bound Metabolites</strong></td>
                <td>{bound}</td>
            </tr><tr>
                <td><strong>Initial Condition</strong></td>
                <td>{ic}</td>
            </tr><tr>
                <td><strong>In {n_reactions} reaction(s)</strong></td>
                <td>{reactions}</td>
            </tr>
        </table>""".format(
            id=self.id,
            name=self.name,
            address="0x0%x" % id(self),
            enzyme=str(self.enzyme_module_id),
            compartment=self.compartment,
            bound=_make_bound_attr_str_repr(self.bound_metabolites),
            ic=self._initial_condition,
            n_reactions=len(self.reactions),
            reactions=", ".join(r.id for r in self.reactions),
        )


def _make_bound_attr_str_repr(attribute_dict):
    """Make the string representation for bound ligands.

    Warnings
    --------
    This method is intended for internal use only.

    """
    return "; ".join(
        ["{0} {1}".format(v, k) for k, v in iteritems(attribute_dict) if v != 0]
    )


__all__ = ("EnzymeModuleForm",)
