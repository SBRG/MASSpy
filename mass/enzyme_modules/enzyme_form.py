# -*- coding: utf-8 -*-
"""TODO Module Docstrings."""
from __future__ import absolute_import

from collections import defaultdict
from itertools import chain

from six import integer_types, iteritems, string_types

from mass.core.massmetabolite import MassMetabolite


class EnzymeForm(MassMetabolite):
    """Class for holding information regarding an enzymatic species.

    Parameters
    ----------
    id: str, mass.MassMetabolite
        The identifier associated with the EnzymeForm, or an exisitng
        MassMetabolite object. If a MassMetabolite object is provided, an 
        EnzymeForm will be instanstiated with stored values identicial to
        the value of those stored in the MassMetabolite.
    name: str, optional
        A human readable name for the enzyme form.

    Attributes
    ----------
    enzyme_id: str, optional
        The identifier of the enzyme represented by the EnzymeForm.
    enzyme_name: str, optional
        A human readable name of the enzyme represented by the EnzymeForm.
    bound_catalytic: dict, optional
        A dict representing the metabolites bound to the enzyme's active 
        site(s), with MassMetabolites as keys and number bound as values.
        If the final coefficient for a metabolite is 0 then it is removed.
    bound_effectors: dict, optional
        A dict representing the metabolites bound to the enzyme's regulatory 
        site(s), with MassMetabolites as keys and the number bound as values.
        If the final coefficient for a metabolite is 0 then it is removed.
    formula: str, optional
        Chemical formula associated with the EnzymeForm, accounting for 
        bound metabolites.
    charge: float, optional
        The charge number associated with the EnzymeForm, accounting for 
        bound metabolites.
    compartment: str, optional
        The compartment where the EnzymeForm is located.

    """

    def __init__(self, id=None, name="", enzyme_id="", enzyme_name="",
                 bound_catalytic=None, bound_effectors=None, 
                 formula=None, charge=None, compartment=None):
        """Initialize the EnzymeForm Object."""
        # Initialize MassMetabolite parent class
        super(EnzymeForm, self).__init__(id, name, compartment=compartment)
        # Set the id of the enzyme represented by the EnzymeForm
        for attr, value in zip(["enzyme_id", "enzyme_name"], 
                               [enzyme_id, enzyme_name]):
            if not isinstance(attr, string_types):
                raise TypeError(attr + " must be a str")
            else:
                setattr(self, attr, value)

        # Set metabolites bound to the active site(s) of the enzyme form
        self._bound_catalytic = {}
        self.bound_catalytic = bound_catalytic

        # Set metabolites bound to the regulatory site(s) of the enzyme form
        self._bound_effectors = {}
        self.bound_effectors = bound_effectors

        # Set formula, charge, and compartment attributes
        for attr, value in zip(["formula", "charge"], [formula, charge]):
            if value is None:
                value = self.__class__.__dict__["get_species_" + attr](self)
            setattr(self, attr, value)
        self.compartment = compartment

    @property
    def bound_catalytic(self):
        """Return metabolites bound to the active site(s) of the enzyme."""
        return getattr(self, "_bound_catalytic")

    @bound_catalytic.setter
    def bound_catalytic(self, value):
        """Set the dictionary of metabolites bound to the active site(s).

        Assigning a dict to this property updates the EnzymeForms's dict of 
        metabolites bound at the active site(s) with the new values.

        Parameters
        ----------
        value : dict
            Dict where keys are MassMetabolites representing metabolites and
            values are the number currently bound to the active site(s). 
            An empty dict will reset the bound metabolites.

        """
        self._set_bound_dict("bound_catalytic", value)

    @property
    def bound_effectors(self):
        """Return metabolites bound to the regulatory site(s) of the enzyme."""
        return getattr(self, "_bound_effectors")

    @bound_effectors.setter
    def bound_effectors(self, value):
        """Set the dictionary of metabolites bound to the regualtory site(s).

        Assigning a dict to this property updates the EnzymeForms's dict of 
        metabolites bound at the regulatory site(s) with the new values.

        Parameters
        ----------
        value : dict
            Dict where keys are MassMetabolites representing metabolites and 
            values are the number currently bound to the regulatory site(s). 
            An empty dict will reset the bound metabolites.

        """
        self._set_bound_dict("bound_effectors", value)

    def generate_enzyme_form_name(self, use_enzyme_name=False, 
                                  update_enzyme=False):
        """Generate a name for the enzymatic species based on bound ligands.

        The enzyme_name (if already set) or enzyme_id, bound_catalytic,
        and bound_effector attributes are used in generating the name of the
        EnzymeForms.

        Parameters
        ----------
        enzyme_name: bool, optional
            If True, then use the enzyme_name attribute in the generation of 
            the EnzymeForm name. Otherwise use the enzyme_id or str "Enzyme"
            if no enzyme_id is defined.
        update_enzyme: bool, optional
            If True, update the name attribute of the enzyme form in
            addition to returning the automatically generated name of the 
            enzyme form as a str. Default is False.

        Returns
        -------
        name: A str representing the name of the enzymatic species.

        Notes
        -----
        If the enzyme_name and the enzyme_id attributes are not set, the string
        "Enzyme" will be used in its place.

        """
        if self.enzyme_name and use_enzyme_name:
            name = self.enzyme_name
        elif self.enzyme_id:
            name = self.enzyme_id
        else:
            name = "Enzyme"

        # Add the ligands bound to the active site(s) 
        catalytic_str = "-".join([met._remove_compartment_from_id_str()
                                  for met in self.bound_catalytic])
        if catalytic_str:
            catalytic_str = "-" + catalytic_str + " complex"
        name += catalytic_str

        # Add the ligands bound to the effector site(s)
        effector_str = self._make_bound_ligand_str_repr(self.bound_effectors)
        if effector_str:
            name += "; " + effector_str

        if update_enzyme:
            self.name = name

        return name

    def get_species_formula(self, update_enzyme=False):
        """Get the chemical formula for the enzyme form and bound species.

        Parameters
        ----------
        update_enzyme: bool, optional
            If True, update the formula attribute of the enzyme form in
            addition to returning the automatically generated name of the 
            enzyme form as a str. Default is False.

        """
        formula = ""
        if self.enzyme_id:
            formula += "[" + str(self.enzyme_id) + "]"

        total_elements = defaultdict(list)
        for dictionary in [self.bound_catalytic, self.bound_effectors]:
            if not dictionary:
                continue
            elem_iters = [iteritems({
                k: v * num_bound for k, v in iteritems(met.elements)}) 
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
        """Get the chemical formula for the enzyme form and bound species.

        Parameters
        ----------
        update_enzyme: bool, optional
            If True, update the charge attribute of the enzyme form in
            addition to returning the automatically generated name of the 
            enzyme form as a str. Default is False.

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

    def _make_bound_ligand_str_repr(self, attribute_dict):
        """Make the string representation for bound ligands.

        Warnings
        --------
        This method is intended for internal use only.

        """
        return "; ".join([
            "{0:d} {1}".format(v, k._remove_compartment_from_id_str())
            for k, v in iteritems(attribute_dict) if v != 0])

    def _set_id_with_model(self, value):
        """Set the id of the EnzymeForm to the associated MassModel.

        Warnings
        --------
        This method is intended for internal use only.

        """
        super(EnzymeForm, self)._set_id_with_model(value)
        self.model.enzyme_forms._generate_index()

    def _repair_bound_pointers(self):
        """Repair object pointer for metabolites in bound dict attributes.

        Requires a model to be associated with the EnzymeForm.

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
        """HTML representation of the overview for the EnzymeForm."""
        return """
        <table>
            <tr>
                <td><strong>EnzymeForm identifier</strong></td>
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
            enzyme=str(self.enzyme_id), compartment=self.compartment,
            catalytic=self._make_bound_ligand_str_repr(self.bound_catalytic),
            effectors=self._make_bound_ligand_str_repr(self.bound_effectors),
            ic=self._initial_condition, n_reactions=len(self.reactions),
            reactions=', '.join(r.id for r in self.reactions))
