## ======== ##
## ======== ##
import cantera as ct

from src.properties.utils import (
    create_standardized_dictionary,
    normalize_dict,
)

## ======== ##
## ======== ##

_std_ea_dict = {
    "C": 1,
    "H": 1,
    "O": 1,
    "N": 1,
    "Ar": 0,
}  # ------------#

# Molecular Weights, g/mol
_element_molecular_weights = {
    "C": 12.011,
    "H": 1.00794,
    "O": 15.9994,
    "N": 14.00674,
    "Ar": 39.948,
}  # ------------#

_std_cond = {"temp": 298.15, "pressure": 100000}  # K,Pa
_std_ea_list = set(_std_ea_dict.keys())

## ======== ##
## ======== ##


class ElementComp:
    """Mass Fraction of Elementary Compositions

    Based on the ultimate analysis of CHON.
    Input is in kg/kg, therefore less one, bigger zero.

    Attributes:
        ea_mass: dictionary with the

    ToDo
        - problem when using 'dict' as basement
        - setter and getter not working with std dict class
    """

    def __init__(self, ea_mass: dict = _std_ea_dict):
        self._ea_mass = create_standardized_dictionary(ea_mass, _std_ea_list)

    @property
    def normalized_ea_mass(self):
        """ Normalized mass compsition """
        return normalize_dict(self._ea_mass)

    @property
    def normalized_ea_molar(self):
        """ Normalized molar compsition """
        ea_mass = self.normalized_ea_mass

        # init
        mean_normal_molecular_weight = 0
        ea_molar = _std_ea_dict

        for i, value in enumerate(_std_ea_list):
            mean_normal_molecular_weight += (
                ea_mass[value] / _element_molecular_weights[value]
            )

        mean_normal_molecular_weight = 1 / mean_normal_molecular_weight

        for i, value in enumerate(_std_ea_list):
            ea_molar[value] = (
                ea_mass[value]
                / _element_molecular_weights[value]
                * mean_normal_molecular_weight
            )

        return normalize_dict(ea_molar)


## ======== ##
## ======== ##


class ElementProp:
    """Base class for all properties

    General properties are:
    Attributes:
        _density: true density, kg/m3
    """

    def __init__(self, density: float = 1000.0):
        """ """
        self._density = density

    def density(self):
        return self._density


## ======== ##
## ======== ##


class GasProp(ElementProp, ElementComp):
    """Gas properties based on cantera

    Attributes:
        species(dict): can contains any species which are
                       part of GRI mechanism
    """

    def __init__(
        self,
        species={"CH4": 1},
    ):
        # -- Init Properties
        self._species = species
        self._gas = ct.Solution("gri30.yaml")
        self._gas.Y = species
        self._ea_mass = {}

        for i, value in enumerate(_std_ea_list):
            self._ea_mass[value] = self._gas.elemental_mass_fraction(value)

        # -- INIT BASE CLASS
        ElementComp.__init__(self, ea_mass=self._ea_mass)

    def density(
        self,
        p: float = _std_cond["pressure"],
        T: float = _std_cond["temp"],
    ) -> float:
        """Density function for gases

        Args:
            p: pressure, Pa
            T: temperature, K
        """
        self._gas.TPY = T, p, self._species
        return self._gas.density


class LiquidProp(ElementProp, ElementComp):
    """Liquid Property class

    Attributes:
        ea_mass: elementar analysis, kg/kg
        crit_values: critical density, kg/m3, critical temperature, K
        density_coef: four coefficents (A,B,C,D), see VDI Wärmeatlas
    """

    def __init__(
        self,
        ea={"O": 0.888093, "H": 0.111907},  # H2O
        crit_values={"rho_crit": 322.0, "T_crit": 647.10},  # H2O
        density_coef=[1094.0233, -1813.2295, 3863.9557, -2479.8130],  # H2O
    ):
        # -- INIT BASE CLASS
        ElementComp.__init__(self, ea_mass=ea)
        self.crit_values: dict = crit_values
        self.density_coef: list = density_coef

    def density(
        self,
        T: float = 293.15,
    ) -> float:
        """Density of common liquids

        Args:
            T: temperature, K
        """
        rho0 = self.crit_values["rho_crit"]
        T0 = self.crit_values["T_crit"]

        den = rho0 + (
              self.density_coef[0] * (1 - T / T0) ** (0.35)
            + self.density_coef[1] * (1 - T / T0) ** (2 / 3)
            + self.density_coef[2] * (1 - T / T0) ** (1)
            + self.density_coef[3] * (1 - T / T0) ** (4 / 3)
        )

        return den

## ------- ##
## ------- ##


class SolidProp(ElementProp, ElementComp):
    """Solid Base Property class for Hydrocarbons
    
    Attributes:
        ea: Elementar Analysis [CHON], kg/kg
        density: constant density, kg/m3
    """

    def __init__(
        self,
        ea=_std_ea_dict,
        density: float = 1000,
    ):
        # -- INIT BASE CLASS
        ElementProp.__init__(self, density=density)
        ElementComp.__init__(self, ea_mass=ea)

