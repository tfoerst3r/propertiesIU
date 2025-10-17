import cantera as ct
import numpy as np

## ======== ##
## ======== ##

std_ea_dict = {"C": 1, "H": 0, "O": 0, "N": 0}
std_ea_list = list(std_ea_dict.keys())

## ======== ##
## ======== ##


def create_standardized_dictionary(input_dict, default_list):
    """
    Creates a new dictionary based on a standard list.

    - Elements in the standard list but not in the input_dict are set to 0.
    - If the input_dict contains elements not in the standard list, it raises an error.

    Args:
        input_dict (dict): The input dictionary with some keys and values.
        default_list (list): The standard list of keys.

    Returns:
        dict: The new standardized dictionary.

    Raises:
        ValueError: If the input_dict contains a key not present in the default_list.
    """

    # 1. Error check for extra elements in input_dict
    for key in input_dict:
        if key not in default_list:
            raise ValueError(
                f"Error: Key '{key}' in input dictionary is not present in the standard list."
            )

    # 2. Create the new standardized dictionary
    standardized_dict = {}
    for key in default_list:
        # Get the value from input_dict if the key exists, otherwise use 0
        standardized_dict[key] = input_dict.get(key, 0)

    # 3. Error handling; if all values are zero the first will be set to 1
    checksum = 0
    first_key = list(standardized_dict.keys())[0]

    for i, val in enumerate(standardized_dict):
        checksum += standardized_dict[val]

    if checksum == 0:
        standardized_dict[first_key] = 1

    return standardized_dict


## ======== ##
## ======== ##


def normalize_dict(inquiry: dict) -> dict:
    """Normalize elements of a dictionary"""
    if not inquiry:
        print(">>> Warning! Dictionary was empty! <<<")
        return {}

    inquiry_new = {}

    try:
        for i, value in enumerate(dict(inquiry)):
            inquiry_new[value] = inquiry[value]
    except TypeError as te:
        print(f"Not a dictionary type class! {te}")

    for i, (valA, valB) in enumerate(inquiry_new.items()):
        if type(valB) is not float and type(valB) is not int:
            raise ValueError(valB, "is not a number")
        if float(valB) < 0:
            print("Input is negative. They are ignored!")
            continue

    sum = 0
    for i, (valA, valB) in enumerate(inquiry_new.items()):
        if valB < 0:
            valB = 0
        sum += valB

    for i, (valA, valB) in enumerate(inquiry_new.items()):
        inquiry_new[valA] = valB / sum

    return inquiry_new


## ======== ##
## ======== ##


def normalize_array(inquiry: list) -> list:
    """Normalize a list
    Returns:
        A normalized list
    """
    sum = 0

    for i, val in enumerate(inquiry):
        sum += val

    for i, val in enumerate(inquiry):
        inquiry[i] = inquiry[i] / sum

    return inquiry


## ======== ##
## ======== ##


def test_type(inquiry: object, type_def: type) -> None:
    """Type testing routine"""
    if not isinstance(inquiry, type_def):
        error_string = f"Not a {type_def.__name__}!"
        raise TypeError(error_string)


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

    def __init__(self, ea_mass: dict = std_ea_dict):  # --------#
        self._ea_list = ["C", "H", "O", "N"]
        self._molecular_weights = {
            "C": 12.011,
            "H": 1.00794,
            "O": 15.9994,
            "N": 14.00674,
        }

        self._ea_mass = create_standardized_dictionary(ea_mass, std_ea_list)

    @property
    def normalized_ea_mass(self):
        return normalize_dict(self._ea_mass)

    @normalized_ea_mass.setter
    def normalized_ea_mass(self, dict_ea_mass: dict):
        """Setter: Elemental Analysis, Mass based"""
        if not list(set(self._ea_list).difference(dict_ea_mass)):
            self._ea_mass = dict_ea_mass
        else:
            diff_ea_mass = set(self._ea_list).difference(dict_ea_mass)
            raise ValueError(
                "<{key}> is/are not part of possible input: {eee}.".format(
                    key=diff_ea_mass, eee=self._ea_list
                )
            )

    @property
    def normalized_ea_molar(self):
        mean_normal_molecular_weight = 0
        ea_mass = self.normalized_ea_mass
        ea_molar = {"C": 0, "H": 0, "O": 0, "N": 0}

        for i, value in enumerate(std_ea_list):
            mean_normal_molecular_weight += (
                ea_mass[value] / self._molecular_weights[value]
            )

        mean_normal_molecular_weight = 1 / mean_normal_molecular_weight

        for i, value in enumerate(std_ea_list):
            ea_molar[value] = (
                ea_mass[value]
                / self._molecular_weights[value]
                * mean_normal_molecular_weight
            )

        return normalize_dict(ea_molar)


## ======== ##
## ======== ##


class ElementProp:
    """Base class for all properties
    **Introduction**

    General properties are:

    - Molar weight
    - density, true density

    Attributes:
        mw: molar weight, kg/kmol

    Todo:
        - adding more properties to general
    """

    def __init__(self, mw: float = 1.0):
        """ """
        self._mole_weight = mw
        self._density = 0.0
        self._std_formation_mass = 0.0
        self._std_formation_mole = 0.0

    @property
    def std_formation_mole(self):
        """Standard molar heat of formation kJ/kmol"""
        try:
            1 / self._mole_weight
        except ValueError:
            print(f"MW not set propertly: {self._mole_weight}")
            raise

        self._std_formation_mole = self.std_formation_mass * self.mole_weight
        return self._std_formation_mole

    @property
    def std_formation_mass(self):
        """Standard mass heat of formation kJ/kg"""
        return self._std_formation_mass

    @std_formation_mass.setter
    def std_formation_mass(self, value):
        """Setter for Standard mass heat of formation kJ/kg"""
        self._std_formation_mass = value

    @property
    def mole_weight(self):
        """Molar weight kg/kmol"""
        return self._mole_weight

    @mole_weight.setter
    def mole_weight(self, value):
        self._mole_weight = value

    @property
    def density(self):
        return self._density

    @density.setter
    def density(self, value):
        self._density = value


## ======== ##
## ======== ##


class GasProp(ElementProp, ElementComp):
    """Gas properties based on cantera

    Gas Class consist of several species properties's

    Attributes:
        species(dict): can containes species which are
                       part of GRI mechanism
    Todo:
        - Replace ElementComp with cantera routines
    @details ----
    """

    def __init__(
        self,
        species={"CH4": 1},
    ):
        # -- INIT BASE CLASS
        ElementProp.__init__(self)
        ElementComp.__init__(self)

        # -- Other properties
        self._gas = ct.Solution("gri30.yaml")
        self._mixture_dh_formation = 0
        self._std_cond = [298.15, 100000]  # K,Pa
        self._ea_mass = self.normalized_ea_mass
        self._ea_molar = self.normalized_ea_molar
        self.__set_prop("mass", species_dict=species)

    @property
    def std_formation_mass(self):
        """Standard Heat of Formation
        return:
            dh0 J/kg
        """
        return self._mixture_dh_formation

    @std_formation_mass.setter
    def std_formation_mass(self):
        raise AttributeError("Setter disabled. Set via gas species composition.")

    @property
    def mole_weight(self):
        """
        @return: kg/mol
        """
        return self._mole_weight

    @property
    def species_mass(self):
        """
        Returns:
            dictionary with species and mass fraction
        """
        return self._species_mass

    @species_mass.setter
    def species_mass(self, species_dict):
        self.__set_prop("mass", species_dict=species_dict)

    @property
    def species_molar(self):
        """
        Returns:
            dictionary with species and mole fraction
        """
        return self._species_molar

    @species_molar.setter
    def species_molar(self, species_dict):
        self.__set_prop("molar", species_dict=species_dict)

    def __set_prop(self, input_type, species_dict):
        """Sets several class properties attributes
        Returns:
            - species_mass
            - species_molar
            - ea_mass
            - ea_molar
            - mixture_dh_formation
            - molar_weight
        """

        try:
            self._gas.TPX = 298.15, 100000, species_dict
        except AttributeError as ae:
            print(f"Wrong input: {ae}")
            raise

        if input_type == "molar":
            self._gas.TPX = self._std_cond[0], self._std_cond[1], species_dict
        elif input_type == "mass":
            self._gas.TPY = self._std_cond[0], self._std_cond[1], species_dict

        self._species_mass = {}
        self._species_molar = {}
        for i, value in enumerate(species_dict):
            index = self._gas.species_index(value)

            self._species_molar[value] = self._gas.X[index]
            self._species_mass[value] = self._gas.Y[index]

        for i, var in enumerate(self._ea_mass):
            self._ea_mass[var] = self._gas.elemental_mass_fraction(var)
            self._ea_molar[var] = self._gas.elemental_mole_fraction(var)

        self._mixture_dh_formation = self._gas.h
        self._mole_weight = self._gas.mean_molecular_weight / 1000


## ======== ##
## ======== ##


class LiquidProp(ElementProp, ElementComp):
    """
    @brief
    @details ----
    Liquid Property class
    Attributes:
        ua: ultimate analysis kg/kg
        mw: molar weight
    Todo:
        - Include species definition as well
    @details ----
    """

    def __init__(self, ua={"H": 0.5, "C": 0.5}, mw=1):
        # -- INIT BASE CLASS
        ElementProp.__init__(self, mw=mw)
        ElementComp.__init__(self, ea_mass=ua)

        self._molar_composition = {}

    def molar_composition(self):
        tar_ea_mass = self.normalized_ea_mass
        mw_elements = standard_values(tar_ea_mass)[2]
        self._molar_composition = {}
        for i, value in enumerate(tar_ea_mass):
            self._molar_composition[value] = (
                tar_ea_mass[value] * self.mole_weight / mw_elements[value]
            )

        return self._molar_composition


## ------- ##
## ------- ##


class SolidProp(ElementProp, ElementComp):
    """Solid Base Property class for Hydrocarbons

    **Materials**

    - Reactive materials containing hydrocarbons
    - Coal
    - Plastics
    - Droplets (technical a liquid)

    Containing properties:

    - Ultimate Analysis [CHON], weight fraction
    - Proximate Analysis [FC Vol Ash], weight fraction daf
    - Higher Heating Value
    - Porosity

    """

    def __init__(
        self,
        ua={"C": 1.0, "H": 0.0, "O": 0.0, "N": 0.0},
        pa={"FC": 1.0, "Vol": 0.0, "Ash": 0.0},
        hhv=0.0,
        porosity=0.0,
        water_mass_fraction=0.0,
        density_bulk=0.0,
        ptemp=298.15,
    ):
        # -- INIT BASE CLASS
        ElementProp.__init__(self)
        ElementComp.__init__(self, ea_mass=ua)

        # -- Setup other properties
        self._pa_list = ["FC", "Vol", "Ash"]

        if not list(set(pa).difference(self._pa_list)):
            self._proximate = normalize_dict(pa)
        else:
            diff_pa = set(pa).difference(self._pa_list)
            raise ValueError(
                "<{key}> is/are not part of possible input: {eee}.".format(
                    key=diff_pa, eee=self._pa_list
                )
            )

        self._hhv = hhv
        self._porosity = porosity
        self._nmr = self._nmr_correlation()
        self._density_bulk = density_bulk
        self._water_mass_fraction = water_mass_fraction
        self._water_volume_fraction = 0
        self._ptemp = ptemp

    @property
    def ptemp(self):
        return self._ptemp

    @ptemp.setter
    def ptemp(self, value):
        self._ptemp = value

    @property
    def water_mass_fraction(self):
        return self._water_mass_fraction

    @water_mass_fraction.setter
    def water_mass_fraction(self, value):
        if isinstance(value, int) or isinstance(value, float):
            self._water_mass_fraction = value
        else:
            raise AttributeError("Only float or integer values allowed!")

    # @property
    # def water_volume_fraction(self):
    #     """Based on the mass fraction of water
    #     """
    #     if self.water_mass_fraction == 0 or self.density_bulk == 0:
    #         return 'Water mass fraction or bulk density not updated!'
    #     else:
    #         density_water = 998.4 # av. @ 298 K
    #         print(density_water)

    #         fraction = ( 1 + ( 1/self.water_mass_fraction -1 ) * density_water/self.density_bulk )**(-1)
    #         return fraction

    # @water_volume_fraction.setter
    @property
    def water_volume_fraction(self):
        """Based on the mass fraction of water"""
        ptemp = self._ptemp

        def density(temp):
            """Calculate water density, based on temperature

            The calculation is based on a simple regression polynom 4th grade
            Arguments:
                temp: temperature, K
            Returns:
                water density, kg/m3
            """
            const = [104.767, 9.09422, -0.0327719, 5.0689e-05, -3.14159e-08]

            return (
                const[0]
                + const[1] * temp
                + const[2] * temp**2
                + const[3] * temp**3
                + const[4] * temp**4
            )

        if self.water_mass_fraction == 0 or self.density_bulk == 0:
            return "Water mass fraction or bulk density not updated!"
        else:
            # this needs some update!
            # density_water = 998.4 #kg/m3 @ 298 K
            density_water = density(ptemp)  # av. @ 343 K

            print("Init density of H2O,kg/m3:", density_water, "at", ptemp, "K.")
            fraction = (
                1
                + (1 / self.water_mass_fraction - 1) * density_water / self.density_bulk
            ) ** (-1)
            return fraction

    @property
    def density_bulk(self):
        return self._density_bulk

    @density_bulk.setter
    def density_bulk(self, value):
        self._density_bulk = value

    @property
    def nmr(self):
        self._nmr = self._nmr_correlation()
        return self._nmr

    @property
    def proximate(self):
        raise AttributeError(
            "Only a full input of: {eee} is supported!".format(eee=self._pa_list)
        )

    @proximate.setter
    def proximate(self, dict_pa: dict):
        if not list(set(dict_pa.keys()).difference(self._pa_list)):
            self._proximate = dict_pa
        else:
            diff_pa = set(dict_pa).difference(self._pa_list)
            raise ValueError(
                "<{key}> is/are not part of possible input: {eee}.".format(
                    key=diff_pa, eee=self._pa_list
                )
            )

    def normalized_proximate(self):
        return normalize_dict(self._proximate)

    @property
    def hhv(self):
        """Higher Heating Value dry ash free!"""
        return self._hhv

    @hhv.setter
    def hhv(self, value):
        print("Setting value")
        self._hhv = value

    def _nmr_correlation(self):
        """Calc parameters using Genetti correlation

        **The Correlation**

        Source
        : https://www.et.byu.edu/~tom/cpd/correlation.html

        Equation to calculate M del, MW, Po and sigma+1:

        Y = c1+c2*C+c3*C^2+c4*H+c5*H^2+c6*O+c7*O^2+c8*VM+c9*VM^2

            where:
            C = Weight Percent Carbon (daf)
            H = Weight Percent Hydrogen (daf)
            N = Weight Percent Nitrogen (daf)
            O = Weight Percent Oxygen (daf)
            VM = ASTM Volatile Matter (daf)

        Returns:
            nmr: NMR parameter as a list
        """

        cpdv = np.array(
            [
                [0.000, 0.000, 0.000, 0.000],
                [421.957, 1301.41, 0.489809, -52.1054],
                [-8.64692, 16.3879, -0.00981566, 1.63872],
                [0.0463894, -0.187493, 0.000133046, -0.0107548],
                [-8.47272, -454.773, 0.155483, -1.23688],
                [1.18173, 51.7109, -0.0243873, 0.0931937],
                [1.15366, -10.0720, 0.00705248, -0.165673],
                [-0.0434024, 0.0760827, 0.000219163, 0.00409556],
                [0.556772, 1.36022, -0.0110498, 0.00926097],
                [-0.00654575, -0.0313561, 0.000100939, -0.0000826717],
            ]
        )

        nmr = {"mdel": 0, "mw": 0, "p0": 0, "sig": 0, "c0": 0}
        # ultimate = [ fcar,fhyd,foxy,fn,fs ]
        ulti = [0.0, 0.0, 0.0, 0.0, 0.0]
        # proxy = [ fc,vm,ash ]
        proxy = [0.45, 0.45, 0.1]

        for i, value in enumerate(self._pa_list):
            proxy[i] = self.normalized_proximate()[value]

        proxy[2] = 0
        proxy = normalize_array(proxy)

        for i, value in enumerate(self.normalized_ea_mass):
            ulti[i] = self.normalized_ea_mass[value]

        ulti = normalize_array(ulti)

        c = np.zeros(10)

        for i, valA in enumerate(nmr):
            if valA == "c0":
                nmr[valA] = min(0.36, max((0.118 * ulti[0] * 100 - 10.1), 0.0)) + min(
                    0.15, max((0.014 * ulti[2] * 100 - 0.175), 0.0)
                )
            else:
                for j, valB in enumerate(cpdv):
                    c[j] = cpdv[j][i]

                nmr[valA] = (
                    c[1]
                    + c[2] * (ulti[0] * 100.0)
                    + c[3] * (ulti[0] * 100) ** 2
                    + c[4] * (ulti[1] * 100)
                    + c[5] * (ulti[1] * 100) ** 2
                    + c[6] * (ulti[2] * 100)
                    + c[7] * (ulti[2] * 100) ** 2
                    + c[8] * (proxy[1] * 100)
                    + c[9] * (proxy[1] * 100) ** 2
                )
        return nmr


def standard_values(species=["CH4"]):
    """
    Standard enthalpies and MW are determined!

    @param species: both dictionary or list can be used
    :return:
    """
    std_enthalpy_formation_mass = {}
    std_enthalpy_formation_molar = {}
    mole_weights_species = {}
    gas = ct.Solution("gri30.yaml")
    gas.TP = 298.15, 100000  # K,Pa, std conditions

    if isinstance(species, dict):
        species_list = species.keys()
    elif isinstance(species, list):
        species_list = species

    for i, value in enumerate(species_list):
        try:
            mole_weights_species[value] = (
                gas.molecular_weights[gas.species_index(value)] / 1000
            )  # kg/mol

            # TODO 'C' and 'C<s>' need to be distingished
            if value in ["N2", "O2", "H2", "C", "H", "O", "N", "S"]:
                std_enthalpy_formation_molar[value] = 0.0

            else:
                # ct.gas_constant, J/kmol.K
                std_enthalpy_formation_molar[value] = (
                    gas.standard_enthalpies_RT[gas.species_index(value)]
                    * gas.T
                    * (ct.gas_constant / 1000)
                )

            std_enthalpy_formation_mass[value] = (
                std_enthalpy_formation_molar[value] / mole_weights_species[value]
            )
        except:
            print(">> Warning: Species <", value, "> unknown in std values list! <<")
            continue

    return (
        std_enthalpy_formation_mass,
        std_enthalpy_formation_molar,
        mole_weights_species,
    )


def diff_enthalpy_values(species={"CH4": 1}, cond=[298.15, 100000]):
    """
    Only enthalpy difference dh
    due to temperature change.
    Reference is standard condition at 298.15 K and 1e5 Pa
    :return sensible dh, J/mol
    """

    gas = ct.Solution("gri30.yaml")
    std_cond = [298.15, 100000]  # K,Pa, std conditions
    dh = {}

    for i, value in enumerate(species):
        dh1 = 0
        dh2 = 0
        try:
            gas.TPY = std_cond[0], std_cond[1], {value: 1}
            dh1 = gas.h
            gas.TPY = cond[0], cond[1], {value: 1}
            dh2 = gas.h
            dh[value] = dh2 - dh1
        except:
            print(">> Warning: Species", value, "unknown! <<")
            continue
    print("-----------------------------")

    return dh


def h2o_saturation(pressure):
    """
    @detail ----
    latent heat of water evaporation
    Not covered here, but this is the motivation:
    This fitting routine is error prone!
    T_bub = 298.15 K
    p_bub = 3170 Pa
    v = 0.00000 m^3/kg
    h = 0.00 kJ/kg
    s = 0.00 J/(kg K)
    h_liq,gas = 2442.71781908 kJ/kg
    h_liq,gas = ~44000 J/mol
    @detail ----

    @param pressure:
    @return: (dh_mass,dh_mole,temp_saturation)
    """

    def fct_latent(coeff, pressure):
        return (
            coeff[0]
            + coeff[1] * pressure
            + coeff[2] * pressure**2
            + coeff[3] * pressure ** coeff[4]
        )

    # 1000 < P < 4.74e5
    lateten1 = np.array(
        [
            -4806214.25985666,
            -0.350219461494522,
            2.65331149894677e-07,
            7624727.13863795,
            -0.00622219768595071,
        ]
    )

    vapT1 = np.array(
        [
            175.609422862424,
            1.49193767247892e-05,
            -6.47740090439212e-12,
            41.0209084682707,
            0.135722389565108,
        ]
    )

    # 4.74e5 < P < 110e5
    lateten2 = np.array(
        [
            -8795956.80499871,
            -0.0514944016026154,
            4.53629713161923e-11,
            12406593.2303227,
            -0.00963508724003428,
        ]
    )

    vapT2 = np.array(
        [
            234.709975060412,
            1.05201444912404e-06,
            -1.79728263568382e-14,
            14.6133399245996,
            0.195378177282553,
        ]
    )

    if pressure < 1000 or pressure > 110e5:
        print("\n>> Warning! <<")
        print("Function routine was fit for 1e3 Pa < pressure < 1.1e7 Pa")

    # Own fit, dh in J/mol
    if pressure < 4.74e5:
        # dh_molar = 42957.9 - 2511.02*P**0.5 + 136.86*P -1.05712*P**2
        dh_mass = fct_latent(lateten1, pressure)
        vapT = fct_latent(vapT1, pressure)
    else:
        # dh_molar = 35449.7 - 117.716*P
        dh_mass = fct_latent(lateten2, pressure)
        vapT = fct_latent(vapT2, pressure)
        # J/kg * kg/mol => J/mol
    dh_molar = dh_mass * 0.01801528

    return dh_mass, dh_molar, vapT


def composition_matrix(species={"CH4": 1}):
    """Generation of the CHONS composition matrix
    @details -----

    A matrix with CHONS composition is generated
    __Example__

    Elemental mass fraction of O2 CO2 H2O N2 SO2
    [ 0     , 12/44 , 0     , 0     , 0     ] [C]
    [ 0     , 0     , 2/18  , 0     , 0     ] [H]
    [ 32/32 , 32/44 , 16/18 , 0     , 32/64 ] [O]
    [ 0     , 0     , 0     , 28/28 , 0     ] [N]
    [ 0     , 0     , 0     , 0     , 32/64 ] [S]

    @details -----
    :param used_species:
    :return:
    """
    if isinstance(species, dict):
        species_list = list(species.keys())
    elif isinstance(species, list):
        species_list = species
    else:
        raise AttributeError("Unknown species input.")

    std_cond = [298.15, 100000]  # K,Pa
    elements = ["C", "H", "O", "N", "S"]
    gas = ct.Solution("gri30.yaml")

    A = np.zeros((len(elements), len(species_list)))

    # -- Generating solving matrix
    for i, value in enumerate(species_list):
        gas.TPY = std_cond[0], std_cond[1], {value: 1}

        for j, element in enumerate(elements):
            A[j][i] = gas.elemental_mass_fraction(element)

    return A


def higher_heating_value(
    fuel_ea={"C": 0.8, "H": 0.05, "O": 0.1, "N": 0.025, "S": 0.025},
    fuel_pa={"FC": 0.5, "Vol": 0.4, "Ash": 0.1},
):
    # Input need to be dry
    # Higher Heating Value defined by Mason and Gandhi (see tomeczek_1994)
    # eq.: a1 * C(wf) + a2 * H(wf) + a3 * S(wf) + a4 * Ash(wf) + a5 * (O+N)(wf)

    fuel_pa = normalize_dict(fuel_pa)
    fuel_ea = normalize_dict(fuel_ea)

    array1 = [34095, 132298, 6848, -1531, -11996]

    array2 = {}

    for i, value in enumerate(fuel_ea):
        array2[value] = fuel_ea[value] * (1 - fuel_pa["Ash"])

    array2["Ash"] = fuel_pa["Ash"]

    array3 = []
    array3.append(array2["C"])
    array3.append(array2["H"])
    array3.append(array2["S"])
    array3.append(array2["Ash"])
    array3.append(array2["O"] + array2["N"])

    hhv = 0
    for i, value in enumerate(array1):
        hhv += array1[i] * array3[i]

    # kJ/kg -> J/kg (SI)
    hhv *= 1000

    return hhv
