from enum import Enum


class ElementGroup(Enum):
    ALKALI_METALS = 0
    ALKALINE_EARTH_METALS = 1
    TRANSITION_METALS = 2
    BORON_GROUP = 3
    CARBON_GROUP = 4
    NITROGEN_GROUP = 5
    OXYGEN_GROUP = 6
    HALOGENS = 7
    NOBLE_GASES = 8


class CompoundType(Enum):
    ACID = 0
    BASE = 1
    AMPHOTERIC = 2
    SALT = 3
    HALIDE = 4
    OXYSALT = 5
    OXIDE = 6
    HYDROXIDE = 7
    PEROXIDE = 8
    SUPEROXIDE = 9
    HYDRIDE = 10
    HYDRATE = 11
    AMMONIATE = 12
    SULPHIDE = 13
    NITRIDE = 14
    CARBIDE = 15
    PHOSPHIDE = 16
    BORIDE = 17
    SILICIDE = 18
    COORDINATION = 19


class Element:
    def __init__(
        self,
        symbol: str,
        name: str,
        group: ElementGroup,
        at_num: int,
        group_num: int,
        oxi: list = None,
        valency: int = 0,
        is_metal: int = 0,
        e_neg: int = 0,
    ):
        self.symbol = symbol
        self.name = name
        self.group = group
        self.atomic_number = at_num
        self.group_number = group_num
        self.oxi_states = oxi if oxi is not None else [0]
        self.valency = valency
        self.is_metal = is_metal
        self.electronegativity = e_neg


class Ion:
    def __init__(
        self,
        fmla: str,
        name: str,
        charge: int,
        ion_type: str,
        elements: list,
        count: int,
    ):
        self.formula = fmla
        self.name = name
        self.charge = charge
        self.ion_type = ion_type
        self.elements = elements
        self.element_count = count


class Compound:
    def __init__(
        self,
        formula: str,
        name: str,
        compound_type: CompoundType,
        ions: list,
        count: int,
        ph: float,
    ):
        self.formula = formula
        self.name = name
        self.compound_type = compound_type
        self.ions = ions
        self.ion_count = count
        self.ph = ph
