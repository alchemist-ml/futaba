from enum import Enum
import csv
import re
from chemlib import Compound


class ElGroup(Enum):
    ALKALI_METALS = 0
    ALKALINE_EARTH_METALS = 1
    TRANSITION_METALS = 2
    BORON_GROUP = 3
    CARBON_GROUP = 4
    NITROGEN_GROUP = 5
    OXYGEN_GROUP = 6
    HALOGENS = 7
    NOBLE_GASES = 8


class ComType(Enum):
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
    GENERIC = 20


class Element:
    def __init__(
        self,
        symbol: str,
        name: str,
        group: ElGroup,
        at_num: int,
        oxi: list = [0],
        valency: int = 0,
        is_metal: int = 0,
        e_neg: int = 0,
    ):
        self.symbol = symbol
        self.name = name
        self.group = group
        self.atomic_number = at_num
        self.oxi_states = oxi
        self.valency = valency
        self.is_metal = is_metal
        self.electronegativity = e_neg


def load_ions(filename: str):
    ions_dict = {}

    with open(filename, "r") as f:
        reader = csv.reader(f)
        next(reader)

        for row in reader:
            if row[0].startswith('#---'):
                continue

            try:
                ion = row[0].strip().strip('"')
                charge = int(row[1])
                ions_dict[ion] = charge
            except ValueError:
                continue

    return ions_dict


ion_charge = load_ions("Ions.csv")


class Comp:
    def __init__(
        self,
        formula: str,
        name: str,
        compound_type=ComType.GENERIC,
        ions: list = [],
        count: int = 0,
        ph: float = 7,
    ):
        self.formula = formula
        self.name = name
        self.compound_type = compound_type
        self.ions = ions
        self.ion_count = count
        self.ph = ph


def extract_ions(formula: str):
    polyatomic_ions = [
        ion
        for ion in ion_charge.keys()
        if len(ion) > 1 and ion not in ["H", "C", "N", "O", "F", "Cl", "Br", "I"]
    ]
    polyatomic_ions.sort(key=len, reverse=True)

    escaped_ions = [re.escape(ion) for ion in polyatomic_ions]

    element_pattern = r"[A-Z][a-z]?"

    all_patterns = "|".join(escaped_ions + [element_pattern])

    def parse_group(group_str, multiplier=1):
        ions = []
        i = 0
        n = len(group_str)

        while i < n:
            if group_str[i] == "(":
                depth = 1
                j = i + 1
                while j < n and depth > 0:
                    if group_str[j] == "(":
                        depth += 1
                    elif group_str[j] == ")":
                        depth -= 1
                    j += 1

                inner_content = group_str[i + 1 : j - 1]

                sub_multiplier = 1
                if j < n and group_str[j].isdigit():
                    k = j
                    while k < n and group_str[k].isdigit():
                        k += 1
                    sub_multiplier = int(group_str[j:k])
                    j = k

                inner_ions = parse_group(inner_content, sub_multiplier)
                ions.extend(inner_ions)
                i = j

            else:
                match = re.match(f"({all_patterns})(\\d*)", group_str[i:])
                if match:
                    ion = match.group(1)
                    subscript = match.group(2)
                    count = int(subscript) if subscript else 1

                    for _ in range(count * multiplier):
                        ions.append(ion)

                    i += len(match.group(0))
                else:
                    i += 1

        return ions

    return parse_group(formula)


def combination(r1: Comp, r2: Comp):
    p = ""
    i1 = ion_charge[r1.ions[0]]
    i2 = ion_charge[r2.ions[0]]

    if i1 == -1 * i2:
        p = r1.ions[0] + r2.ions[0]
    else:
        p = r1.ions[0] + str(i2) + r2.ions[0] + str(i1)

    return tidy(p)


def decomposition(r: Comp): ...


def single_displacement(r1: Comp, r2: Comp):
    p1, p2 = "", ""
    r1_i1 = ion_charge[r1.ions[0]]
    r2_i1 = ion_charge[r2.ions[0]]
    r2_i2 = ion_charge[r2.ions[1]]

    if r1_i1 == -1 * r2_i2:
        p1 = r1.ions[0] + r2.ions[1]
        p2 = r2.ions[0]
    else:
        p1 = r1.ions[0] + str(r2_i2) + r2.ions[1] + str(r1_i1)
        p2 = r2.ions[0]

    return tidy(p1, p2)


def oxidation(r1: Comp):
    o = Comp("O2", "Oxygen", ComType.GENERIC, ["O--"], 1, 7)
    p = ""
    i1 = ion_charge[r1.ions[0]]
    i2 = ion_charge[o.ions[0]]
    if ion_charge[r1.ions[0]] == 2:
        p = r1.ions[0] + o.ions[0]
        return tidy(p)
    else:
        p = r1.ions[0] + str(i2) + o.ions[0] + str(i1)
        return tidy(p)


def double_displacement(r1: Comp, r2: Comp):
    p1, p2 = "", ""
    r1_i1 = ion_charge[r1.ions[0]]
    r1_i2 = ion_charge[r1.ions[1]]
    r2_i1 = ion_charge[r2.ions[0]]
    r2_i2 = ion_charge[r2.ions[1]]

    if r1_i1 == -1 * r2_i2:
        p1 = r1.ions[0] + r2.ions[1]
    else:
        p1 = r1.ions[0] + str(r2_i2) + r2.ions[1] + str(r1_i1)

    if r2_i1 == -1 * r1_i2:
        p2 = r2.ions[0] + r1.ions[1]
    else:
        p2 = r2.ions[0] + str(r1_i2) + r1.ions[1] + str(r2_i1)

    return tidy(p1, p2)


def tidy(*args: str):
    products = []
    for p in args:
        products.append(p.replace("+", "").replace("-", "").replace("1", ""))

    products = ["H2O" if x == "HOH" else x for x in products]

    return products


def shuffle(ions: list): ...


if __name__ == "__main__":

    test_formulas = [
        "H2SO4",
        "NaCl",
        "Ca(OH)2",
        "(NH4)2SO4",
        "Al2(SO4)3",
        "Fe2(SO4)3",
        "KMnO4",
        "Na2Cr2O7",
        "CH3COOH",
        "NaHCO3",
        "Mg(NO3)2",
        "Ca3(PO4)2",
        "NH4Cl",
        "NaOH",
        "CuSO4",
        "FeCl3",
        "Na2CO3",
        "CaCO3",
    ]

    for formula in test_formulas:
        ions = extract_ions(formula)
        for i in ions:
            if ions.count(i) > 1:
                ions.append(i + str(ions.count(i)))
                ions = list(set(ions))

        print(formula, ions)
