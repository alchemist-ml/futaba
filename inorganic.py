from enum import Enum
import csv


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
    ELEMENT = 20


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
            if row[0].startswith("#"):
                continue

            try:
                ion = row[0].strip().strip('"')
                charge = int(row[1])
                ions_dict[ion] = charge
            except ValueError:
                continue

    return ions_dict


ion_charge = load_ions("Ions.csv")


# Conjugate acids for common neutral bases
Conjugate_Acid = {
    "NH3": "NH4+",
}


class Compound:
    def __init__(
        self,
        formula: str,
        name: str,
        compound_type: ComType,
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


def extract_ions(c: Compound): ...


def combination(r1: Compound, r2: Compound):
    p = ""
    i1 = ion_charge[r1.ions[0]]
    i2 = ion_charge[r2.ions[0]]

    if i1 == -1 * i2:
        p = r1.ions[0] + r2.ions[0]
    else:
        p = r1.ions[0] + str(i2) + r2.ions[0] + str(i1)

    return tidy(p)


def decomposition(r: Compound): ...


def single_displacement(r1: Compound, r2: Compound):
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


def oxidation(r1: Compound):
    o = Compound("O2", "Oxygen", ComType.ELEMENT, ["O--"], 1, 7)
    p = ""
    i1 = ion_charge[r1.ions[0]]
    i2 = ion_charge[o.ions[0]]
    if ion_charge[r1.ions[0]] == 2:
        p = r1.ions[0] + o.ions[0]
        return tidy(p)
    else:
        p = r1.ions[0] + str(i2) + o.ions[0] + str(i1)
        return tidy(p)


def double_displacement(r1: Compound, r2: Compound):
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


def acid_base_test():

    r1 = Compound("HCl", "Sulphuric Acid", ComType.ACID, ["H+", "Cl-"], 2, 1)
    r2 = Compound("NaOH", "Sodium Chloride", ComType.BASE, ["Na+", "OH-"], 2, 14)

    if r1.compound_type == ComType.ACID and r2.compound_type == ComType.BASE:
        products = double_displacement(r1, r2)
        print(f"{r1.formula} + {r2.formula}")
        print(products)


def combination_test():
    r1 = Compound("Mg", "Magnesium", ComType.ELEMENT, ["Mg++"], 1, 7)
    r2 = Compound("O2", "Oxygen", ComType.ELEMENT, ["O--"], 1, 7)

    products = combination(r1, r2)
    print(f"{r1.formula} + {r2.formula}")
    print(products)


def oxidation_test():
    r = Compound("Ca", "Calcium", ComType.ELEMENT, ["Ca++"], 1, 7)
    product = oxidation(r)
    print(f"{r.formula} + O2")
    print(product)


if __name__ == "__main__":
    oxidation_test()
