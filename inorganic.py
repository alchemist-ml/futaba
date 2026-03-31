import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

import csv
from collections import defaultdict
import math
import re
import json
import sys


def load_ions(filename: str):
    ions_dict = defaultdict(list)

    with open(filename, "r") as f:
        reader = csv.reader(f)
        next(reader)

        for row in reader:
            if row[0].startswith('"---'):
                continue

            try:
                ion = row[0].strip().strip('"')
                charge = int(row[1])
                if ion in ions_dict:
                    ions_dict[ion].append(charge)
                else:
                    ions_dict[ion].append(charge)
            except ValueError:
                continue

    return ions_dict


ion_charge = load_ions("Ions.csv")


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


def decomposition(compound):
    ions = extract_ions(compound)

    # carbonates decompose to metal oxide + CO2
    if "CO3" in ions:
        cations = [ion for ion in ions if ion != "CO3"]
        if cations:
            metal = cations[0]
            return [f"{metal}O", "CO2"]

    # metal hydroxides decompose to metal oxide + water
    if "OH" in ions:
        cations = [ion for ion in ions if ion != "OH"]
        if cations:
            metal = cations[0]
            return [f"{metal}O", "H2O"]

    # hydrate decomposition
    # e.g., CuSO4·5H2O → CuSO4 + 5H2O
    if "·" in compound:
        parts = compound.split("·")
        anhydrous = parts[0]
        water_count = ""
        if len(parts) > 1:
            water_count = parts[1][0] if parts[1][0].isdigit() else "1"
        return [anhydrous, "H2O"]

    # chlorates decompose to chloride + oxygen (e.g. 2KClO3 → 2KCl + 3O2)
    if "ClO3" in ions or "ClO4" in ions:
        anion_type = "ClO4" if "ClO4" in ions else "ClO3"
        cations = [ion for ion in ions if anion_type not in ion]
        if cations:
            metal = cations[0]
            return [f"{metal}Cl", "O2"]

    # metal nitrates decompose to metal nitrite + oxygen
    if "NO3" in ions:
        cations = [ion for ion in ions if ion != "NO3"]
        if cations:
            metal = cations[0]
            return [f"{metal}NO2", "O2"]

    # metal sulfites decompose to metal oxide + SO2
    if "SO3" in ions and "SO4" not in ions:
        cations = [ion for ion in ions if ion != "SO3"]
        if cations:
            metal = cations[0]
            return [f"{metal}O", "SO2"]

    decomposed_ions = []
    for i in ions:
        if len(i) > 1 and i.isupper():
            decomposed_ions += i.split()
        decomposed_ions += i
    for d in decomposed_ions:
        if d.isdigit():
            decomposed_ions.remove(d)

    decomposed_ions = list(set(decomposed_ions))
    print(decomposed_ions)

    results = []
    for i in range(len(decomposed_ions)):
        for j in range(i + 1, len(decomposed_ions)):
            results += shuffle([decomposed_ions[i], decomposed_ions[j]])

    for r in results:
        if r == compound:
            results.remove(r)

    return results


def tidy(arg):
    arg = arg.replace("+", "").replace("-", "").replace("1", "")
    product = "H2O" if arg == "HOH" else arg
    return product


def shuffle(reactants: list):

    reactant_ions = []
    for r in reactants:
        reactant_ions += extract_ions(r)

    species = [{i: ion_charge[i][0]} for i in reactant_ions]
    no_of_species = len(species)

    products = []
    for i in range(no_of_species):
        for j in range(no_of_species):
            if species[i] == species[j]:
                continue

            i1, c1 = next(iter(species[i].items()))
            i2, c2 = next(iter(species[j].items()))
            if (c1 > 0 and c2 > 0) or (c1 < 0 and c2 < 0):
                continue

            cation, anion = (i1, i2) if c1 > c2 else (i2, i1)
            if c1 < c2:
                c1, c2 = c2, c1

            gcd = math.gcd(c1, -1 * c2)
            c1 = int(c1 / gcd)
            c2 = int(c2 / gcd)

            if len(cation) > 1 and cation.isupper() and abs(c2) > 1:
                cation = f"({cation})"
            if len(anion) > 1 and anion.isupper() and abs(c1) > 1:
                anion = f"({anion})"

            products.append(cation + str(c2) + anion + str(c1))
            products = list(set(products))

    # special handling for metal-acid reactions to produce H2 (eg. Zn + HCl -> ZnCl + H2)
    metals = ["Zn", "Mg", "Fe", "Al", "Ca", "Cu", "Pb", "Ni", "Sn", "Cr", "Mn"]
    has_metal = any(metal in reactants for metal in metals)
    has_acid = any("H" in r for r in reactants)
    has_hydrogen_in_products = any("H" in p for p in products)

    if has_metal and has_acid and has_hydrogen_in_products:
        products_to_remove = []
        for p in products:
            if "H" in p and p != "H2":
                products_to_remove.append(p)

        for p in products_to_remove:
            if p in products:
                products.remove(p)

        if "H2" not in products:
            products.append("H2")

    # metals combust to form metal oxide
    if len(reactants) == 2 and "O2" in reactants:
        non_oxygen = (
            [r for r in reactants if r != "O2"][0]
            if any(r != "O2" for r in reactants)
            else None
        )
        if non_oxygen and non_oxygen in metals:
            products.clear()
            products.append(f"{non_oxygen}O")

    # nonmetal combustion (e.g. S + O2 → SO2, C + O2 → CO2)
    nonmetals = ["C", "S", "P", "N"]
    if len(reactants) == 2 and "O2" in reactants:
        non_oxygen = (
            [r for r in reactants if r != "O2"][0]
            if any(r != "O2" for r in reactants)
            else None
        )
        if non_oxygen and non_oxygen in nonmetals:
            products.clear()
            if non_oxygen == "C":
                products.append("CO2")
            elif non_oxygen == "S":
                products.append("SO2")
            elif non_oxygen == "P":
                products.append("P2O5")
            elif non_oxygen == "N":
                products.append("NO2")

    # metal + nonmetal → salt
    if len(reactants) == 2:
        r1, r2 = reactants[0], reactants[1]
        r1_is_metal = r1 in metals
        r2_is_metal = r2 in metals
        r1_is_nonmetal = r1 in nonmetals
        r2_is_nonmetal = r2 in nonmetals

        # metal + monmetal combination
        if (r1_is_metal and r2_is_nonmetal) or (r2_is_metal and r1_is_nonmetal):
            metal = r1 if r1_is_metal else r2
            nonmetal = r2 if r1_is_metal else r1

            if nonmetal == "Cl":
                products = [f"{metal}Cl"]
            elif nonmetal == "Br":
                products = [f"{metal}Br"]
            elif nonmetal == "I":
                products = [f"{metal}I"]
            elif nonmetal == "S":
                products = [f"{metal}S"]
            elif nonmetal == "O":
                products = [f"{metal}O"]
            elif nonmetal == "N":
                products = [f"{metal}N"]

    # halogen displacement
    halogens = ["F", "Cl", "Br", "I"]
    halogen_order = ["F", "Cl", "Br", "I"]  # reactivity order

    for halogen in halogens:
        if halogen in reactants:
            other = [r for r in reactants if r != halogen]
            if other:
                compound = other[0]

                for other_halogen in halogens:
                    if other_halogen in compound:
                        if halogen_order.index(halogen) < halogen_order.index(
                            other_halogen
                        ):
                            products.clear()
                            new_compound = compound.replace(other_halogen, halogen)
                            displaced = other_halogen
                            products.append(new_compound)
                            products.append(displaced)
                        break

    final_products = []
    for p in products:
        if tidy(p) in reactants:
            continue
        final_products.append(tidy(p))
    return final_products


def test():
    acid_base_tests = [
        ["HCl", "NaOH"],
        ["H2SO4", "KOH"],
        ["HNO3", "Ca(OH)2"],
        ["H3PO4", "NaOH"],
        ["H2CO3", "Mg(OH)2"],
    ]
    metal_acid_tests = [
        ["Zn", "HCl"],
        ["Mg", "H2SO4"],
        ["Fe", "HCl"],
        ["Al", "HNO3"],
        ["Ca", "H2SO4"],
    ]
    double_displacement_tests = [
        ["NaCl", "AgNO3"],
        ["BaCl2", "Na2SO4"],
        ["KBr", "Pb(NO3)2"],
        ["CaCl2", "Na2CO3"],
        ["Na2S", "Cd(NO3)2"],
    ]
    synthesis_tests = [
        ["H", "O"],
        ["Na", "Cl"],
        ["N", "H"],
        ["Ca", "O"],
    ]
    combustion_tests = [
        ["C", "O2"],
        ["S", "O2"],
        ["Na", "O2"],
        ["Mg", "O2"],
    ]
    halogen_displacement_tests = [
        ["Cl2", "NaBr"],
        ["Br2", "NaCl"],
        ["F2", "NaCl"],
    ]

    decomposition_tests = [
        "CaCO3",
        "Mg(OH)2",
        "KClO3",
        "Ca(NO3)2",
    ]

    print("\nAcid base test:")
    for reactants in acid_base_tests:
        print(
            f"{reactants[0]} + {reactants[1]} -> possible results: {shuffle(reactants)}"
        )

    print("\nMetal acid test:")
    for reactants in metal_acid_tests:
        print(
            f"{reactants[0]} + {reactants[1]} -> possible results: {shuffle(reactants)}"
        )

    print("\nDouble displacement test:")
    for reactants in double_displacement_tests:
        print(
            f"{reactants[0]} + {reactants[1]} -> possible results: {shuffle(reactants)}"
        )

    print("\nSynthesis/Combination test:")
    for reactants in synthesis_tests:
        print(
            f"{reactants[0]} + {reactants[1]} -> possible results: {shuffle(reactants)}"
        )

    print("\nCombustion test:")
    for reactants in combustion_tests:
        print(
            f"{reactants[0]} + {reactants[1]} -> possible results: {shuffle(reactants)}"
        )

    print("\nHalogen displacement test:")
    for reactants in halogen_displacement_tests:
        print(
            f"{reactants[0]} + {reactants[1]} -> possible results: {shuffle(reactants)}"
        )

    print("\nDecomposition test:")
    for compound in decomposition_tests:
        print(f"{compound} -> {react([compound])}")


def react(reactant: list):
    if len(reactant) == 1:
        return decomposition(reactant[0])
    else:
        return shuffle(reactant)


if __name__ == "__main__":
    print(json.dumps({"products": react(sys.argv[1:])}))
