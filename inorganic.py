import csv
from collections import defaultdict
import math
import re
from chemlib import Compound


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


def react(reactant: list):
    if len(reactant) == 1:
        return decomposition(reactant[0])
    else:
        return shuffle(reactant)


if __name__ == "__main__":
    # print(react(["Mg(OH)2", "H2SO4"]))
    ...
