import numpy as np
import thermochem as tc
import periodictable as pt

#### TEMP STRUCT

class molecule:
    class atom:
        def __init__(self):
            self.num: float = 0
            self.pt_data = None        

    def __init__(self):
        self.formula: list[molecule.atom] = []
        self.num_elements: int = 0
        self.molar_mass: float = 0 # g/mol
        self.__all_elements_list__: list = []
        self.__defined__: bool = False
        for el in pt.elements:
            self.__all_elements_list__.append(el.symbol)

    def define_molecule(self, atom_list: list) -> bool:
        if self.__defined__ is True:
            return False
        for atom, num in atom_list:
            new_atom: molecule.atom = molecule.atom()
            new_atom.num = num
            if atom not in self.__all_elements_list__:
                return False
            new_atom.pt_data = pt.elements.symbol(atom)
            self.formula.append(new_atom)
            self.num_elements = self.num_elements + 1
            self.molar_mass += new_atom.num * new_atom.pt_data.mass
        self.__defined__ = True
        return True
    
    def get_element_list(self) -> list:
        return self.__all_elements_list__
    
    def check_for_element(self,element_symbol) -> bool:
        for atom in self.formula:
            if atom.pt_data.symbol is element_symbol:
                return True
        return False
    
class chem_eq:
    def __init__(self):
        self.reactants: list[tuple(molecule,float)] = []
        self.products: list[tuple(molecule,float)] = []
        self.num_elements: int = 0
        self.__defined__ = False
        self.__solved__ = False

    def define_formula(self, reactants: list, products: list):
        if self.__defined__ is True:
            return False
        self.reactants = self.separate_elements(reactants)
        self.products = self.separate_elements(products)
        if not self.check_elements():
            self.reactants = []
            self.products = []
            return False
        return True

    def separate_elements(self, el_list: list):
        sep_el_list = []
        for mol, coef in el_list:
            el = ""
            num = 0
            atom_list: list[tuple(str,float)] = []
            i = 0
            while i <= len(mol):
                if i >= len(mol):
                    num = float(1)
                    atom_list.append((el,num))
                elif mol[i].isupper() is True and len(el) > 0:
                    num = float(1)
                    atom_list.append((el,num))
                    el = mol[i]
                elif mol[i].isnumeric() is True or mol[i] == '.':
                    it = 1
                    if i + it < len(mol):
                        while mol[i + it].isnumeric() is True or mol[i + it] == '.':
                            it += 1
                            if i + it >= len(mol):
                                break
                    num = float(mol[i:i+it])
                    i += it - 1
                    atom_list.append((el,num))
                    if i >= len(mol):
                        break
                    el = ""
                    num = 0
                else:
                    el = el + mol[i]
                i += 1
            new_mol = molecule()
            new_mol.define_molecule(atom_list)
            pair = (new_mol,coef)
            sep_el_list.append(pair)
        return sep_el_list

    def check_elements(self):
        reactants = []
        for mol, n in self.reactants:
            for atom in mol.formula:
                if atom.pt_data.symbol not in reactants:
                    reactants.append(atom.pt_data.symbol)
        products = []
        for mol, n in self.products:
            for atom in mol.formula:
                products.append(atom.pt_data.symbol)
        for el in reactants:
            if el not in products:
                return False
        self.num_elements = len(reactants)
        return True

    def solve_eq(self):
        unknown = np.zeros([self.num_elements, len(self.products)])
        known = np.zeros(self.num_elements)
        el_list = []
        unknown_eq_index = 0
        
        for molecule_pair in self.reactants:
            for el in molecule_pair[0].formula:
                if el.pt_data.symbol not in el_list:
                    el_list.append(el.pt_data.symbol)
                el_index = el_list.index(el.pt_data.symbol)
                known[el_index] += molecule_pair[1] * el.num
        for molecule_pair in self.products:
            for el in el_list:
                el_index = el_list.index(el)
                if not molecule_pair[0].check_for_element(el):
                    continue
                formula_index = 0
                for i in range(len(molecule_pair[0].formula)):
                    if molecule_pair[0].formula[i].pt_data.symbol is el:
                        formula_index = i
                        break
                unknown[el_index][unknown_eq_index] = molecule_pair[0].formula[formula_index].num
            unknown_eq_index += 1

        solution = np.linalg.solve(unknown, known)

        print("Known")
        print(known)
        print("Unknown")
        print(unknown)
        print("Solution")
        print(solution)

        for i in range(len(solution)):
            self.products[i] = (self.products[i][0],solution[i])

        return




#### Functions


# TODO: add define_formula() w/ [(molecule,coef)] input instead of string

#### TEST MAIN

reactants = [("Na6Cl12",12)]
products = [("Na123",0),("Cl12",0)]

# reactants = [("C213.8H323O4.6N2.3",1),("NH4ClO4",102.114)]
# # reactants = [("C2H4")]
# products = [("HCl",0),("CO",0),("H2",0),("H2O",0),("N2",0)]


test = chem_eq()
test.define_formula(reactants,products)
# print(test.reactants[1][0].molar_mass)
test.solve_eq()
None
