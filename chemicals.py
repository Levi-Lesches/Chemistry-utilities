from argparse import ArgumentParser
from collections import Counter
from functools import reduce
from operator import add 
from textwrap import dedent

from roman import convert as to_roman
from elements import ELEMENTS as periodic_table
from my_stuff.misc import pause, init
from my_stuff.strs import find_closing_paren, consume_int
from my_stuff.lists import VOWELS

PREFIXES: {int, str} = {
	1: "Mono",
	2: "Di",
	3: "Tri",
	4: "Tetra",
	5: "Penta",
	6: "Hexa",
	7: "Hepta",
	8: "Octa",
	9: "Nona",
	10: "Deca"
} 

POLYATOMIC_IONS: {str: (str, int)} = { # formula: (name, charge)
	"C2H3O2": ("Acetate", -1),
	"CO3": ("Carbonate", -2),
	"PO4": ("Phosphate", -3),
	"OH": ("Hydroxide", -1),
	"NO2": ("Nitrite", -1),
	"NO3": ("Nitrate", -1),
	"HCO3": ("Hydrogen carbonate", -1),
	"CrO4": ("Chromate", -2),
	"CrO7": ("Dichromate", -2),
	"HPO4": ("Hydrogen phosphate", -2),
	"H2PO4": ("Dihydrogen phosphate", -1),
	"NH4": ("Ammonium", 1),
	"SO3": ("Sulfite", -2),
	"SO4": ("Sulfate", -2), 
	"HSO3": ("Hydrogen sulfite", -1),
	"HSO4": ("Hydrogen sulfate", -1),
	"MnO4": ("Permanganate", -1),
	"CN": ("Cyanide", -1),
	"O2": ("Peroxide", -2),
	"Hg2": ("Mercury(I)", 2)
}

# some special polyatomic ions
_flourine = "F", "fluor"
_chlorine = "Cl", "chlor"
_bromine = "Br", "brom"
_iodine = "I", "iod"
_astatine = "At", "ast"

for symbol, base_name in (_flourine, _chlorine, _bromine, _iodine, _astatine):
	POLYATOMIC_IONS [f"{symbol}O"] = (f"Hypo{base_name}ite", -1)
	POLYATOMIC_IONS [f"{symbol}O2"] = (f"{base_name.title()}ite", -1)
	POLYATOMIC_IONS [f"{symbol}O3"] = (f"{base_name.title()}ate", -1)
	POLYATOMIC_IONS [f"{symbol}O4"] = (f"Per{base_name}ate", -1)

args = ArgumentParser()
args.add_argument ("formula", nargs = "?", help = "Formula to parse")
args = vars (args.parse_args())

def counter_mul (counter: Counter, mul: int):
	for element in counter: counter [element] *= mul

def get_english_elements (elements: dict) -> dict:
	try: english_elements = {
		periodic_table [element].name: count 
		for element, count in elements.items()
	}
	except Exception as error:
		invalid = error.__context__.__context__
		raise KeyError (f"{invalid} is not on the periodic table!") from None
	else: return english_elements

def get_gen_index (generator, index: int):
	tup = tuple (generator)
	return tup [index]

class Molecule:
	@init
	def __init__ (self, formula): 
		self.elements = {
			periodic_table [name]: count 
			for name, count in self.get_elements().items()
		}
		# self.english_elements: dict = get_english_elements (self.elements)
		self.english_elements:dict = {
			element.name: count 
			for element, count in self.elements.items()
		}
		self.elements_list: tuple = tuple (self.elements.keys())
		self.refined_formula = self.refine_formula()
		self.mass: int = self.get_mass()
		# self.ionic: bool = self.get_ionic()
		self.type: str = self.get_type()
		self.name: str = self.get_name()
		self.coefficient, index = consume_int (self.formula, 0)
		if self.coefficient is None: self.coefficient = 1
		self._base_molecule = self.formula [index:]


	def __str__(self): return dedent (f"""
		Molecule ({self.formula}):
			Elements: {self.english_elements}
			Mass: {self.mass}
			Type of Molecule: {self.type.title()}
			Name of Molecule: {self.name}
	""")

	def __repr__(self): return f"Molecule ({self.formula})"

	def __eq__(self, other): return type (other) is Molecule and self.formula == other.formula
	def __hash__(self): return hash (self.formula)

	def get_elements (self, formula = None):
		if formula is None: formula = self.formula
		last_element = ""
		last_num = None
		elements = Counter()
		sub_elements = Counter()
		skips = []
		mul = 1
		open_paren_index = -1 # -1 + 1 = 0
		# handle nested parens
		for _ in range (formula.count ("(")):
			open_paren_index = formula.find ("(", open_paren_index + 1)
			if open_paren_index in skips: continue
			elif open_paren_index != -1: 
				close_paren_index = find_closing_paren (formula, open_paren_index)
				if close_paren_index is None: 
					raise ValueError (f"Could not find closing parenthesis for the {_} opening parenthesis -- {formula}")
				new_elements: Counter = self.get_elements (
					formula [open_paren_index + 1:close_paren_index]
				)
				# find coefficient
				temp_mul, index = consume_int (formula, close_paren_index + 1)
				if temp_mul is None: temp_mul = 1
				counter_mul (new_elements, temp_mul)
				# prep
				skips.extend (range (open_paren_index, index))
				sub_elements += new_elements
			else: raise ValueError (f"ERROR: Please tell Levi you got this error with this chemical formula: {formula}")

		# parse the formula
		for index, letter in enumerate (formula):
			if index in skips: continue
			elif index == 0 and letter.isdigit(): 
				mul, end_index = consume_int (formula, 0)
				skips.extend (range (end_index))
				continue
			elif letter.islower(): 
				if not last_element: #starts with a lowercase?!
					raise ValueError (f"Caught an element that starts with a lowercase letter: {formula}")
				elif last_element [-1].islower(): #already a lowercase letter in the last element!
					raise ValueError (f"Cannot have two lowercase letters (letters {index} and {index + 1}) together: {formula}")
				else: last_element += letter
			elif letter.isupper():
				if last_element: elements [last_element] += 1 if last_num is None else last_num

				last_element = letter
				last_num = None
			elif letter.isdigit(): # subscripts
				if last_num is None: last_num = int (letter)
				else: last_num = int (str (last_num) + letter)
			else: raise ValueError (f"{letter} at pos {index + 1} was unexpected: {formula}")

		# cleanup
		if last_element: 
			elements [last_element] += last_num if last_num is not None else 1

		if sub_elements is not None: elements += sub_elements # add paren contents
		counter_mul (elements, mul) # distribute coefficient
		return elements

	def get_mass (self): return reduce (add, 
		(
			element.mass * count
			for element, count in self.elements.items()
		)
	)

	def refine_formula (self): 
		formula = self.formula[:]
		while formula [0].isdigit() or formula [0] == "(":
			formula = formula [1:]
		while formula [-1] == ")": formula = formula [:-1]
		return formula

	def get_type(self): 
		first = self.elements_list [0]
		if len (self.elements) == 1: return "element"
		elif first.series in (3, 4, 7, 8) or self.refined_formula.startswith ("NH4"): 
			return "ionic"
		else: return "molecular"

	def get_name (self):
		if self.type == "molecular": return self.name_molecular()
		elif self.type == "ionic": return self.name_ionic()
		elif self.type == "acid": return self.name_acid()
		elif self.type == "element": return self.name_element()
		else: raise ValueError (f"{self!r} has an invalid type -- {self.type}")

	def name_molecular (self): 
		assert len (self.elements) != 1, (
			f"Misclassified {self!r} as {self.type} instead of 'element'"
		)

		if len (self.elements) != 2: 
			return "Unknown (Cannot name compounds with more than 2 elements)"
		element1, element2 = self.elements
		count1 = self.elements [element1]
		count2 = self.elements [element2]
		assert count1 > 0 and count2 > 0, f"{self!r} got incorrect count for elements"
		if count1 > 10 or count2 > 10:
			return "Unknown (Element count greater than 10 not supported)"
		prefix1 = PREFIXES [count1]
		prefix2 = PREFIXES [count2]
		# prime names for prefixes
		name1 = element1.name.lower().lstrip (prefix1 [-1])
		name2 = element2.name.lower().lstrip (prefix2 [-1])
		base_name = get_base_name (name2)
		if prefix1 == "Mono": 
			prefix1 = ""
			name1 = name1.title()

		if self.refined_formula in POLYATOMIC_IONS: 
			special_name = f" ({POLYATOMIC_IONS [self.refined_formula] [0]})"
		else: special_name = ""

		return f"{prefix1}{name1} {prefix2}{base_name}ide{special_name}"

	def name_ionic (self): 
		cation = self.elements_list [0]
		anion = self.elements_list [1:]
		# get cation
		charge = ""
		if self.refined_formula.startswith ("NH4"): 
			cation_name = "Ammonium"
			if len (self.elements_list) <= 2: return "Ammonium"
			else: anion = self.elements_list [-1]
		else: 
			cation_charges = get_charges (cation) 
			if len (cation_charges) > 1: 
				anion_charges = get_charges (anion)
				if (
						anion_charges is None or 
						len (anion_charges) > 1 or 
						len (anion) > 1
				): # rest is a polyatomic ion

					polyatomic_ion = self.get_polyatomic_ion (cation)
					if polyatomic_ion is None:
						return "Unknown (Cannot determine charge of cation)"
					anion_name, polyatomic_charge, polyatomic_count = polyatomic_ion
					try: charge = calc_charge (
						cation_charges,
						self.elements [cation],
						polyatomic_charge,
						polyatomic_count
					)
					except IndexError: return "Unknown (Molecule is not neutral)"
					if charge is None: 
						return "Unknown (Cannot determine charge of cation)"
					else: charge = f"({to_roman (charge)})"
				else: 
					anion = anion [0]
					anion_charge = anion_charges [0]
					try: cation_charge = calc_charge (
						cation_charges,
						self.elements [cation], # cation count
						anion_charge,
						self.elements [anion] # anion count
					)
					except IndexError: return "Unknown (Molecule is not neutral)"
					if cation_charge is None: 
						return "Unknown (Cannot determine charge of cation)"
					charge = f"({to_roman (cation_charge)})"
			cation_name = f"{cation.name.title()}{charge}"

		# get anion
		if type (anion) is not tuple or len (anion) == 1:
			if type (anion) is tuple: anion = anion [0]
			anion_name = f"{get_base_name (anion.name).title()}ide"
		else: # Must be a polyatomic Ionic
			polyatomic_ion = self.get_polyatomic_ion (cation)
			if polyatomic_ion is None: 
				return f"Unknown (Unrecognized polyatomic ion -- {self.formula})"
			else: anion_name = polyatomic_ion [0]

		return f"{cation_name} {anion_name}"

	def name_acid (self): return "Unknown (Acids not supported)"

	def name_element(self): 
		if self.refined_formula in POLYATOMIC_IONS: 
			special_name = f" ({POLYATOMIC_IONS [self.refined_formula] [0]})"
		else: special_name = ""

		return self.elements_list [0].name + special_name

	def get_polyatomic_ion (self, cation):
		index = self.formula.find (cation.symbol) + len (cation.symbol)  
		polyatomic_ion = ""
		break_next = False
		while True: 
			try: letter = self.formula [index]
			except IndexError: 
				if break_next: break
				else: return None
			else: 
				if letter == "(": 
					index += 1
					continue
				elif not polyatomic_ion and letter.isdigit(): 
					index += 1
					continue
				elif break_next and letter.isdigit(): polyatomic_ion += letter
				elif break_next: break
				else: polyatomic_ion += letter
				index += 1
				if polyatomic_ion in POLYATOMIC_IONS: 
					break_next = True
		if len (self.formula) > index and self.formula [index] == ")": index += 1
		count = consume_int (self.formula, index) [0]
		if count is None: count = 1
		return POLYATOMIC_IONS [polyatomic_ion] + (count,)

def get_charges (element): 
	# for class purposes
	if type (element) is tuple: 
		if len (element) > 1: return None
		else: element = element [0]
	if element.symbol == "H": return [1, -1]
	if element.group == 1: return [1]
	elif element.group == 2: return [2]
	elif element.symbol in ("Sc", "Y", "La", "Ac", "Al"): return [3]
	elif element.symbol in ("N", "P", "As"): return [-3]
	elif element.symbol in ("O", "S", "Se", "Te"): return [-2]
	elif element.symbol in ("F", "Cl", "Br", "I", "At"): return [-1]
	elif element.series == 2: return [] # the noble gasses
	else: 
		result = []
		append = False
		for charge in element.oxistates.split (", "):
			if "*" in charge: append = True
			if append: result.append (int (charge.strip ("*")))
		return [charge for charge in result if charge]

def calc_charge (
		possible_charges: list, 
		count: int, 
		other_charge: int, 
		other_count: int,
		silent: bool = False
) -> int: 
	target_charge = other_charge * other_count #O2 -> -2
	target_charge *= -1 # gotta cancel out that charge!
	for charge in possible_charges:
		if target_charge / charge == count: return charge
	else: 
		if silent: return None
		else: raise IndexError (
			f"No possible charges can be computed for {count} ions to cancel out a"
			f" charge of {target_charge * -1}"
		)

def get_base_name (name: str) -> str: 

	if name.endswith ("orus"): base = name [:-4]
	elif name.endswith ( ("ium", "ine", "gen") ): base = name [:-3]
	elif name.endswith ( ("on", "ur", "ic") ): base = name [:-2]
	else: return None

	if base [-1] in VOWELS + ["y"]: base = base [:-1]

	return base


if __name__ == "__main__":
	if not args ["formula"]: # we were double clicked!
		formula:str = input ("What is the chemical formula? ")
		should_pause:bool = True
	else: #run on cmd
		formula:str = args ["formula"]
		should_pause:bool = False

	print (Molecule (formula))
	if should_pause: pause()