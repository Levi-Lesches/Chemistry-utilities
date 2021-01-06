from chemicals import Molecule
from collections import Counter
from my_stuff.misc import init

def lcm (nums) -> int: 
    max_denom = float("-inf")
    denoms = set()
    for fraction in nums: 
        denoms.add(fraction.denom)
        if fraction.denom > max_denom: max_denom = fraction.denom

    result = max_denom
    while any (result % denom for denom in denoms): 
        result += max_denom

    return result

def expand_fractions(nullspace: list): 
    lcd = lcm (nullspace)
    return [
        abs (val.num * (lcd // val.denom))
        for val in nullspace
    ]

class Fraction: 
    @init
    def init (self, num, denom, negative): pass
    def __repr__(self): return f"{'-' if self.negative else ''}({self.num}/{self.denom})"
    def __new__(cls, num, denom): 
        # if num == 0: return 0
        if denom == 0: raise ZeroDivisionError()
        negative: bool = (num >= 0) != (denom >= 0)
        num = abs (num)
        denom = abs (denom)
        for n in range (num, 1, -1):  # simplify
            if not num % n and not denom % n: 
                num //= n
                denom //= n

        fraction: Fraction = object.__new__(cls)
        fraction.init(num, denom, negative)
        return fraction

    def to_float(self): return self.num / self.denom * (-1 if self.negative else 1)


class Matrix: 
    @init
    def __init__(self, matrix): self.rows, self.cols = self.set_shape()
    def __getitem__(self, p): return self.matrix [p[0]] [p[1]]
    def __eq__(self, other): return self.matrix == other.matrix
    def __len__(self): return self.rows * self.cols
    def __iter__(self): return iter (self.matrix)
    def __str__(self): return repr (self)
    def __repr__(self): return (
        "Matrix(" + 
        "".join (
            [
                f"\n\t{', '.join (map (str, row))}"
                for row in self
            ]
        ) + 
        "\n)"
    )

    def flatten(self): return [num for row in self for num in row]

    def set_shape(self): 
        rows = 0
        cols = None
        for row in self.matrix:
            if cols is None: cols = len (row)
            elif len (row) != cols: raise SyntaxError("Inconsistent dimensions")
            rows += 1
        return rows, cols

    def fromDimensions(rows, cols, List):
        result = []
        index = 0
        for row in range (rows):
            values = []
            for col in range (cols): 
                values.append (List [index])
                index += 1
            result.append (values)
        if index != len (List): raise SyntaxError("fromDimensions failed")
        return Matrix (result)

    def get_pivot (self, col: list) -> (int, int): 
        result_index = None
        result = None
        for index, value in enumerate (col): 
            if value != 0: return index, value
            else: 
                result_index = index
                result = value

        else:
            if all (num == 0 for num in col): return None, None
            else: return result_index, result

    def rref(self):
        cols = self.cols
        matrix = self.flatten()

        def get_col (index): return matrix [index::cols]

        def swap_rows (row1, row2): 
            index1_1 = row1 * cols
            index1_2 = (row1 + 1) * cols
            index2_1 = row2 * cols
            index2_2 = (row2 + 1) * cols
            matrix [index1_1 : index1_2], matrix [index2_1 : index2_2] = (
                matrix [index2_1 : index2_2], matrix [index1_1 : index1_2]
            )

        def cross_cancel (row, value, pivot_row, pivot_value): 
            offset = (pivot_row - row) * cols
            for n in range (row * cols, (row + 1) * cols): 
                matrix [n] = pivot_value * matrix [n] - value * matrix [n + offset]

        pivot_row = 0
        pivot_col = 0
        pivots = []
        while pivot_col < cols and pivot_row < self.rows:
            offset, value = self.get_pivot (get_col (pivot_col) [pivot_row:])

            if offset is None: 
                pivot_col += 1
                continue

            pivots.append (pivot_col)
            if offset != 0: swap_rows (pivot_row, offset + pivot_row)

            for row in range (self.rows): 
                if row == pivot_row: continue

                val = matrix [row * cols + pivot_col]
                if val == 0: continue
                else: cross_cancel (row, val, pivot_row, value)

            pivot_row += 1

        for index, col in enumerate (pivots):
            temp = index * cols + col
            value = matrix [temp]
            matrix [temp] = 1
            for index2 in range (temp + 1, (index + 1) * cols): 
                matrix [index2] = Fraction (matrix [index2], value)

        return Matrix.fromDimensions(self.rows, self.cols, matrix)

    def nullspace(self, simplify = True): 
        rref = self.rref()
        nullspace = [
            rref [n, -1] 
            for n in range (rref.rows)
        ]
        
            
        for index, value in enumerate (nullspace): 
            if value == 0: 
                nullspace [index] = (
                    Fraction (1, 1)
                    if simplify
                    else Fraction (0, 1)
                )
                    
        if not simplify: return nullspace
        nullspace.extend ([
            Fraction (1, 1) 
            for _ in range (rref.cols - len (nullspace))
        ])


        return expand_fractions(nullspace)


class Side: 
    def __init__ (self, formula: str):
        self.molecules: [Molecule] = self.get_molecules(formula)
        self.molecules_list = self.get_molecule_list(formula)
        self.elements: {"Element": int} = self.get_elements()

    def __repr__ (self): return f"Side ({self})"
    def __str__ (self): return " + ".join ([
        f"{count if count != 1 else ''}{molecule._base_molecule}" 
        for molecule, count in self.molecules.items()
    ])

    def get_molecule_list(self, formula): return list (self.molecules.keys())

    def get_molecules (self, formula: str) -> [Molecule]: return {
        Molecule (molecule._base_molecule): molecule.coefficient
        for molecule in map (Molecule, formula.split (" + "))
    }

    def get_elements (self) -> {"Element": int}: 
        elements: Counter = Counter()
        for molecule, coefficient, in self.molecules.items():
            for element, count in molecule.elements.items():
                elements [element] += count * coefficient
        return elements


class Equation: 
    def __init__(self, formula: str): 
        sides: [str] = formula.split (" --> ")
        self.left: Side = Side (sides [0])
        self.right: Side = Side (sides [1])
        self.verify()
        self.matrix = self.get_matrix()

    def __repr__(self): return f"Equation ({self})"
    def __str__(self): return f"{self.left} --> {self.right}"

    def is_balanced(self): 
        self.left.elements = self.left.get_elements()
        self.right.elements = self.right.get_elements()
        return self.left.elements == self.right.elements

    def verify(self) -> None: 
        if (
            any (
                element not in self.right.elements
                for element in self.left.elements
            ) or any (
                element not in self.left.elements
                for element in self.right.elements
            )
        ): raise SyntaxError (f"There is an inconsistency in {self}")

    def get_matrix(self) -> Matrix:
        matrix = []
        for element in self.left.elements: 
            row = []
            for molecule in self.left.molecules: 
                if element in molecule.elements: 
                    row.append (-molecule.elements [element])
                else: row.append (0)
            for molecule in self.right.molecules: 
                if element in molecule.elements: 
                    row.append (molecule.elements [element])
                else: row.append (0)
            matrix.append (row)
        return Matrix (matrix)

    def set_coefficients(self, nullspace: Matrix): 
        for index in range (len (self.left.molecules)):
            molecule = self.left.molecules_list [index]
            self.left.molecules [molecule] = nullspace [index]

        for index2 in range (len (self.right.molecules)):
            molecule = self.right.molecules_list [index2]
            self.right.molecules [molecule] = nullspace [index + index2 + 1]

    def balance(self) -> None: 
        self.set_coefficients(self.matrix.nullspace())
        if not self.is_balanced(): raise Exception (self)


def balance (input_: str) -> Equation: 
    equation: Equation = Equation (input_) 
    equation.balance()
    return equation

if __name__ == '__main__': 
    # from argparse import ArgumentParser
    # parser = ArgumentParser()
    # parser.add_argument ("formula", nargs = "*", help = "Equation to balance")
    # args = parser.parse_args()
    # eq = args.formula [0]
    eq = input("Enter a formula: ")
    # print(f"Balancing {eq}")
    print (balance (eq))