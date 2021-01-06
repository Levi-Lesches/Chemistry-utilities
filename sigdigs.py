# BUG: Eye on main for loop 

from operator import add, sub, mul, truediv as div
from my_stuff.nums import is_even
from my_stuff.strs import find_closing_paren
import re


OPERATORS = {
	"+": add,
	"-": sub,
	"*": mul,
	"/": div
}

def get_sigdigs (num):
	num = num.lstrip ("0")
	decimal = False
	count = 0
	last_sigdig: int = [
		index
		for index, digit in enumerate (num)
		if digit.isdigit()
	] [-1]

	for index, digit in enumerate (num):
		if digit == ".": decimal = True
		elif digit != "0" or (count and (index < last_sigdig or decimal)): count += 1

	return count

def refine_sigdigs_after_decimal (num: str, sigdigs = None, decimals = None):
	if num.find ("e") == -1 and decimals is not None: 
		return f"{float (num):.{decimals}f}"

	trim = False
	decimal = False
	count = 0
	int_count = 0
	for index, digit in enumerate (num): 
		if (
			decimals is not None and 
			count > decimals or 
			sigdigs is not None and 
			(count + int_count) > sigdigs
		): trim = index
		if digit == ".": decimal = True
		elif digit == "e": break
		elif decimal: count += 1
		else: int_count += 1
	else: index += 1

	if trim: 
		index_of_e = num.find ("e")
		if index_of_e == -1: 
			if not decimals: 
				return num [:num.find (".")]
			else: 
				return num [:trim]
		else: 
			if not decimals: 
				return num [:num.find (".")] + num [index_of_e:]
			else: return num [:trim] + num [index_of_e:]

	num = list (num)
	decimals = sigdigs - (count + int_count) if sigdigs else decimals
	if decimals and not count: 
		assert not decimal
		num.insert (index, ".")
		index += 1
	for _ in range (decimals):
		num.insert (index, "0")
		index += 1
	return "".join (num)


def calculate (a: str, b: str, operator: str) -> str:
	if operator in ("*", "/"): 
		sigdigs = min (map (get_sigdigs, (a, b)))
		result: float = OPERATORS [operator] (float (a), float (b))
		return refine_sigdigs_after_decimal (
			f"{result:.{sigdigs}g}", # 10 ^ 0
			sigdigs = sigdigs
		)
	else: 
		decimals: [int] = [0, 0]
		for index, num in enumerate ((a, b)):
			index_of_decimal: int = num.find (".")
			if index_of_decimal != -1: 
				decimal_spots: int = len (num [index_of_decimal + 1:])
				decimals [index] = decimal_spots
		decimal_places: int = min (decimals)
		result: float = OPERATORS [operator] (float (a), float (b))
		return refine_sigdigs_after_decimal (
			# f"{result:.{decimal_places}g}",
			str (result),
			decimals = decimal_places
		)


regex_parens = re.compile(r"\((.+)\) ([\+\-\*\/]) \((.+)\)")
regex_expression = re.compile(r"([\d\.]+) ([\+\-\*\/]) ([\d\.]+)")

def parse_expression (expression: str) -> str: 
	if not all (
		digit.isdigit() or 
		digit.isspace() or 
		digit in OPERATORS or
		digit in (".", "(", ")")
		for digit in expression
	): raise SyntaxError ("Invalid expression")

	match_parens = regex_parens.match(expression)
	if match_parens is not None: 
		left = match_parens [1]
		operator = match_parens [2]
		right = match_parens [3]
		return calculate(parse_expression(left), parse_expression(right), operator)
	else: 
		match_expression = regex_expression.match(expression)
		if match_expression is None: 
			raise SyntaxError ("Invalid expression")
		left = match_expression [1]
		operator = match_expression [2]
		right = match_expression [3]
		return calculate(left, right, operator)

if __name__ == "__main__": 
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument ("expression", nargs = "?", help="Expression to evaluate")
	args = parser.parse_args()
	expression = args.expression
	if not expression: expression = input ("Enter the expression: ")
	
	print (parse_expression (expression))