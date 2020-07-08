from numericalunits import eval, kg, m, s
assert eval("kg") == kg
assert eval("kg * m / s ** 2") == kg * m / s ** 2
