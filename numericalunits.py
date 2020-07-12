# -*- coding: utf-8 -*-
"""
For information and usage see README, or http://pypi.python.org/pypi/numericalunits
"""
# Copyright (C) 2012-2020 Steven J. Byrnes
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from math import pi

__version__ = 1.24


########## Set all variables, to help introspection libraries ################
# This part is functionally pointless, it only helps IDE autocompletion
# / introspection libraries know that these variables exist. The actual
# values are set below, using the "global" keyword inside functions.

m = kg = s = C = K = 0.
cm = mm = um = nm = pm = fm = km = angstrom = Å = lightyear = \
    astro_unit = pc = kpc = Mpc = Gpc = inch = foot = mile = thou = 0.
L = mL = uL = nL = pL = fL = aL = kL = ML = GL = 0.
ms = us = ns = ps = fs = minute = hour = day = week = year = 0.
Hz = mHz = kHz = MHz = GHz = THz = PHz = rtHz = rpm = 0.
Hz·2π = mHz·2π = kHz·2π = MHz·2π = GHz·2π = THz·2π = PHz·2π = rpm·2π = 0.
g = mg = ug = ng = pg = fg = tonne = amu = Da = kDa = lbm = 0.
J = mJ = uJ = nJ = pJ = fJ = kJ = MJ = GJ = erg = eV = meV = keV = MeV = GeV = \
    TeV = btu = smallcal = kcal = Wh = kWh = 0.
NA = mol = mmol = umol = nmol = pmol = fmol = M = mM = uM = nM = pM = fM = 0.
N = mN = uN = nN = pN = fN = kN = MN = GN = dyn = lbf = 0.
Pa = hPa = kPa = MPa = GPa = bar = mbar = cbar = dbar = kbar = Mbar = atm = \
    torr = mtorr = psi = 0.
W = mW = uW = nW = pW = kW = MW = GW = TW = horsepower_imperial = \
    horsepower_metric = 0.
Gal = mGal = uGal = eotvos = 0.
degFinterval = degCinterval = mK = uK = nK = pK = 0.
mC = uC = nC = Ah = mAh = 0.
A = mA = uA = nA = pA = fA = 0.
V = mV = uV = nV = kV = MV = GV = TV = 0.
ohm = mohm = kohm = Mohm = Gohm = Ω = mΩ = kΩ = MΩ = GΩ = S = mS = uS = nS = 0.
T = mT = uT = nT = G = mG = uG = kG = Oe = Wb = 0.
F = uF = nF = pF = fF = aF = H = mH = uH = nH = 0.
c0 = mu0 = μ0 = eps0 = ε0 = Z0 = hPlanck = hbar = ħ = kB = GNewton = sigmaSB = \
    σSB = alphaFS = αFS = 0.
Rgas = e = uBohr = uNuc = aBohr = me = mp = mn = Rinf = Ry = Hartree = \
    ARichardson = Phi0 = KJos = RKlitz = 0.
REarth = g0 = Msolar = MEarth = 0.

########################### Main code #######################################


def reset_units(seed=None):
    """
    Set all units to new, self-consistent, floating-point values. See package
    documentation for detailed explanation and examples:
    http://pypi.python.org/pypi/numericalunits

    reset_units() --> units are randomized. This is the suggested use, and is
    done automatically the first time the module is imported. So you don't need
    to call this function explicitly; just do your calculation, display the
    final answer, then repeat in a fresh Python session. If you get the same
    answer both times, then your calculations are almost guaranteed to be free of
    dimensional-analysis-violating errors.

    reset_units('SI') --> Set units so that all values are given in standard SI
    units (meters-kilograms-seconds) by default. In this mode, there is no way
    to test for dimensional-analysis-violating errors.

    reset_units(x) --> If you pass any other argument x, it's used as the seed
    for the random number generator.
    """
    import random

    global m, kg, s, C, K

    if seed == 'SI':
        m = 1.
        kg = 1.
        s = 1.
        C = 1.
        K = 1.
    else:
        prior_random_state = random.getstate()

        if seed is None:
            random.seed()
        else:
            random.seed(seed)

        m = 10 ** random.uniform(-2,2) # meter
        kg = 10 ** random.uniform(-2,2) # kilogram
        s = 10 ** random.uniform(-2,2) # second
        C = 10 ** random.uniform(-2,2) # coulomb
        K = 10 ** random.uniform(-2,2) # kelvin

        # Leave the random generator like I found it, in case something else is
        # using it.
        random.setstate(prior_random_state)

    set_derived_units_and_constants()

def set_derived_units_and_constants():
    """
    Assuming that the base units (m, kg, s, C, K) have already been set as
    floating-point values, this function sets all other units and constants
    to the appropriate, self-consistent values.
    """

    # Length
    global cm, mm, um, nm, pm, fm, km, angstrom, Å, lightyear, \
        astro_unit, pc, kpc, Mpc, Gpc, inch, foot, mile, thou
    cm = 1e-2 * m
    mm = 1e-3 * m
    um = 1e-6 * m
    nm = 1e-9 * m
    pm = 1e-12 * m
    fm = 1e-15 * m
    km = 1e3 * m
    angstrom = 1e-10 * m
    Å = angstrom # easier-to-read alias (see https://sjbyrnes.com/unicode.html )
    lightyear = 9460730472580800. * m
    astro_unit = 149597870700. * m # astronomical unit
    pc = (648000./pi) * astro_unit # parsec
    kpc = 1e3 * pc
    Mpc = 1e6 * pc
    Gpc = 1e9 * pc
    inch = 2.54 * cm
    foot = 12. * inch
    mile = 5280. * foot
    thou = 1e-3 * inch  # thousandth of an inch; also called mil

    # Volume
    global L, mL, uL, nL, pL, fL, aL, kL, ML, GL
    L = 1e-3 * m**3 # liter
    mL = 1e-3 * L
    uL = 1e-6 * L
    nL = 1e-9 * L
    pL = 1e-12 * L
    fL = 1e-15 * L
    aL = 1e-18 * L
    kL = 1e3 * L
    ML = 1e6 * L
    GL = 1e9 * L

    # Time
    global ms, us, ns, ps, fs, minute, hour, day, week, year
    ms = 1e-3 * s
    us = 1e-6 * s
    ns = 1e-9 * s
    ps = 1e-12 * s
    fs = 1e-15 * s
    minute = 60. * s
    hour = 60. * minute
    day = 24. * hour # solar day
    week = 7. * day
    year = 365.256363004 * day # sidereal year

    # Frequency
    global Hz, mHz, kHz, MHz, GHz, THz, PHz, rtHz, rpm
    Hz = 1./s
    mHz = 1e-3 * Hz
    kHz = 1e3 * Hz
    MHz = 1e6 * Hz
    GHz = 1e9 * Hz
    THz = 1e12 * Hz
    PHz = 1e15 * Hz
    rtHz = Hz**0.5 # "root Hertz"
    rpm = 1/minute # revolutions per minute

    # Angular frequency
    # Example: ω = 3 * kHz·2π means that ω is the angular frequency
    #   corresponding to a rotation whose *ordinary* frequency is 3 kHz.
    global Hz·2π, mHz·2π, kHz·2π, MHz·2π, GHz·2π, THz·2π, PHz·2π, rpm·2π
    Hz·2π  =  Hz * 2*pi
    mHz·2π = mHz * 2*pi
    kHz·2π = kHz * 2*pi
    MHz·2π = MHz * 2*pi
    GHz·2π = GHz * 2*pi
    THz·2π = THz * 2*pi
    PHz·2π = PHz * 2*pi
    rpm·2π = rpm * 2*pi

    # Mass
    global g, mg, ug, ng, pg, fg, tonne, amu, Da, kDa, lbm
    g = 1e-3 * kg
    mg = 1e-3 * g
    ug = 1e-6 * g
    ng = 1e-9 * g
    pg = 1e-12 * g
    fg = 1e-15 * g
    tonne = 1e3 * kg
    amu = 1.66053906660e-27 * kg # atomic mass unit
    Da = amu # Dalton
    kDa = 1e3 * Da
    lbm = 0.45359237 * kg # pound mass (international avoirdupois pound)

    # Energy
    global J, mJ, uJ, nJ, pJ, fJ, kJ, MJ, GJ, erg, eV, meV, keV, MeV, GeV, \
           TeV, btu, smallcal, kcal, Wh, kWh
    J = (kg * m**2)/s**2
    mJ = 1e-3 * J
    uJ = 1e-6 * J
    nJ = 1e-9 * J
    pJ = 1e-12 * J
    fJ = 1e-15 * J
    kJ = 1e3 * J
    MJ = 1e6 * J
    GJ = 1e9 * J
    erg = 1e-7 * J
    eV = 1.602176634e-19 * J
    meV = 1e-3 * eV
    keV = 1e3 * eV
    MeV = 1e6 * eV
    GeV = 1e9 * eV
    TeV = 1e12 * eV
    btu = 1055.06 * J  # British thermal unit
    smallcal = 4.184 * J # small calorie ("gram calorie")
    kcal = 4184. * J  # kilocalorie ("large Calorie", "dietary Calorie")
    Wh = 3600. * J # watt-hour
    kWh = 1e3 * Wh # kilowatt-hour

    # Moles, concentration / molarity
    global NA, mol, mmol, umol, nmol, pmol, fmol, M, mM, uM, nM, pM, fM
    NA = 6.02214076e23  # Avogadro's number
    mol = NA  #1 mole (see README)
    mmol = 1e-3 * mol
    umol = 1e-6 * mol
    nmol = 1e-9 * mol
    pmol = 1e-12 * mol
    fmol = 1e-15 * mol
    M = mol/L # molar
    mM = 1e-3 * M
    uM = 1e-6 * M
    nM = 1e-9 * M
    pM = 1e-12 * M
    fM = 1e-15 * M

    # Force
    global N, mN, uN, nN, pN, fN, kN, MN, GN, dyn, lbf
    N = (kg * m)/s**2 # newton
    mN = 1e-3 * N
    uN = 1e-6 * N
    nN = 1e-9 * N
    pN = 1e-12 * N
    fN = 1e-15 * N
    kN = 1e3 * N
    MN = 1e6 * N
    GN = 1e9 * N
    dyn = 1e-5 * N # dyne
    lbf = lbm * (9.80665 * m/s**2) # pound-force (international avoirdupois pound)

    # Pressure
    global Pa, hPa, kPa, MPa, GPa, bar, mbar, cbar, dbar, kbar, Mbar, atm, \
           torr, mtorr, psi
    Pa = N/m**2 # pascal
    hPa = 1e2 * Pa # hectopascal
    kPa = 1e3 * Pa
    MPa = 1e6 * Pa
    GPa = 1e9 * Pa
    bar = 1e5 * Pa
    mbar = 1e-3 * bar
    cbar = 1e-2 * bar # centibar
    dbar = 0.1 * bar # decibar
    kbar = 1e3 * bar
    Mbar = 1e6 * bar
    atm = 101325. * Pa
    torr = (1./760.) * atm
    mtorr = 1e-3 * torr
    psi = lbf / inch**2

    # Power
    global W, mW, uW, nW, pW, kW, MW, GW, TW, \
           horsepower_imperial, horsepower_metric
    W = J/s
    mW = 1e-3 * W
    uW = 1e-6 * W
    nW = 1e-9 * W
    pW = 1e-12 * W
    kW = 1e3 * W
    MW = 1e6 * W
    GW = 1e9 * W
    TW = 1e12 * W
    horsepower_imperial = 33000 * foot * lbf / minute
    horsepower_metric = (75 * kg) * (9.80665 * m/s**2) * (1 * m/s)

    # Acceleration and related
    global Gal, mGal, uGal, eotvos
    Gal = 1*cm/s**2
    mGal = 1e-3 * Gal
    uGal = 1e-6 * Gal
    eotvos = 1e-9 / s**2

    # Temperature
    global degFinterval, degCinterval, mK, uK, nK, pK
    degFinterval = (5./9.) * K # A temperature difference in degrees Fahrenheit
    degCinterval = K # A temperature difference in degrees Celsius
    mK = 1e-3 * K
    uK = 1e-6 * K
    nK = 1e-9 * K
    pK = 1e-12 * K

    # Charge
    global mC, uC, nC, Ah, mAh
    mC = 1e-3 * C
    uC = 1e-6 * C
    nC = 1e-9 * C
    Ah = 3600. * C # amp-hour
    mAh = 1e-3 * Ah

    # Current
    global A, mA, uA, nA, pA, fA
    A = C/s
    mA = 1e-3 * A
    uA = 1e-6 * A
    nA = 1e-9 * A
    pA = 1e-12 * A
    fA = 1e-15 * A

    # Voltage
    global V, mV, uV, nV, kV, MV, GV, TV
    V = J/C
    mV = 1e-3 * V
    uV = 1e-6 * V
    nV = 1e-9 * V
    kV = 1e3 * V
    MV = 1e6 * V
    GV = 1e9 * V
    TV = 1e12 * V

    # Resistance and conductivity
    global ohm, mohm, kohm, Mohm, Gohm, Ω, mΩ, kΩ, MΩ, GΩ, S, mS, uS, nS
    ohm = V / A
    mohm = 1e-3 * ohm
    kohm = 1e3 * ohm
    Mohm = 1e6 * ohm
    Gohm = 1e9 * ohm
    Ω = ohm # easier-to-read alias (see https://sjbyrnes.com/unicode.html )
    mΩ = mohm # easier-to-read alias
    kΩ = kohm # easier-to-read alias
    MΩ = Mohm # easier-to-read alias
    GΩ = Gohm # easier-to-read alias
    S = 1./ohm # siemens
    mS = 1e-3 * S
    uS = 1e-6 * S
    nS = 1e-9 * S

    # Magnetic fields and fluxes
    global T, mT, uT, nT, G, mG, uG, kG, Oe, Wb
    T = (V * s) / m**2 # tesla
    mT = 1e-3 * T
    uT = 1e-6 * T
    nT = 1e-9 * T
    G = 1e-4 * T # gauss
    mG = 1e-3 * G
    uG = 1e-6 * G
    kG = 1e3 * G
    Oe = (1000./(4.*pi)) * A/m # oersted
    Wb = J/A # weber

    # Capacitance and inductance
    global F, uF, nF, pF, fF, aF, H, mH, uH, nH
    F = C / V # farad
    uF = 1e-6 * F
    nF = 1e-9 * F
    pF = 1e-12 * F
    fF = 1e-15 * F
    aF = 1e-18 * F
    H = m**2 * kg / C**2 # henry
    mH = 1e-3 * H
    uH = 1e-6 * H
    nH = 1e-9 * H

    # Constants--general
    global c0, mu0, μ0, eps0, ε0, Z0, hPlanck, hbar, ħ, kB, GNewton, sigmaSB, σSB, alphaFS, αFS
    c0 = 299792458. * m/s  # speed of light in vacuum
    mu0 = 1.25663706212e-6 * N/A**2  # magnetic constant, permeability of vacuum
    μ0 = mu0 # easier-to-read alias (see https://sjbyrnes.com/unicode.html )
    eps0 = 1./(mu0 * c0**2) # electric constant, permittivity of vacuum
    ε0 = eps0 # easier-to-read alias
    Z0 = mu0 * c0  # vacuum impedance, 377 ohms
    hPlanck = 6.62607015e-34 * J*s  # planck constant
    hbar = hPlanck / (2.*pi)  # reduced planck constant
    ħ = hbar # easier-to-read alias
    kB = 1.380649e-23 * J/K  # Boltzmann constant
    GNewton = 6.67430e-11 * m**3 / (kg * s**2) # Gravitational constant
    sigmaSB = (pi**2 / 60.) * kB**4 / (hbar**3 * c0**2)  # Stefan-Boltzmann constant
    σSB = sigmaSB # easier-to-read alias
    alphaFS = 7.2973525693e-3  # fine-structure constant
    αFS = alphaFS  # easier-to-read alias

    # Constants--chemistry, atomic physics, electrons
    global Rgas, e, uBohr, uNuc, aBohr, me, mp, mn, Rinf, Ry, Hartree, \
           ARichardson, Phi0, KJos, RKlitz
    Rgas = kB # ideal gas constant (see README)
    e = 1.602176634e-19 * C  # charge of proton
    uBohr = 9.2740100783e-24 * J/T  # Bohr magneton
    uNuc = 5.0507837461e-27 * J/T # nuclear magneton
    aBohr = 5.29177210903e-11 * m  # Bohr radius
    me = 9.1093837015e-31 * kg  # electron mass
    mp = 1.67262192369e-27 * kg # proton mass
    mn = 1.67492749804e-27 * kg # neutron mass
    Rinf = 10973731.568160 / m # Rydberg constant
    Ry = 2.1798723611035e-18 * J  # Rydberg energy, approximately 13.6 eV
    Hartree = 2*Ry # Hartree energy, approximately 27.2 eV
    ARichardson = (4.*pi*e*me*kB**2) / hPlanck**3  # Richardson constant
    Phi0 = hPlanck / (2*e) # magnetic flux quantum
    KJos = (2*e) / hPlanck # Josephson constant
    RKlitz = hPlanck / e**2 # von Klitzing constant

    # Constants--astronomical and properties of earth
    global REarth, g0, Msolar, MEarth
    REarth = 6371. * km # radius of earth
    g0 = 9.80665 * m / s**2 # standard earth gravitational acceleration
    Msolar = 1.98847e30 * kg # mass of the sun
    MEarth = 5.9722e24 * kg # mass of earth

# Set units randomly when this module is initialized. (Don't worry: If the
# module is imported many times from many places, this command will only
# execute during the first import.)
reset_units()


def nu_eval(expression):
    """
    Evaluates a string expression in the context of this module, so that you
    can make APIs that don't require their users to import numericalunits;
    instead the API user can run a function like load_data(data, unit='km')

    For example:

        import numericalunits as nu
        x = nu.nu_eval('kg * m / s**2')

    ...is exactly equivalent to...

        import numericalunits as nu
        x = nu.kg * nu.m / nu.s**2

    Input strings are required to be of the form of stereotypical unit
    expressions—e.g. addition and subtraction are banned—to catch user errors.
    """
    # Based on https://stackoverflow.com/a/9558001
    import ast
    import operator as op
    operators = {ast.Mult: op.mul, ast.Div: op.truediv, ast.Pow: op.pow, ast.USub: op.neg}

    def _eval(node):
        if isinstance(node, ast.Num):
            return node.n
        elif isinstance(node, ast.BinOp):  # <left> <operator> <right>
            return operators[type(node.op)](_eval(node.left), _eval(node.right))
        elif isinstance(node, ast.UnaryOp):  # <operator> <operand> e.g., -1
            return operators[type(node.op)](_eval(node.operand))
        elif isinstance(node, ast.Name):
            return globals()[node.id]
        else:
            raise TypeError(node)

    return _eval(ast.parse(expression, mode='eval').body)
