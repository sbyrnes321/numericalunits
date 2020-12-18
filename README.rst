=========================================================================
numericalunits: Units and dimensional analysis compatible with everything
=========================================================================

`Package homepage at PyPI <http://pypi.python.org/pypi/numericalunits>`_ -- 
`Source code at github <http://github.com/sbyrnes321/numericalunits>`_ -- 
Written by `Steve Byrnes <http://sjbyrnes.com>`_

This package implements units and dimensional analysis in an unconventional 
way that has the following unique advantages:

* **Compatible with everything:** Compatible with virtually any numerical 
  calculation routine, including numpy and scipy, and even including routines 
  not written in Python! That means, for example, if you have a decades-old 
  closed-source C routine for numerical integration, you can pass it a 
  quantity with units of velocity and an integration range with units of 
  time, and the final answer will magically have units of distance. This 
  extreme compatibility is possible because if the variable ``x`` represents 
  a quantity with dimensions (like "3.3 kg"), ``x`` is actually stored 
  internally as an ordinary floating-point number. The dimension is 
  encoded in the value as a multiplicative factor. When two numbers are 
  multiplied, their dimensions are automatically multiplied, and so on. 


* **Modular and non-intrusive:** When you input data, you say what units 
  they are in. When you display results, you say what units you want to 
  display them in. These steps are very little trouble, and in fact help you 
  create nice, self-documenting code. Other than that, you have to do nothing 
  at all to pass dimensionful quantities into and out of any already-written 
  programs or routines.

* **Powerful tool for debugging:** Not *all* calculation mistakes cause 
  violations of dimensional analysis, but *most* do--for example, if you 
  accidentally multiply two lengths instead of adding them, the result will 
  have the wrong dimension. If you use this package, it will alert you to 
  these sorts of mistakes.

* **Zero storage overhead**

* **Zero calculation overhead**

These great features come with the disadvantage that the interface is less  
*slick* than other unit packages. If you have a quantity with units, you 
cannot directly see what the units are. You are supposed to already know 
what the units are, and then the package will tell you whether you made a 
mistake. Even worse, you only get alerted to the mistake after running a 
calculation all the way through twice.

Therefore the package is *not* suggested for students exploring how units work.
It *is* suggested for engineering and science professionals who want to make
their code more self-documenting and self-debugging.

Installation
============

You can install from PyPI: ::

    pip install numericalunits

Alternatively---since it's a single module that requires no setup or 
compilation---you can download ``numericalunits.py`` from `PyPI 
<http://pypi.python.org/pypi/numericalunits>`_ or `github 
<http://github.com/sbyrnes321/numericalunits>`_ and use it directly.

Usage and examples
==================

To assign a unit to a quantity, **multiply** by the unit, e.g.
``my_length = 100 * mm``. (In normal text you would write "100 mm", but
unfortunately Python does not have "implied multiplication".)

To express a dimensionful quantity in a certain unit, **divide** by that unit,
e.g. when you see ``my_length / cm``, you pronounce it "my_length expressed
in cm".

Unit errors, like trying to add a length to a mass, will not *immediately*
announce themselves as unit errors. Instead, you need to run the whole
calculation twice (in a new Python session each time). If you get the
same final answers both times, then congratulations, all your calculations
are almost guaranteed to pass dimensional analysis! If you get different
answers every time you run, then you made a unit error! It is up to you to
figure out where and what the error is.

**Example 1:** What is 5 mL expressed in cubic nanometers?::

    from numericalunits import mL, nm
    x = 5 * mL  # "Read: x is equal to 5 milliliters"
    print(x / nm**3)   # "Read: x expressed in cubic nanometers is..." --> 5e21

**Example 2:** An electron is in a 1e5 V/cm electric field. What is its
acceleration? (Express the answer in m/s².) ::

    from numericalunits import V, cm, e, me, m, s
    Efield = 1e5 * (V / cm)
    force = e * Efield # (e is the elementary charge)
    accel = force / me # (me is the electron mass)
    print(accel / (m / s**2)) # Answer --> 1.7588e18

**Example 3:** You measured a voltage as a function of the position of dial: 
10 volts when the dial is at 1cm, 11 volts when the dial is at 2cm, etc. 
etc. Interpolate from this data to get the expected voltage when the dial is 
at 41mm, and express the answer in mV. ::

    from numericalunits import cm, V, mm, mV
    from numpy import array
    from scipy.interpolate import interp1d
    voltage_data = array([[1 * cm, 10 * V],
                          [2 * cm, 11 * V],
                          [3 * cm, 13 * V],
                          [4 * cm, 16 * V],
                          [5 * cm, 18 * V]])
    f = interp1d(voltage_data[:,0], voltage_data[:,1])
    print(f(41 * mm) / mV) # Answer --> 16200
	

**Example 4:** A unit mistake ... what is 1 cm expressed in atmospheres? ::

    from numericalunits import cm, atm
    print((1 * cm) / atm) # --> a randomly-varying number
    # The answer randomly varies every time you run this (in a new Python
    # session), indicating that you are violating dimensional analysis.

Parallelism
-----------

Whenever you spawn multiple independent workers solving a common task you have
to ensure that the units are the same. You have three possible options here:

* Initialize numericalunits with a common seed. ::

      # multiprocessing is given for the example
      from numericalunits import reset_units
      from multiprocessing import Pool

      reset_units(42)
      pool = multiprocessing.Pool(initializer=reset_units, initargs=(42,))

  The above code ensures that all dependent processes in the Pool together with
  are master process are initialized with the same seed 42. It works only if
  you have control of the first import of numericalunits.

* A more robust way is to serialize the full content of the numericalunits package. ::

      import numericalunits
      from multiprocessing import Pool

      def serialize_nu():
          return {k: v for k, v in numericalunits.__dict__.items() if isinstance(v, float) and not k.startswith("_")}
      def load_nu(data):
          numericalunits.__dict__.update(data)

      # ... some code here ...

      state = serialize_nu()
      pool = multiprocessing.Pool(initializer=load_nu, initargs=(state,))

* The final option is to ensure a proper conversion. I.e. the data passed to
  and returned from workers has specific units applied. numericalunits can
  still be used inside workers. ::

      from numericalunits import eV, J, kB, K
      from multiprocessing import Pool
      
      energies = eV, 2 * eV, 3 * eV

      def worker(energy):
          # Lift units
          energy *= J

          # Compute something
          result = energy / kB

          # Apply proper units before returning
          return result / K
      
      pool = multiprocessing.Pool()
      temperatures = pool.map(worker, (i / J for i in energies))  # Apply units before passing to worker ...
      temperatures = tuple(i * K for i in temperatures)  # ... and lift units after obtaining the result

How it works
============

A complete set of independent base units (meters, kilograms, seconds, 
coulombs, kelvins) are defined as randomly-chosen positive floating-point 
numbers. All other units and constants are defined in terms of those. In a 
dimensionally-correct calculation, the units all cancel out, so the final 
answer is deterministic, not random. In a dimensionally-incorrect 
calculations, there will be random factors causing a randomly-varying final 
answer.

Included units and constants
============================

Includes a variety of common units, both SI and non-SI, everything from 
frequency to magnetic flux. Also includes common physical constants like 
Planck's constant and the speed of light. Browse the source code to see a 
complete list. It is very easy to add in any extra units and constants that
were left out.

Notes
=====

Notes on implementation and use
-------------------------------

* **What does it mean to "run the calculation again in a new Python
  session?"** You know that you've started a new Python session if all
  the variable definitions have been forgotten. Three examples: In Spyder, each "Console"
  tab is its own session. In Jupyter, make a new Python session by selecting
  "Restart kernel". From the command line, each time you type
  ``python blah.py``, you are opening a new Python session. 

* For little, self-contained calculations (a few lines that are all within a
  single module), it is possible to check the units without opening a new Python
  session: Run the function ``numericalunits.reset_units()`` at the beginning of
  the calculation before any variables are defined; then check for
  dimensional errors by re-running the whole calculation (including the
  ``reset_units()`` part). Note that if you are using ``from``-style imports,
  like ``from numericalunits import cm``, you need to put them *after*
  ``reset_units()`` in the code.

* While debugging a program, it may be annoying to have intermediate values 
  in the calculation that randomly vary every time you run the program. In 
  this case, you can use ``reset_units('SI')`` instead of the normal 
  ``reset_units()``. This puts all dimensionful variables in standard (MKS)
  SI units: All times are in seconds, all lengths are in meters, all forces
  are in newtons, etc. Alternatively, ``reset_units(123)`` uses ``123`` as
  the seed for the random-number generator. Obviously, in these modes, you
  will *not* get any indication of dimensional-analysis errors. As above,
  if you are going to use any version of ``reset_units()``, make sure you do
  it before any dimensionful variable is defined in any module.

* There are very rare, strange cases where the final answer does not seem to 
  randomly vary even though there was a dimensional-analysis violation: For 
  example, the expression ``(1 + 1e-50 * cm / atm)`` fails dimensional 
  analysis, so if you calculate it the answer is randomly-varying. But, it is 
  only randomly varying around the 50th decimal point, so the variation is
  hidden from view. You would not notice it as an error.

* Since units are normal Python ``float``-type numbers, they follow the normal
  casting rules. For example, ``2 * cm`` is a python ``float``, not an ``int``.
  This is usually what you would want and expect.

* You can give a dimension to complex numbers in the same way as real 
  numbers--for example ``(2.1e3 + 3.9e4j) * ohm``.

* Requires Python 3. (For Python 2 compatibility, install numericalunits
  version 1.23 or earlier.)

* If you find bugs, please tell me by `email <http://sjbyrnes.com>`_ or 
  `github issue board <https://github.com/sbyrnes321/numericalunits/issues>`_.

* If you get overflows or underflows, you can edit the unit initializations.
  For example, the package sets the meter to a random numerical value between 0.1
  and 10. Therefore, if you're doing molecular simulation, most lengths you
  use will be tiny numbers. You should probably set the meter instead to be
  between, say, a random numerical value between 1e8 and 1e10.

* Some numerical routines use a default *absolute* tolerance, rather than
  relative tolerance, to decide convergence. This can cause the calculation
  result to randomly vary even though there is no dimensional analysis error.
  When this happens, you should set the absolute tolerance to a value with the
  appropriate units. Alternatively, you can scale the data before running the
  algorithm and scale it back afterwards. Maybe this sounds like a hassle, but
  it's actually a benefit: If your final result is very sensitive to some
  numerical tolerance setting, then you really want to be aware of that.

Notes on unit definitions
-------------------------

* For electromagnetism, all units are intended for use in SI formulas. If 
  you plug them into cgs-gaussian electromagnetism formulas, or cgs-esu 
  electromagnetism formulas, etc., you will get nonsense results.

* The package does not keep track of "radians" as an independent unit 
  assigned a random number. The reason is that the "radians" factor does not 
  always neatly cancel out of formulas.

* The package does not keep track of "moles" as an independent unit assigned 
  a random number; instead ``mol`` is just a pure number (~6e23), like you
  would say "dozen"=12. That means: (1) gram/mol is exactly the same as amu,
  and Boltzmann constant is exactly the same as the ideal gas constant, and so
  on. (2) You should rarely need to use Avogadro's number ``NA`` -- it is just a
  synonym of ``mol`` (``NA = mol ~ 6e23``). Here are a few examples using moles: ::
  
      from numericalunits import um, uM, kcal, mol, fmol, J
      
      # There are eight copies of a protein inside a yeast nucleus of volume
      # 3 cubic microns. What is the concentration of the protein, in micromolar (uM)?
      print((8 / (3 * um**3)) / uM)   # Answer --> 0.0044
      
      # 5 kcal / mol is how many joules?
      print((5 * kcal / mol) / J)   # Answer --> 3.47e-20
      
      # How many molecules are in 2.3 femtomoles?
      print(2.3 * fmol)   # Answer --> 1.39e9

* The package cannot convert temperatures between Fahrenheit, Celsius, and 
  kelvin. The reason is that these scales have different zeros, so the units 
  cannot be treated as multiplicative factors. It is, however, possible to 
  convert temperature *intervals*, via the units ``degCinterval`` (which is a 
  synonym of kelvin, ``K``) and ``degFinterval``.
