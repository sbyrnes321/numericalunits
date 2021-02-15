=========================================================================
numericalunits README appendix
=========================================================================

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
