

========================================================================================
Bioelectrics: Monodomain problem on a 2D domain with simple HodgkinHuxley1952 cell model
========================================================================================
This example demonstrates solving a monodomain problem on a 2D domain using Hodgkin and Huxley cell model passed as a .cellml model. The geometry is discretised with bilinear Lagrange finite elements with numberOfXElements in the X direction and numberOfYElements in the Y direction. The material parameters for the domain are controlled with the Am parameter for the cell membrane area, Cm for the membrane capacitance and conductivity for continuum conducticity.

To start the simulation a stimulation is applied to the left half of the bottom row of nodes. This stimulation lasts from time zero until time stimStop. The magnitude of the stimulus is given by stimValue parameter. After time stimStop the stimulation is turned of and the simulation is continued until time timeStop. The time step for the spatial PDE problem is given by the pdeTimeStep parameter and the time step for the ODE integration is given by the odeTimeStep parameter. The final parameter, outputFrequency, controls how many time steps pass before the solution is output to file.

Running the example
===================
Activate the OpenCMISS virtual environment for python bindings::

  source <path-to-opencmiss>/install/virtual_environments/oclibs_venv_py27_release/bin/activate

Run the python script::

  cd ./src/Python/
  python monodomain_2D_HH.py

To run in parallel::

  cd ./src/Python/
  mpiexec -np <number-of-processors> python monodomain_2D_HH.py

Verifying the example
=====================

Results can be visualised by running `visualise.cmgui <./src/Python/visualise.cmgui>`_ with the `Cmgui visualiser <http://physiomeproject.org/software/opencmiss/cmgui/download>`_.
The expected results from this example are available in `expected_results <./src/Python/expected_results>`_ folder.

Prerequisites
=============
There are no additional input files required for this example as it is self-contained.

License
=======
License applicable to this example is described in `LICENSE <./LICENSE>`_.
