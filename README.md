# Correction Vector Method
The [correction vector method](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.60.335) allows to iterativly calculate the spectral function of an aribitrary Hamiltonian. <br/> We implement the algorithm using [ITensor](https://itensor.org/).

## Prerequisites
To install [ITensor](https://itensor.org/) please follow the instructions on the [ITensor GitHub repository](https://github.com/ITensor/ITensor).

## Running the example
To calculate the spectral function of the Heisenberg model:

* Run ```make``` to compile the program.<br/>
* Generate a plot of the spectral function by running ```python main.py``` followed by ```python plot.py```.

![spectral_func](https://user-images.githubusercontent.com/45107198/61218087-1fb4d600-a709-11e9-85e8-4a102e9b1c34.png)
