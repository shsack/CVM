# Correction Vector Method
The [correction vector method](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.60.335) allows to iterativly calculate the spectral function of an aribitrary Hamiltonian. <br/> We implement the algorithm using [ITensor](https://itensor.org/).

## Prerequisites
To install [ITensor](https://itensor.org/) please follow the instructions on the [ITensor GitHub repository](https://github.com/ITensor/ITensor).

## Running the program
To calculate the spectral function of the 1D Heisenberg model:

![heisenberg](https://user-images.githubusercontent.com/45107198/61218669-6ce57780-a70a-11e9-9abd-a148c7538c66.png)

* Run ```make``` to compile the program.<br/>
* Generate a plot of the spectral function by running ```python main.py``` followed by ```python plot.py```.

![spectral_func](https://user-images.githubusercontent.com/45107198/61218087-1fb4d600-a709-11e9-85e8-4a102e9b1c34.png)

The chain length as well as other parameters can be adjusted in the ```main.py``` file.

```run(omega=o, N=5, eta=0.05, max_iter=50, tol=1E-5, i=1, j=1, maxm=20, cut=1E-6)
```
