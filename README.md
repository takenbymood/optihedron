# optihedron
An engine of algorithms for design optimisation of membrane reshaping nanoparticles.

The preprint of the study is found at the following address:

[Designing membrane-reshaping nanostructures through artificial evolution](https://www.biorxiv.org/content/10.1101/2020.02.27.968149v1)

This forked version of the larger (and messier) original repo is the one which is linked for paper submission, it includes everything needed to install and run the main codebase but does not contain the data analysis scripts used to generate the figures from the paper. These files are found in the following repository:

[optihedron-analysis](https://github.com/takenbymood/optihedron-analysis)

The following commands will download the Feb2016 lammps release and make install it along with the Yuan membrane potential used in the study and connect it to python. The install script will also get all required pip packages from requirements.txt and install them in its own virtual environment, activate that environment and run a genetic algorithm on a 22 ligand particle with 11kT affinity. The output will be saved to db/datastore.db as a sqlite file, which is required for analysis of the results.

```bash
source ./install.sh
source ./activate.sh
python run.py -n 10 -p 50 -d 1 -
```
You will need wget, python2.7, pip and virtualenv to be installed for all of this to work. Please ensure you do that before running install.sh

The main script is unimaginatively titled "run.py" and has the following flags as options

#### Genetic Algorithm Options

| Flag        | Default           | Description  |
| ------------- |:-------------| :-----|
| -n (--ngen) | 1 | The number of generations to run for before stopping |
| -d (--demes) | 1 | The number of demes/subpopulations |
| -p (--population) | 1 | The number of individuals in each deme |
| -f (--migfreq) | 1 | The number of generations between each migration event |
| -mg (--migrations) | 1 | The number of migrations to do each migfreq generations |
| -c (--cxpb) | 0.5 | The independent probability that an individual is selected for crossover |
| -m (--mutpb) | 0.2 | The independent probability that an individual is selected for mutation |
| -mpb (--mindpb) | 0.05 | The independent probability that a bit in the genome is chosen for mutation |
| -t (--tournsize) | 3 | The number of individuals in each selection tournament |
| -g (--graph) | islands | The topology of the migration network, options: singlet, islands, star, megastar |
| -mo (--mutation) | defaultMut | The mutation regime to use, options: defaultMut, fixedActivationMut |
| -xo (--mate) | defaultGeneWiseTwoPoint | The crossover regime to use, options: defaultGeneWiseTwoPoint, fixedActivationGeneWiseTwoPoint |
| -sg (--startinggen) | 0 | The generation to start at (for continuing checkpointed runs) |
| -ff (--fitnessfunction) | budding | The fitness function to use for assigning scores, options: budding, diameter, smallworld |
| -a (--algorithm) | eaSimple | The algorithm to apply during the run, options: eaSimple |
| -gs (--genomesize) | 40 | The number of genes in each individual genome (can be overwritten) |

#### Model Options

| Flag        | Default           | Description  |
| ------------- |:-------------| :-----|
| -s (--seed) | ```python int(time.time())``` | The number of generations to run for before stopping |

#### Concurrency Options

| Flag        | Default           | Description  |
| ------------- |:-------------| :-----|
| -mpi (--mpi) | False | Run LAMMPs using MPI, this is highly reccomended if you are running only 1 worker |

#### Data Options

| Flag        | Default           | Description  |
| ------------- |:-------------| :-----|
| -v (--verbose) | False | Flag to run the code in verbose mode (more printing) |
