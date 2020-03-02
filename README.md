# optihedron
An engine of algorithms for design optimisation of membrane reshaping nanoparticles.

The preprint of the study is found at the following address:

[Designing membrane-reshaping nanostructures through artificial evolution](https://www.biorxiv.org/content/10.1101/2020.02.27.968149v1)

This forked version of the larger (and messier) original repo is the one which is linked for paper submission, it includes everything needed to install and run the main codebase but does not contain the data analysis scripts used to generate the figures from the paper. These files are found in the following repository:

[optihedron-analysis](https://github.com/takenbymood/optihedron-analysis)

The following commands will download the Feb2016 lammps release and make install it along with the Yuan membrane potential used in the study and connect it to python. The install script will also get all required pip packages from requirements.txt and install them in its own virtual environment, activate that environment and run a genetic algorithm on a 22 ligand particle with 11kT affinity. The output will be saved to db/datastore.db as a sqlite file, which is required for analysis of the results. It will also save the best individuals from each generation as lammps files along with their trajectories.

```bash
source ./install.sh
source ./activate.sh
python run.py -n 10 -p 50 -d 1 -fl 22 -epmn 11 -epmx 11 -br 100 -tw 10 -sr -kb -w 4
```
You will need wget, python2.7, pip and virtualenv to be installed for all of this to work. Please ensure you do that before running install.sh

## Options

The main script is unimaginatively titled "run.py" and has the following flags as options. Some of these have default values which do not make much sense for the fixed ligand number model described in the paper. If you are wanting to reproduce those results be sure to set -fl, -epmn, and -epmx accordingly.

#### Genetic Algorithm Options

| Flag | Long | Default | Description  |
| ----- | :-------- | :--- | :- |
| -n | --ngen | 1 | The number of generations to run for before stopping |
| -d | --demes | 1 | The number of demes/subpopulations |
| -p | --population | 1 | The number of individuals in each deme |
| -f | --migfreq | 1 | The number of generations between each migration event |
| -mg | --migrations | 1 | The number of migrations to do each migfreq generations |
| -c | --cxpb | 0.5 | The independent probability that an individual is selected for crossover |
| -m | --mutpb | 0.2 | The independent probability that an individual is selected for mutation |
| -mpb | --mindpb | 0.05 | The independent probability that a bit in the genome is chosen for mutation |
| -t | --tournsize | 3 | The number of individuals in each selection tournament |
| -g | --graph | islands | The topology of the migration network, options: singlet, islands, star, megastar |
| -mo | --mutation | defaultMut | The mutation regime to use, options: defaultMut, fixedActivationMut |
| -xo | --mate | defaultGeneWiseTwoPoint | The crossover regime to use, options: defaultGeneWiseTwoPoint, fixedActivationGeneWiseTwoPoint |
| -sg | --startinggen | 0 | The generation to start at (for continuing checkpointed runs) |
| -ff | --fitnessfunction | budding | The fitness function to use for assigning scores, options: budding, diameter, smallworld |
| -a | --algorithm | eaSimple | The algorithm to apply during the run, options: eaSimple |
| -gs | --genomesize | 40 | The number of genes in each individual genome (can be overwritten) |

#### Model Options

| Flag | Long | Default | Description  |
| ----- | :-------- | :--- | :- |
| -s | --seed | ```python int(time.time())``` | The number of generations to run for before stopping |
| -hof | --hofsize | 5 | The size of the hall of fame (saves best individuals) |
| -expr | --exprplaces | 1 | The number of bits in each gene for a ligand to switched on or off (off if all bits=0) |
| -eps | --epsplaces | 8 | Number of bits in each gene to encode the value of epsilon (affinity) |
| -epmn | --epsmin | 0 | Minimum value for epsilon for each ligand (affinity) |
| -epmx | --epsmax | 15 | Maximum value for epsilon for each ligand (affinity) |
| -r | --runtime | 25000 | Number of timesteps to run the LAMMPs simulation for |
| -ts | --timestep | 0.01 | LAMMPs timestep size |
| -rs | --repeats | 4 | Number of repeat simulations to run for each unique individual |
| -br | --buddingreward | 400 | Reward for crossing the membrane |
| -tw | --timeweight | 25 | Reward for fast budding (mutiplied by runtime/tb) |
| -fl | --fixedligands | -1 | Number of ligands allowed on the particle (-1 allows any number) |
| -pp | --partialpacking | False | Runs the algorithm in partial packing mode which allows free placement of ligands |
| -polang | --polangplaces | 8 | Number of bits in each gene to encode polar position of a ligand (in pp mode) |
| -aziang | --aziangplaces | 8 | Number of bits in each gene to encode azimuthal position of a ligand (in pp mode) |


#### Concurrency Options

| Flag | Long | Default | Description  |
| ----- | :-------- | :--- | :- |
| -mpi | --mpi | False | Run LAMMPs using MPI, this is highly reccomended if you are running only 1 worker |
| -w | --workers | 10 | Number of parallel workers to run (probably shouldn't exceed the number of cores you have) |
| -np | --nodes | 4 | Number of nodes to run each MPI task on |
| -tm | --timeout | 1800 | How long to wait before considering an mpi run broken |

#### Data Options

| Flag | Long | Default | Description  |
| ----- | :-------- | :--- | :- |
| -v | --verbose | False | Flag to run the code in verbose mode (more printing) |
| -sr | --saveresults | False | Save results in a database |
| -ko | --keepoutput | False | Keep all .xyz and .xyza trajectories |
| -ki | --keepinput | False | Keep all .in and .data LAMMPs files |
| -kb | --keepbest | False | Keep the input and outputs of the best individual each generation |
| -kw | --keepworst | False | Keep the input and outputs of the worst individual each generation |
| -wd | --wdir | ```python os.path.dirname(os.path.realpath(__file__))``` | Sets the working directory |
| -i | --input | None | Read parameters from a json file |
| -db | --database | None | Path to the database file to read state from (for checkpointing) |
| -dbs | --dbsession | -1 | Set the id of the session to reload from the specified database |
| -z | --zoo | None | Set the input to a zoo file for benchmarking performance |

#### Deprecated Options

| Flag | Long | Default | Description  |
| ----- | :-------- | :--- | :- |
| -lw | --ligandweight | 10 | Reward for being at the target ligand number |
| -pw | --penaltyweight | 10 | Penalty for higher than target affinity |
| -tl | --targetligands | 1 | Ideal number of ligands |
| -q | --qsub | False | Submit jobs to a cluster using qsub |



