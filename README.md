# Optihedron
An engine of algorithms for design optimisation of membrane reshaping nanoparticles.

The preprint of the study is found at the following address:

[Designing membrane-reshaping nanostructures through artificial evolution](https://www.biorxiv.org/content/10.1101/2020.02.27.968149v1)

This forked version of the larger (and messier) original repo is the one which is linked for paper submission, it includes everything needed to install and run the main codebase but does not contain the data analysis scripts used to generate the figures from the paper. These files are found in the following repository:

[optihedron-analysis](https://github.com/takenbymood/optihedron-analysis)

## Installing and Running

The following commands will download the Feb2016 lammps release and make install it along with the Yuan membrane potential used in the study and connect it to python. The install script will also get all required pip packages from requirements.txt and install them in its own virtual environment, activate that environment and run a genetic algorithm on a 22 ligand particle with 11kT affinity. The output will be saved to db/datastore.db as a sqlite file, which is required for analysis of the results. It will also save the best individuals from each generation as lammps files along with their trajectories.

```bash
source ./install.sh
source ./activate.sh
python run.py -w 4 -n 30 -p 16 -d 1 -fl 22 -expr 1 -eps 0 -polang 0 -aziang 0 -epmn 11 -epmx 11 -br 100 -tw 10 -sr -kb
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


## Acknowledgements

We acknowledge support from EPSRC (JCF), MRC (BB and AS), the ERC grant NEPA (JK and AS), the Royal Society (AS), the UK Materials and Molecular Modelling Hub for computational resources, which is partially funded by EPSRC (EP/P020194/1)

Peter Wirnsberger for the implementation of the Yuan membrane model in LAMMPs https://github.com/pw359

## Citations

1. Huajian Gao, Wendong Shi, and Lambert B Freund. Mechanics of receptor-mediated endocytosis
2. Tine Curk, Peter Wirnsberger, Jure Dobnikar, Daan Frenkel, and Andela Šaric. Controlling cargo trafficking in multicomponent membranes. Nano Letters, 18(9):5350–5356, 2018. doi:10.1021/acs.nanolett.8b00786. PMID: 29667410.
3. Hongyan Yuan, Changjin Huang, Ju Li, George Lykotrafitis, and Sulin Zhang. One-particlethick, solvent-free, coarse grained model for biological and biomimetic fluid membranes. Phys. Rev. E, 82:011905, Jul 2010. doi: 10.1103/PhysRevE.82.011905.
4. Steve Plimpton, Aidan Thompson, Stan Moore, and Axel Kohlmeyer. Lammps documentation, May 2017.
5. Michael Affenzeller, Stephan Winkler, Stefan Wagner, and Andreas Beham. Genetic Algorithms and Genetic Programming: Modern Concepts and Practical Applications. Chapman & Hall/CRC, 1st edition, 2009. ISBN 1584886293, 9781584886297.
6. Gerhard Goos, Juris Hartmanis, Jan van Leeuwen, David Hutchison, Takeo Kanade, Josef Kittler, Jon M Kleinberg, Friedemann Mattern, John C Mitchell, Moni Naor, and et al. Lecture notes in computer science. page 960.
7. Thomas M. J. Fruchterman and Edward M. Reingold. Graph drawing by force-directed placement. Software: Practice and Experience, 21(11):1129–1164, Nov 1991. ISSN 00380644, 1097024X. doi: 10.1002/spe.4380211102.
8. Gang Ren, Gabby Rudenko, Steven J Ludtke, Johann Deisenhofer, Wah Chiu, and Henry J Pownall. Model of human low-density lipoprotein and bound receptor based on cryoem. Proceedings of the National Academy of Sciences, 107(3):1059–1064, 2010.
9. Xiangxi Wang, Jingshan Ren, Qiang Gao, Zhongyu Hu, Yao Sun, Xuemei Li, David J Rowlands, Weidong Yin, Junzhi Wang, David I Stuart, et al. Hepatitis a virus and the origins of picornaviruses. Nature, 517(7532):85, 2015.
10. Jodi A Hadden, Juan R Perilla, Christopher John Schlicksup, Balasubramanian Venkatakrishnan, Adam Zlotnick, and Klaus Schulten. All-atom molecular dynamics of the hbv capsid reveals insights into biological function and cryo-em resolution limits. Elife, 7:e32478, 2018.
11. Bärbel Kaufmann, Alan A Simpson, and Michael G Rossmann. The structure of human parvovirus b19. Proceedings of the National Academy of Sciences, 101(32):11628–11633, 2004.
12. Laura Filion and Marjolein Dijkstra. Prediction of binary hard-sphere crystal structures. Physical Review E, 79(4):046714, 2009.
13. Gernot J Pauschenwein and Gerhard Kahl. Clusters, columns, and lamellae—minimum energy configurations in core softened potentials. Soft Matter, 4(7):1396–1399, 2008.
14. Jure Dobnikar, J Fornleitner, and G Kahl. Ground states of model core-softened colloids. Journal of Physics: Condensed Matter, 20(49):494220, 2008.
15. Kathrin Muller, Natan Osterman, Dusan Babic, Christos N Likos, Jure Dobnikar, and Arash Nikoubashman. Pattern formation and coarse-graining in two-dimensional colloids driven by multiaxial magnetic fields. Langmuir, 30(18):5088–5096, 2014.
16. Jian Qin, Gurdaman S Khaira, Yongrui Su, Grant P Garner, Marc Miskin, Heinrich M Jaeger, and Juan J de Pablo. Evolutionary pattern design for copolymer directed self-assembly. Soft Matter, 9(48):11467–11472, 2013.
17. Sam Kriegman, Douglas Blackiston, Michael Levin, and Josh Bongard. A scalable pipeline for designing reconfigurable organisms. Proceedings of the National Academy of Sciences, 2020.
18. Babji Srinivasan, Thi Vo, Yugang Zhang, Oleg Gang, Sanat Kumar, and Venkat Venkatasubramanian. Designing dna-grafted particles that self-assemble into desired crystalline structures using the genetic algorithm. Proceedings of the National Academy of Sciences, 110(46):18431–18435, 2013.
19. Marc Z Miskin and Heinrich M Jaeger. Adapting granular materials through artificial evolution. Nature materials, 12(4):326, 2013.
20. Robert Vácha, Francisco J. Martinez-Veracoechea, and Daan Frenkel. Receptor-mediated endocytosis of nanoparticles of various shapes. Nano Letters, 11(12):5391–5395, Dec 2011. ISSN 1530-6984, 1530-6992. doi: 10.1021/nl2030213.
21. Kai Xiong, Jiayin Zhao, Daowen Yang, Qingwen Cheng, Jiuling Wang, and Hongbing Ji. Cooperative wrapping of nanoparticles of various sizes and shapes by lipid membranes. Soft Matter, 13(26):4644–4652, 2017. ISSN 1744-683X, 1744-6848. doi: 10.1039/C7SM00345E.
22. Drew R. Elias, Andrei Poloukhtine, Vladimir Popik, and Andrew Tsourkas. Effect of ligand density, receptor density, and nanoparticle size on cell targeting. Nanomedicine: Nanotechnology, Biology and Medicine, 9(2):194–201, Feb 2013. ISSN 15499634. doi:10.1016/j.nano.2012.05.015.
23. Veronika Schubertová, Francisco J. Martinez-Veracoechea, and Robert Vácha. Influence of ligand distribution on uptake efficiency. Soft Matter, 11(14):2726–2730, 2015. ISSN 1744-683X, 1744-6848. doi: 10.1039/C4SM02815E.
24. N. Sloane. Tables of sphere packings and spherical codes. IEEE Transactions on Information Theory, 27(3):327–338, May 1981. ISSN 0018-9448. doi: 10.1109/TIT.1981.1056351.
25. Félix-Antoine Fortin, François-Michel De Rainville, Marc-André Gardner, Marc Parizeau, and Christian Gagné. DEAP: Evolutionary algorithms made easy. Journal of Machine Learning Research, 13:2171–2175, jul 2012.
