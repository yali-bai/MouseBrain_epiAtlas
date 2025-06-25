### 00_split_h5ad_by_chromosome
The ipynb file in "00_split_h5ad_by_chromosome" is first used to split the h5ad generated from mcds by chromosome.
This allows parallelization, and reduces the requirement for memory during DMR calling.

### 01_threeclass_1vsOthers, 02_subclasses_in_neuron_1vsOthers, 03_subclasses_in_NN_1vsOthers
Then the files in "01_threeclass_1vsOthers", "02_subclasses_in_neuron_1vsOthers", "03_subclasses_in_NN_1vsOthers" are used to conduct two-sided Wilcoxon rank-sum test for DMR-calling in 3 stratified levels: 
(1) among excitatory neurons, inhibitory neurons, and non-neurons; 
(2) among all neuronal subclasses; 
(3) among all non-neuron subclasses.

### run_RankSumTest_1vsOthers.py
The 'run_RankSumTest_1vsOthers.py' script submits jobs in slurm by calling the 'to_RankSumTest_1vsOthers.sh' script, which calls the 'work_RankSumTest_1vsOthers.py' to actually conduct DMR calling in Python.