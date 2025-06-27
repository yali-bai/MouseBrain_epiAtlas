
### number_statstics_for_DMRs.xlsx  :  
Statistical results of DMR number.

Three sheets are included, corresponding to the three levels on which DMRs are called.
The sheet '01_threeclass' records the number of DMRs called among excitory neurons, inhibitory neurons, and non-neurons.
The sheet '02_subclasses_in_neuron' records the number of DMRs called among the neuronal subclasses.
The sheet '03_subclasses_in_NN' records the number of DMRs called among the subclasses of non-neurons.

In every sheet,
the rows represent the (sub)classes, with a final row 'total' recording the numbers for DMRs in any subclass;
the columns, '5mCG hyper', '5mCG hypo', '5mCG total', '5hmCG hyper', '5hmCG hypo', '5hmCG total',
represent the DMR types,
with '5mCG' signifying 5hmCG+5mCG, '5hmCG' signifying 5hmCG,
'hyper' and 'hypo' indicating the state of the DMRs, and 'total' indicating hyper or hypo.

Note that one segment can be in the same or different states in different subclasses,
(for example, one segment is a hyper-DMR in 2 subclasses, or is a hyper-DMR in one subclass and a hypo-DMR in another subclass,)
which causes differences between the sum of DMR numbers in subclasses and the DMR number in the 'total' row.