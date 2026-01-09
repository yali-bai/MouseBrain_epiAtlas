
The preprocess of TSO-joint datasets

![image](https://github.com/yali-bai/MouseBrain_epiAtlas/blob/main/00.TSO-joint.pipeline/data/TSO-joint.png)
### 1. Separate single cell 
Separate single cell from a mixture plate by single cell barcode

There are two cell barcode libraries, including 96 barcode and 384 barcode. The corresponding barcode was selected according to the experimental scheme and split
 
```
run.tso.384barcode.sh 
run.tso.96barcode.sh
```

### 2. Preprcoss single cell
Processing of the splited single cell raw read data, including trim, Align, call methylation et. al.
 

```
run.forDNA.PE_SE.sh
run.forRNA.PE.sh
```



