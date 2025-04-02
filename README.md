PhyMapMethyl_v0.1.1  
(Copyright 2025, Authors and University of Mississippi; see license below)

Updated April 1, 2025
==================

The PhyMapMehyl maps methylation shift for each branch of a tumor tree [1]. PhyMapMethyl is written in Python and works smoothly on Linux at present.  It is under continuous development, so anticipate releasing a version for Linux in the near future.  You are free to download, modify, and expand this program under a permissive license similar to the BSD 2-Clause License (see below). 


Dependencies
==================
1. R (version 4.3.3 was tested)

    l1ou (https://github.com/khabbazian/l1ou)
   
    phytools

3. python (version 3.11.5 was tested)


How to prepare input files
==================
1. beta value.
   
    A single file should be created for each sample. All samples should be ID-Sample.txt. Normal samples should be ID-N.txt and ID-N1.txt. The following format should be used:

    #ID_REF = 		
  
    #VALUE = Average Beta		
  
    ID_REF	VALUE	Detection Pval
  
    cg19252052	0.450634308381794	0
  
    cg19252053	0.126666703029043	0
  
    cg19252054	0.86149334281374	0
  
    cg19252055	0.950406048332555	0
  
    ...
  
    Example file: Example\input\Tree6-T1.txt. ID is Tree6 and Sample is T1 in this case.

2. Tumor tree (nwk).
   
    Tree should be time tree without normal samples.

    Example file: Example\Tree.nwk


How to run PhyMapMethyl
==================
Open command prompt and run the command, 

python3 PhyMapMethyl.py [path to input beta value files] [path to tree file] "ID" 

For example: 

To perform the example data analysis, try:

python3 PhyMapMethyl.py Example/input Example/Tree.nwk Tree6


Output file
==================
All the output files can be found found in the same folder as the input files. 

1. PhyMapMethyl inference
   
In PhyMap_Supported.txt, methylation shits that are detected are listed for each branch.

2. PhyMapMethyl inference (figure)
   
Heatmap of beta values at detected positions are produced for each branch.

3. Tree node ID table
   
Node map tables are saved in the same directory as the tree input file. 

 

--------
Reference:
[1] Deyana Tabatabaei, Sayaka Miura and Sudhir Kumar, The Pattern of Gain and Loss of DNA Methylation DuringTumor Evolution (2025) Under Review

--------
Copyright 2025, Authors and University of Mississippi
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
