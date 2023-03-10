# PredLLPS_PSSM
## About PredLLPS_PSSM

we proposed a novel predictor (PredLLPS_PSSM) for LLPS protein identifying only based on sequence evolution information.

The datasets can be found in `./data/`. The PredLLPS_PSSM models is available in `./model/`. The prediction code can be found in `mian.py`.

Your sequences must be submitted in fasta format and the csv file of the AB PSSM features you got from POSSUM. Please ensure that the fasta file submitted to POSSUM is the same as the fasta file submitted this time. POSSUM's web site is https://possum.erc.monash.edu.
1. The submitted sequence length should be no less than 50 and no longer than 5000.
2. The number of sequences submitted is within 500.
3. Submit A sequence of one and only 20 kinds of amino acids, including 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'.
