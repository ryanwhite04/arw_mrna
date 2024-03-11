# arw_mrna
Adaptive Random Walk for the mRNA design problem

Depends on the ViennaRNA package: `pip install ViennaRNA`
Tested using ViennaRNA==2.6.4

Run `run_awalk.py` from the `src/` directory to design an mRNA sequence
Example usage: `python3 run_awalk.py --aa_seq=MVSKGEELFTGVVPILVELDGDVNGH --verbose --stability=efe --cai_threshold=0.9`