# GWIC-with-An-Efficient-Data-Structure
GWIC with An Efficient Data Structure


This is the improvement implementation of the Generalized-Weak-incNGs-Consistency (GWIC) on the incNGs global constraint proposed in the paper "Jimmy H.M. Lee and Zichen Zhu. Filtering Nogoods Lazily in Dynamic Symmetry Breaking During Search, Proceedings of the 24th International Joint Conference on Artificial Intelligence (IJCAI 2015), pages 339-345, Buenos Aires, Argentina, July, 2015". 

Please refer "Speeding Up Symmetry Breaking During Search and Partial Symmetry Breaking in Constraint Programming, Zichen Zhu, PhD Thesis, The Chinese University of Hong Kong, 2016" for detailed description.

Please use Gecode Solver 4.2.0 to run these files.

Put GOG folder into gecode folder and efpa_lresbds_GOG.cpp file into example folder.

To run the EFPA problem (e.g. (3 4 6 4)) using one of the symmetry breaking methods NONE/DoubleLex/SnakeLex/LDSB/LReSBDS:

./efpa_lresbds_GOG -search none 3 4 6 4

./efpa_lresbds_GOG -search double 3 4 6 4

./efpa_lresbds_GOG -search snake 3 4 6 4

./efpa_lresbds_GOG -search ldsb 3 4 6 4

./efpa_lresbds_GOG -search lresbds 3 4 6 4
