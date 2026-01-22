- **shuffle.py**: Run python shuffle.py human and python shuffle.py mouse to perform 1000 label permutations on the benchmark datasets. The pre-generated permutation labels used in this study are stored in `human_shuffled` and `mouse_shuffled`, respectively.

- **svm.py**: Run this script to conduct permutation tests using the SVM model. Example: python svm.py human heart

- **mlp.py**: Run this script to conduct permutation tests using the MLP model. Example: python mlp.py human heart

- **permutation_test.py**: Run this script to compute and summarize the results of the permutation tests, which are used to assess whether model performance is significantly better than random expectation.