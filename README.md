# Hypergraph 1-Spectral Clustering with General Submodular Weights

This repo contains code for *Zhu et al., Hypergraph 1-Spectral Clustering with General Submodular Weights, Asilomar 2022*.

## 0. Prerequisites
The implementation is based on the following software and libraries:

1. MATLAB R2022b - https://www.mathworks.com/products/matlab.html

2. Python 3.9 & conda - https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-python.html

3. Eigen 3.40 - https://eigen.tuxfamily.org/index.php?title=Main_Page

## I. Datasets and prepare the data
### i) RCV1
1. Download from https://github.com/XifengGuo/DEC-keras/blob/2438070110b17b4fb9bc408c11d776fc1bd1bd56/data/reuters/get_data.sh, save to `./python/data/RCV1/`.

2. Set sparsity parameters and number of keywords (hyperedges) in `main_rcv1.py`, or use the default ones.

3. Run `python main_rcv1.py`, and the processed dataset will be saved as `./matlab/data/rcv1/data.mat`.

### i) Covtype
1. Dataset and description can be found via https://archive.ics.uci.edu/ml/datasets/covertype.

2. For your convenience, this dataset is available at `./matlab/data/covtype/covtype.mat`.
  

## II. Test the proposed and baseline methods

1. Launch MATLAB in the `./matlab/` folder.

2. In `setup.m`, change the path to your installed Eigen library.

3. Run `setup` in MATLAB command window. This should give you mex-ed functions. 

4. Try running `test_rcv1(2, 0.2, 0)`, or see all used p (corresponding to the parameter alpha in our paper) and delta (corresponding to the parameter beta in our paper) choices in `main.m`.

5. Checkout saved results in `./matlab/data/DATASET_NAME/NAMED_AFTER_YOUR_CONFIGS/...`.

## III. Contact

If you have any questions, please contact yz126@rice.edu.
