# Readme
This is the accompanying code for the paper *Continuous-time control synthesis under nested signal temporal logic specifications*  by Pian Yu, Xiao Tan, and Dimos V. Dimarogonas. In this code, we implement the the sTLT tree structure in the paper by automating the parsing of given STL formulas, the construction of an sTLT, and the synthesis of corresponding CBFs from a given sTLT. In the current code, only dynamical systems of type `singleIntegrator` or `unicycle`are implemented, but it should be easily adaptable to other dynamics by following the example code. There are several MATLAB objects defined in this repo: `setNodeObj`,`operatorNodeObj`,`stltObj`,`cbfObj`.

# How to run
You need to install both [CORA](https://tumcps.github.io/CORA/) and [helperOC](https://github.com/HJReachability/helperOC) to run this code. CORA is not necessary for the core part of the code, but is only used in drawing some circles. After installing them, add `\utility_func` to the path.

To reproduce Fig. 3 and Fig. 5 in the paper, run `main_singleIntegrator_Fig35.m`. The result will only show one trajectory. To see the other one, you need to manually change the initial state manually.

To reproduce Fig. 4 and Fig.6, run `main_unicycle_Fig4.m` and `main_singleIntegrator_Fig6.m`, respectively. It can take a while when reproducing unicycle examples.

In this code we used the reachability toolbox [helperOC](https://github.com/HJReachability/helperOC) to do reachability analysis and control barrier function synthesis, which unfortunately is not well documented. A brief introduction to this toolbox is given in `\docs`

# If you want to go deeper
The constructed CBFs are MATLAB objects. For single integrator dynamics, only circular regions are considered as regions of interests (ROI). `obj = singleIntegratorCBF(timeInterval,c,r,vMax,obs)`. The `obs` flag is needed only if the circular region is an obstacle. `timeInterval ` is  $$[\underline{t}_{\mathfrak{b}_i}, \underline{t}(\mathbb{X}_{f_i}), \bar{t}_{\mathfrak{{b}_i}}]$$ for a temporal fragment. `c` and `r` are the center and radius of the circular region, and `vMax` is the velocity bound. For unicycle models, there are two different ways to construct a CBF, given in `main_unicycle_Fig4.m` and `main_singleIntegrator_Fig6.m`. It is recommended to follow `main_singleIntegrator_Fig6.m`. The procedure is to first construct the set (represented by a grid and a super-level set from the data). This step is done by conducting the reachability analysis. Then `uniCBF = unicycleCBF(timeInterval,grid,data0,vRange,wMax)`. More details on the `uniCBF` and reachability analysis are given in two tutorials `tutorial_unicyleObj.m` and `tutorial_uni_ReachableOperation_RM_Rm.m`.



