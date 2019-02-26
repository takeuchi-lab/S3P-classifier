# Description

This code was used in the experiments in the following papers. For details of the algorithm, please see the following papers.

Efficient learning algorithm for sparse subsequence pattern-based classification and applications to comparative animal trajectory data analysis  
[https://www.tandfonline.com/doi/full/10.1080/01691864.2019.1571438](https://www.tandfonline.com/doi/full/10.1080/01691864.2019.1571438)

Note that this code is created based on Nakagawa's code ([https://github.com/takeuchi-lab/SafePatternPruning](https://github.com/takeuchi-lab/SafePatternPruning)).  


# Verified Environmental

* gcc version 5  
* GNU Make 4.1  

# Setup

* make

# Usage
`./train [-options] input_file`

options:

-   -u : problem type (default 1)
    -   1 -- regularization path computation for L1-reguralized L2-SVM  
    -   2 -- regularization path computation for Lasso  
-   -t : learning lambda index(when do not cv) (default:most minimum lambda)  
-   -m : minimum supportSum (default 1)  
-   -L : maximum length of pattern (default 10)  
-   -T : the number of regularization parameter (default 100)  
-   -r : lambda min ratio (default 2.0)  
-   -i : max outer iteration in optimization (default 1000000)  
-   -f : frequency of calculate duality gap and convergence check (default 50)  
-   -e : convergence criterion of duality gap (default 1e-6)  
-   -F : name of reslut file (default output/result.csv)  
-   -p : maximum interval of event (default 0|-1:none)  
-   -c : whether to do cross validation (default 0:do not|1:do)  
-   -k : k-fold of k (when do cross validation)(default:10)  
-   -a : times to do cross validation(default:10)  
-   -C : whether to do CloSpan (default 0:do not|1:do)  
-   -M : whether to do Multiprocess  (default 0:do not|1:do)  
-   -P : whether to do cross validation for Performance evaluation (default 0:do not|1:do)  
-   -s : the Mode of counting supportSum(default 0),0 is 0 or 1 per record,1 is the number of pattern  

## Example
`./train -L 100 -F ./output/result.csv ./data/sequence.txt`

If you decide the best lambda by cross-validation and want to learn the model, give the option as follows.  
Note: by default, 10-fold cross-validation is executed 10 times.   
(Since cross-validation takes time, we recommend enabling M option and executing in multi-process.)

`./train -L 100 -c 1 -M 1 -F ./output/result.csv ./data/sequence.txt`

# Demo

##Input data

As shown in the following example, input data should have a label at the beginning and a sequence of symbols behind it.
Note that separate each column with space, and the symbol must a positive integer.  

>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 2 2  
>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 3  
>1 1 1 1 1 2 2  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 2 2  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 3  
>1 1 1 1 1 1 2 2  
>-1 1 1 1 1 1 2 2 3  
>-1 1 1 1 1 1 2 2 3  
>-1 2 2 3  
>-1 2 2 3  
>-1 2 2 2  
>-1 1 1 1 1 1 2 2 3  
>-1 1 1 1 1 1 2 2 3  
>-1 2 2 3  
>-1 2 2 3  
>-1 2 2 2  
>-1 1 1 1 1 1 2 2 3  
>-1 1 1 1 1 1 2 2 3  
>-1 2 2 3  
>-1 2 2 3  
>-1 2 2 2  
>-1 1 1 1 1 1 2 2 3  
>-1 1 1 1 1 1 2 2 3  
>-1 2 2 3  
>-1 2 2 3  
>-1 2 2 2  

For example, the line  
>1 1 1 1 1 3  

indicates that the label is "1" and the sequence "1 1 1 1 3".

## Output

For the output result, the total support, the lambda when that the pattern was added, the weight, the ID and support of the data including that pattern are output in CSV for each extracted pattern as follows.

>mPattern,supportSum,add lambda,w,list[id:sup]  
>1 1 1 1,28,15,0.741554,[0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1,14:1,15:1,16:1,17:1,18:1,19:1,20:1,21:1,25:1,26:1,30:1,31:1,35:1,36:1,]  
>1 1 1 1 1 2 2 3,8,35,-0.538109,[20:1,21:1,25:1,26:1,30:1,31:1,35:1,36:1,]  
>1 1 1 1 3,16,1,0.03,[0:1,1:1,2:1,3:1,5:1,6:1,7:1,8:1,10:1,11:1,12:1,13:1,15:1,16:1,17:1,18:1,]  
>2 2 2,4,35,-1.17845,[24:1,29:1,34:1,39:1,]  
>2 2 3,16,1,-1.40189,[20:1,21:1,22:1,23:1,25:1,26:1,27:1,28:1,30:1,31:1,32:1,33:1,35:1,36:1,37:1,38:1,]  


-----------------------------------------------------------------------
# Licence
GNU General Public License v3.0
