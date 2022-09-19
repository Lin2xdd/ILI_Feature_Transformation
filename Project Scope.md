## Project Description
The **`ILI_Transformation`** project is to develop a tool to assist multi-year ILI feature matching by applying affine transformation on feature data


## Objectives
1. given a merged list as input, the program will perform affine transformation on multiple years of ILI feature data in matched joints iteratively. output an updated merge list with modified `distance to us girth weld` and `orientation`  
 

## Work Decomposition 
1. Test on real datasets
2. Automate the parameter tunning without manual intervention
3. Standardize the input and output of the transformation
   1. [input format](../../data/input/featurematch_input.csv)
   2. output format will be the same as the input but give an updated `DS_TO_US` and `Orientation`
4. Evaluate the transformation and identify anomalies in transformation
5. Add iteration process to perform the transformation on multiple joints
6. Add iteration process to perform the transformation on multiple years 


## Deliverables
1. code file which can be integrated to the main ILI project
2. documents on how each function works 
3. transformation evaluation and anomaly detection proccedure
4. Reports on test dataset

## Milestones


## Resources
