To automate feature transformation with a merged list containing multiple ili years of joints, we built an interconnected program of three layers. As long as the data we are testing is also a merged list, we only need to modify the global parameters in Layer 1. This will maximize processing speed by reducing repitions and avoid creating errors. 

# Input Data
{Joint Number,DS_TO_US(m),Orientation(deg/clockwise position)}

# Parameters
iliyr(fix & move),  diameter

# Output Data
Input Data + {Orientation_after_trans(m),DS_TO_US_after_trans(m)}

# Structure

- Layer 1: Iteration over multiple yeras and global parameter setup
- Layer 2: Parameter setup for each ili yera
- Layer 3: Transformation, main script
