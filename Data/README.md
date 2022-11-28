This file contains all five(5) datasets we used to test the performance of the algorithm. 
Fix year is usually the later iliyear and is used as reference; Move year is the opposite.

# Simulated Data:
- Move: EE_ILI1_031
- Fix: EE_ILI2_037
  
# Real Data:
  Multi-year: Test_multiple year ILI
  Two years: Test13_1_2022_2020 (fix yera 2022, move year 2020); Test8_1_2020_2019 (fix year 2020, move year 2019)
  
# Pre-processing
Before throwing the data in the algorithm, we need to transform
- Orientation

  Clockwise - Degree - Arch
  
  Convert clock position to Degrees and then calculate the circumferential distance in meters. The diameter of the pipeline is 168.28mm. below is a reference on how to convert hrs:mins to degrees/
"Convert clock position in hours / minutes or time format to circumferential degrees. Copy the entire clock position column and paste into Notepad. Insert a column to the right to receive the conversion calculation. Set the format of the original clock position column to "Text" and the new column for orientation degrees to "Number - 0 decimals". Select All and Copy from Notepad and Paste back to the original clock position column. Filter out the blanks from the clock position column to only rows with values. Use the formula:

=IF(ROUND((LEFT(A1,FIND(":",A1)-1)*30)+(RIGHT(A1, LEN(A1) - FIND( ":",A1 ) )/2),0)<=360,ROUND((LEFT(A1,FIND(":",A1)-1)*30)+(RIGHT(A1, LEN(A1) - FIND( ":",A1 ) )/2),0),ROUND((LEFT(A1,FIND(":",A1)-1)*30)+(RIGHT(A1, LEN(A1) - FIND( ":",A1 ) )/2),0)-360)

*Edit the “A1” to the appropriate clock position source row and column in the formula. Copy / paste the formula to bottom. Unfilter the column. Copy the full column and paste “as values” back to the same column. Find all 360 values and Replace with 0.

"
- Distance to US

- In first try, we also used $\Delta$Distance
  Reset the start and calculate the difference in distance (Dist_i - Dist_0). Dist_0 should be the first row in the raw dataset
