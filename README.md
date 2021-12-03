# lotek_processing
Functions and workflows for processing Lotek GPS collar data

### gps_cleaning.ipynb
A quick overview of the current output from gps_cleaning.ipynb:
1. For the two files described below, I made the following calculations FIRST, in order:
* Determined which pasture the fix is in from the shapefile (saved as ‘past_mask’)
* Calculated the distance (m) of each fix from the assigned pasture (saved as ‘distance_nearest’)
* Calculated the coordinates of the point within the assigned pasture that is nearest to the fix (saved as ‘UTM_X_nearest’, ‘UTM_Y_nearest’)
* Removed all fixes that where ‘distance_nearest’ > 100 m
* If the pasture of the fix location (‘past_mask’) was NOT the same as the assigned pasture, I changed the fix location to match ‘UTM_X/Y_nearest’. Otherwise I kept the coordinates as-is. These final coordinates are saved as ‘UTM_X/Y_final’
* Calculated step length (‘steplength’) in meters and turning angle (‘turnangle’) in degrees from straight line (0-180) using the ‘UTM_X/Y_final’ coordinates
* Calculated movement rate (‘moverate’) in meters per minute based on step length and step duration
2. The script then outputs the following two files in the output directory specified (outDIR)
  * *_flagged.csv: this has the following five flags (where 1 is ‘bad’ and 0 is ‘good’) included for each fix
    * jump_flag: suspected jump based on movement rate > 42 m/min and turn angle > 120 degrees
    * fast_flag: suspected error based on movement rate > 84 m/min
    * missingfix_flag: flag for entire day’s fixes based on > 20 missing 5-min fixes for the day
    * badfix_flag: flag indicating if any of the previous three flags exist. This is used for calculating grazing hrs per day (if any of the previous three flags are present, the fix is not included when calculating grazing hours per day)
    * grazinghrs_flag: flag for the entire day’s fixes based on grazing hours < 6 or > 13.
  * *_cleaned.csv: this file has removed all fixes with any of the flags above. This was done in order to be able to calculate grazing bouts and associated statistics. After removing flagged fixes, the following were calculated:
    * grazing_bout: this simply numbers each unique bout (both grazing and non-grazing), with a change in bout defined as: 
      * Change bouts when grazing activity changes, unless the two fixes before and two fixes after are the same. In detail, change to new bout if (1) activity is not the same as the previous row, (2) AND activity is the same as one of the next two rows (3) AND activity is not the same as one of the two rows after it. NOTE: this only works for bouts of 4+ fixes. See below for defining as ‘Transition’ bouts (see below).
    * bout_mins: calculates the duration of a bout in minutes from the first/last time stamps
    * bout_maj: the ‘GrazingAct’ value associated with the majority of fixes within a bout
    * bout_act: converting ‘GrazingAct’ values to ‘Nongrazing’ (0) and ‘Grazing’ (1). Bouts with a duration < 20 mins are manually reassigned to be labelled ‘Transition’. **I am not yet sure what happens in the rare case that a bout is less than 4 fixes long, but more than 20 mins long (i.e., due to missing fixes)! I need to investigate this.
    * act_bout_count_daily: calculates the number of bouts per day each activity. NOTE: this value gets repeated across all bouts of the same type within a day
    * act_budget_daily: calculates the proportion of total bout durations in each activity. NOTE: so this is not based on 24 hrs, but based on available fix data, which will be at least 23 hrs because of missingfix_flag above.
