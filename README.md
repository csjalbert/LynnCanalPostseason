# LynnCanalPostseason

**ARCHIVE ONLY - this code was used for prior to creation of an 'all-encompassing' MAGMA script. It is here for archive purposes but is no longer used or maintained. CJ 6/16/2020**

This repository holds all scripts associated with analyzing sockeye salmon mixtures from Lynn Canal, D115 commercial samples as part of our postseason genetics project from 2018 forward. Specifically, it houses our MAGMA (Mark, Age, Genotype, Mixture Analysis) pipeline.

<b> Project Overview </b> 

This is a post-season analysis of Lynn Canal D115 fishery. We want an age-specific stock composition for all major contributing age classes (>5%; 1.2, 1.3, 2.2, 2.3, other) using mark- and age-enhanced genetic mixed-stock analysis (MAGMA).    

For the purposes of this project, we are concerend about the breakdown by age, across all statweeks, and statareas. I will assign <b>all</b> stat areas as 115-<i>00</i>, since this subdistrict doesn't exist and is currently how stat area is recorded.    
 
Total season estimates are to be provided by age group. The MAGMA algorithm will be run for 40 000 repititions, and the first 20 000 repititions will be discarded.     

Final deliverables: 
    
    Point estimates, Credibility intervals

There are 4 steps to this pipeline which involve data preparation, cleaning, model setup, model running, and data summary. 
Very briefly, the steps are as follows: 
# 1- Data Cleanup 
This takes creates the directory structure, analyizes ASL, Genetic, and harvest data, and creates population grouping objects. 
# 2- Workspace Setup
This creates the workspaces required for running the MAGMA model. 
# 3- Model Run
Runs the MAGMA model...
# 4- Summary
Converts model output into something more usable (R or Excel format) and performs QC steps to check for model issues. Finally, summary tables are created for distribution / reporting. 
