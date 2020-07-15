# Super-Effects

Lead: Kelly Garner
Co-authors: Abbey Nydam, Zoie Nott 
Working Title: The contribution of sample size to the reproducability of individual differences methods in cognitive psychology 

Project Overview:

Th aim is to assess the generalizability of effect sizes across different sample sizes. 
Inspired by Howard Bowman's paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6018568/
Also see this paper by Tal Yarkoni: https://psyarxiv.com/jqw35effect 

Backround: 

Data taken from the Executive Function and Implicit Learning (EFIL) project run by Abbey Nydam and Paul Dux for the Team Honours Thesis in 2019.
It included data from ~350 1st-year psychology participants on a battery of seven cognitive tasks: 

Tasks were:
1. Operation Span
2. Attentional Blink
3. Go No-Go Task
4. Single Dual Task
5. Contextual Cuing
6. Serial Reaction Time Task
7. Visual Statistical Learning

These were organised into 4 tasks that tap executive functions (OS, AB, GNG, SD) and 3 that tap implicit learning (CC, SRT, VSL). The aim of the EFIL study was to look for relationships between and across EF and IL task, using an individual differnces approach (i.e., correlating individual scores across tasks). We found correlations between OS and AB, AB and VSL (I think - Abbey to check).

Planned Analysis:

Effect sizes will be sampled across Ns, and we will compare these distributions when the task effects were analysed using anova/t-tests versus linear mixed effects modelling. 

Correlation effect sizes between tasks will be sampled across Ns.

Data Structures:

Original data is in .mat files (hosted by Abbey - 1 file per 7 tasks)
Trial level data is in .xls or .csv format (analysed by MATLAB scripts from Abbey)
Condition level data is in .xls or .csv format  (analysed by MATLAB scripts from Abbey)

Analysis Scripts:

ANOVA & t-test scripts are in MATLAB scripts from Abbey
mixed-effects analysis - yet to be generated 
effect-size resampling scripts - example in R written by Kelly (for SRT task)
