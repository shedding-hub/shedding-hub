Data and R syntax for the following publication:
Tim K. Tsang, Vicky J. Fang, Kwok-Hung Chan, Dennis K. M. Ip, Gabriel M. Leung, J. S. Malik Peiris, Benjamin J. Cowling, Simon Cauchemez. Individual correlates of infectivity of influenza A virus infections in households. PLoS ONE 2016.


+-----------------+
|  List of files  |
+-----------------+

A. [data.csv] This data file contains information on the information recrutied from the household transmission study in Hong Kong, from 2008-2012.
Definition of column headings >>>
1. Household ID
2. Member ID (0 is index case)
3. Number of household member in the household
4. PCR-confirmed infection
5. Date of symptom onset ( 1/1/2008 is set as Day 1)
6. Blank column for adding variable in analysis
7. Delay between home visit and symptom onset in index cases
8. Age (years)
9. Sex (1: male)
10. Vaccination
11. Date of symptom onset of the index case
12. The end of the follow-up date
13. Intervention for the household (1: control, 3: hand hygiene, 4: hand hygiene + facemask)
14. Antiviral treatment
15. Asymptomatic infection
16-28. Observed viral shedding from Day 0 to Day 12. (Day 0 refers to the symptom onset date in index cases)
29. Chronic disease 
30. Smoking
31. PCR-confimred infection of influenza A virus
32. Influenza virus subtype (0:Unsubtyple influenz A 1: Seasonal A(H1N1) 2: Seasonal A(H3N2) 3:Pandemic A(H1N1) ) 
33. Households with co-index 
34. Presence of fever
35. Presence of sore throat
36. Presence of cough
37. Presence of runny nose
38. Presence of phlegm
39. Presence of muscle pain
40. Presence of headache
41. Blank column for adding variable in analysis
*** in all data file, -1 represents missing

B. [ILILAB.csv] The influenza proxy used in the analysis. ( 1/1/2008 is set as Day 1)

C. [lmec regression.r] R syntax to use log-linear random effects censored regression model to fix the observed viral shedding data and divide index cases into three viral shedding groups. This will create a data file [grouped_data.csv] for further analysis.
Definition of column headings >>>
All of the same from [data.csv] except the following:
16-28. Predicted viral shedding trajectories from the log-linear censored regression model
33. The viral shedding groups of index cases (1. The lower viral shedding groups 2. The medium viral shedding groups 3. The higher viral shedding groups)
41. Households with co-index

D. [Model_1.m] Contain the matlab function to generate results on the model in the main analysis

E. [Model_2.m] Contain the matlab function to generate results on the additional analysis of exploring the effect of presence of symptom. (a to i respents fever, sore throat, cough, runny nose, phlegm, muscle pain, headache, any three or more of above symptoms and ILI (fever + cough or sough) )

F. [Figure_1.r] R syntax to reproduce Figure 1.

G. [Table_1.r] R syntax to reproduce results in Table 1.

H. [Table_2.r] R syntax to reproduce results in Table 2.

I. [Table_3.r] R syntax to reproduce results in Table 3.



