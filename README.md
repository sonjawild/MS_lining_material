# TITS (PARIDAE SP.) USE SOCIAL INFORMATION WHEN LOCATING AND CHOOSING NEST LINING MATERIAL  

## lining material.R 
contains the R code needed to replicate all analyses and plot figures

## coordinates_boxes 
contain GPS locations of all nest boxes and the five wool dispensers for calculating distance matrices

## dispenser.data.RDA 
contains raw data of all visits to dispensers 1-5
 columns:
- date.time (yymmddHHMMSS)
- Antenna (A for auxiliary, M for main)
- PIT: PIT tag of visiting female (10-digit alphanumeric code)
- Location: D1-D5 for each of the 5 dispensers
- visit: cumulative number of visits by birds

## wool choice.txt 
contains information on the wool females incorporated into their nest as a first colour and the majority colour they used
- Box: nest box name
- first_color: the first colour they incorporated
- second_color: the second colour they incorporated (empty if only one colour)
- PIT: PIT tag of breeding female (10-digit alphanumeric code)
- Species: GRETI (great tit), BLUTI (blue tit), MARTI (marsh tit), unknown
- Age: adult or first year
- Closest_disp: which dispenser was closest (Euclidean distance)
- Initial_col_provided: initally seeded colour in the respective dispenser area
- Matched: 1 if first colour matched seeded colour, 0 if no match
- Demos: 'yes' if built nest before access to both colour was granted, 'no' if otherwise

## ILVs.combined.RDA
contains individual-level variables for NBDA analyses
 - Box: nestbox number
 - PIT_f: PIT tag of female
 - Species: GRETI for great tit, MARTI for marsh tit, BLUTI for blue tit
 - Age: either 'adult' or 'first.year'
 - Lay.Date: Date of first egg laid in number of days after the first of April
 - Hatch.Date: date of first chicks hatching in number of days after first of April
 - First.visits: Date (yymmddHHMMSS) of first visit to a dispenser
 - D1-D5: distance (in m) to each dispenser
 - D1.visited-D5.visited: 0 if not visited dispenser, 1 if registered on the respective dispenser
 - closest.dispenser: distance (in m) to the closest dispenser
 
 ## gmm.spring
 contains raw data for calculating foraging associations 
 - slot $gbi = group by individual matrix (each row is one individual, each column is a group. Entry of 1 if present in group, 0 if absent)
 
