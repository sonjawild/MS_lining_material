# MS_lining_material

# testing-things
## Sonja Wild

### lining material.R contains the R code needed to replicate all analyses and plot figures

### coordinates_boxes contain GPS locations of all nest boxes and the five wool dispensers for calculating distance matrices

### dispenser.data.RDA contains raw data of all visits to dispensers 1-5
# columns:
- date.time (yymmddHHMMSS)
- Antenna (A for auxiliary, M for main)
- PIT: PIT tag of visiting female (10-digit alphanumeric code)
- Location: D1-D5 for each of the 5 dispensers
- visit: cumulative number of visits by birds

### wool choice.txt contains information on the wool females incorporated into their nest as a first colour and the majority colour they used
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
