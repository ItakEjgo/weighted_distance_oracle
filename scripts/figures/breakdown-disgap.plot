#The template for generating multiple plots
#sharing the same legend in the same figure 

#The figure setting (applied to all plots)
#=========================================

set out "out/breakdown-disgap.eps"
set terminal postscript portrait enhanced mono "Helvetica" 40

set size 3.9000, 0.870000
# set size 2.50000, 0.800000
set pointsize 5.000000

set multiplot layout 2,4

set size 3.900000,0.100000;  
# set origin -0.15,0.62;
set origin -0.45, 0.77

set key inside top center horizontal samplen 2 width 1
# set key top center samplen 1.5
# set key vertical spacing 0.8

set yrange [0.0000000001:0.0000000002]

unset border
unset tics
unset label
unset xlabel
unset ylabel
set bmargin 1

plot NaN title 'FS(Fixed Scheme)' with linespoints linetype 1 pointtype 1, \
NaN title 'US(Unfixed Scheme)' with linespoints linetype 1 pointtype 4, \
NaN title 'K-Algo' with linespoints linetype 1 pointtype 8, \
NaN title 'SE-Oracle' with linespoints linetype 1 pointtype 10, \
NaN title 'EAR-Oracle' with linespoints linetype 1 pointtype 6#, \ 
# NaN title 'EAR-Oracle' with linespoints linetype 1 pointtype 6#, \
# # NaN title 'Theoretical bound ({/Symbol e})' with lines linetype 2
# #NaN title 'Hybrid' with linespoints linetype 1 pointtype 6, \
# #, \NaN title 'RN-Adapt' with linespoints linetype 1 pointtype 4, \
# #NaN title 'Cao-Appro2-Adapt2' with linespoints linetype 1 pointtype 12, \
# #NaN title 'Long-Appro-Adapt2' with linespoints linetype 1 pointtype 14

set border
set tics
set label

#Plot default query time
#=========================================
set size 1.90000,0.7000000;  
set origin 0.0,0.0;
# set origin 0.0, 0.0;

set xlabel  "dataset"
# set xlabel font ", 30"
# set xtics font ", 25"
set ylabel  "Time Cost (ms)"
set key top right vertical
set label

set yrange [0.1 : 100000]
set xrange [0 : 80]
set label 11 center at graph 0.5,char 1 "(a) Breakdown Time" 
set bmargin 5
set format x "%g"
set format y "10^{%T}"

# set xtics ("AL" 5, "" 10, "HE" 15, "" 20, "HO" 25, "" 30,  "BH" 35, "" 40,  "EP" 45, "" 50,  "SF" 55)
set xtics ("HM" 5, "" 10, "BM" 15, "" 20, "HL" 25, "" 30,  "RM" 35, "" 40,  "GF" 45, "" 50,  "LM" 55, "" 60, "BH" 65, "" 70, "EP" 75)
set log y

# plot ARG1 using ($1):($2 * $3):(1.5) title 'T_{G^P}' lt 1 with boxes fs solid 0.6, \
# ARG1 using ($1):($2):(1.5) title 'T_{G^C}' lt 1 with boxes fs solid 0.25
plot ARG1 using ($1-1.5):($2):(1.5) title 'graph construction' lt 1 with boxes fs solid 0.6, \
ARG1 using ($1):($3):(1.5) title 'graph processing' lt 1 with boxes fs solid 0.25

# ARG1 using ($2 + 1.5):11:(1.5) title 'EAR-Oracle' lt 1 with boxes fs pattern 3
# ARG1 using ($2):9:(1.5) title 'K-Algo' lt 1 with boxes fs pattern 1, \
# ARG1 using ($2 + 1.5):10:(1.5) title 'SE-Oracle' lt 1 with boxes fs pattern 5, \

#The 2nd plot (notitle)
#=========================================
# set size 1.25000,0.6000000;  
set size 1.90000,0.7000000;  
set origin 1.9, 0.0;

set xlabel  "Query Group"
set ylabel  "Running Time (ms)"

set xrange [0.020000: 0.53000000]
set yrange [1: 1500000]
set label 11 center at graph 0.5,char 1 "(b) Effect of Query Distance" 
set bmargin 5
set format x "%g"
set format y "10^{%T}"
# set xtics font ", 25"
set xtics ( "Q0" 0.05, "Q1" 0.1, "Q2" 0.15, "Q3" 0.2, "Q4" 0.25, "Q5" 0.30, "Q6" 0.35, "Q7" 0.40, "Q8" 0.45, "Q9" 0.50)
# set key font ", 20"
unset key
# set key inside top center samplen 2
# set key horizontal spacing 1
set log y

plot ARG2 using 1:3 title "FS" with linespoints linetype 1 pointtype 1, \
ARG2 using 1:4 title "US" with linespoints linetype 1 pointtype 4, \
ARG2 using 1:5 title "K-Algo" with linespoints linetype 1 pointtype 8, \
ARG2 using 1:6 title "SE" with linespoints linetype 1 pointtype 10, \
ARG2 using 1:2 title "EAR" with linespoints linetype 1 pointtype 6

unset multiplot