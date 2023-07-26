#The template for generating multiple plots
#sharing the same legend in the same figure 

#The figure setting (applied to all plots)
#=========================================

#Origin file format
# col1  col2  col3  col4       col5      col6
#      fixedS KAlgo EAR-Oracle SE-Oracle unfixedS
# eps
# eps
# eps
# eps
# eps

set out "out/spnum-2row.eps"
set terminal postscript portrait enhanced mono "Helvetica" 32 

set size 2.5800, 1.240000
set pointsize 3.000000

set multiplot layout 2,4
#The first plot (which generates the common legend)
#===========================================================
#Note the NaN variable.

set size 2.50000,0.100000;  
set origin -0.150,1.18;

set key center top horizontal samplen 2 font "Helvetica, 28"
set yrange [0.0000000001:0.0000000002]

unset border
unset tics
#unset label

plot NaN title 'FS(Fixed Scheme)' with linespoints linetype 1 pointtype 1, \
NaN title 'US(Unfixed Scheme)' with linespoints linetype 1 pointtype 4, \
NaN title 'K-Algo' with linespoints linetype 1 pointtype 8, \
NaN title 'SE-Oracle' with linespoints linetype 1 pointtype 10, \
NaN title 'EAR-Oracle' with linespoints linetype 1 pointtype 6#, \
# NaN title 'Theoretical bound ({/Symbol e})' with lines linetype 2
#NaN title 'Hybrid' with linespoints linetype 1 pointtype 6, \
#, \NaN title 'RN-Adapt' with linespoints linetype 1 pointtype 4, \
#NaN title 'Cao-Appro2-Adapt2' with linespoints linetype 1 pointtype 12, \
#NaN title 'Long-Appro-Adapt2' with linespoints linetype 1 pointtype 14
 
set border
set tics
set label

#The 2nd plot (notitle)
#=========================================
set size 1.400000,0.6000000;  
set origin 0.0,0.55;

set xlabel  "Steiner points m"
set ylabel  "Building Time (s)"
set key left

set xrange [0.020000: 0.28000000]
set yrange [10: 10000]
set label 11 center at graph 0.5,-0.5, char 1 "(a) Building Time" 
set bmargin 5
set format x "%g"
set format y "10^{%T}"
# set xtics font "Helvetica, 32"

# set xtics ( "0.05" 0.05, "0.1" 0.1, "0.15" 0.15, "0.2" 0.2, "0.25" 0.25)
set xtics ( "3" 0.05, "4" 0.1, "5" 0.15, "6" 0.2, "7" 0.25)

set log y

plot ARG1 using 1:2  notitle with linespoints linetype 1 pointtype 10, \
ARG1 using 1:4 notitle with linespoints linetype 1 pointtype 6

#The 4th plot (notitle)
#=========================================
set size 1.400000,0.600000;  
set origin 0.0,0.0;

set xlabel  "Steiner points m"
set ylabel  "Size (MB)"
set key above

set xrange [0.020000: 0.28000000]
set yrange [100: 1000000]
set label 11 center at graph 0.5,-0.5, char 1 "(b) Space Consumption" 
set bmargin 5
set format x "%g"
#set format y "10^{%T}"

# set xtics ( "0.05" 0.05, "0.1" 0.1, "0.15" 0.15, "0.2" 0.2, "0.25" 0.25)
set xtics ( "3" 0.05, "4" 0.1, "5" 0.15, "6" 0.2, "7" 0.25)

set log y

plot ARG1  using 1:3 notitle with linespoints linetype 1 pointtype 10, \
ARG1  using 1:5 notitle with linespoints linetype 1 pointtype 6

#The 5th plot (notitle)
#=========================================
set size 1.400000,0.6000000;  
set origin 1.28,0.55;

set xlabel  "Steiner points m"
set ylabel  "Query Time (ms)"
set key above

set xrange [0.020000: 0.28000000]
set yrange [0.1: 100000]
set label 11 center at graph 0.5,-0.5, char 1 "(c) Query Time" 
set bmargin 5
set format x "%g"
#set format y "10^{%T}"

# set xtics ( "0.05" 0.05, "0.1" 0.1, "0.15" 0.15, "0.2" 0.2, "0.25" 0.25)
set xtics ( "3" 0.05, "4" 0.1, "5" 0.15, "6" 0.2, "7" 0.25)

set log y

plot ARG1  using 1:6 notitle with linespoints linetype 1 pointtype 1, \
ARG1  using 1:7 notitle with linespoints linetype 1 pointtype 4, \
ARG1  using 1:8 notitle with linespoints linetype 1 pointtype 8, \
ARG1  using 1:9 notitle with linespoints linetype 1 pointtype 10, \
ARG1  using 1:10 notitle with linespoints linetype 1 pointtype 6

#The 6th plot (notitle)
#=========================================
set size 1.400000,0.6000000;  
set origin 1.28,0.0;

set xlabel  "Steiner points m"
set ylabel  "Relative Error"
set key above

unset log y
set xrange [0.020000: 0.28000000]
# set yrange [0.000: 0.25]
set yrange [0.000: 0.06001]
set label 11 center at graph 0.5,-0.5, char 1 "(d) Relative Error" 
set bmargin 5
set format x "%g"
#set format y "10^{%T}"
# set xtics ( "0.05" 0.05, "0.1" 0.1, "0.15" 0.15, "0.2" 0.2, "0.25" 0.25)
set xtics ( "3" 0.05, "4" 0.1, "5" 0.15, "6" 0.2, "7" 0.25)
# set ytics ("0.01" 0.0, "0.05" 0.05, "0.10" 0.1, "0.15" 0.15, "0.20" 0.2, "0.25" 0.25)
# set ytics ("0.01" 0.01, "0.02" 0.02, "0.03" 0.03, "0.04" 0.04, "0.05" 0.05, "0.06" 0.06)
set ytics ("0.00" 0.00, "0.02" 0.02, "0.04" 0.04, "0.06" 0.06, "0.08" 0.08, "0.10" 0.10)
# set label 12 'Default Error Bound ({/Symbol e} = 0.2)' center at graph 0.5, 0.7 font "Helvetica, 32"

plot ARG1 using 1:11 notitle with linespoints linetype 1 pointtype 1, \
ARG1 using 1:12 notitle with linespoints linetype 1 pointtype 4, \
ARG1 using 1:13 notitle with linespoints linetype 1 pointtype 8, \
ARG1 using 1:14 notitle with linespoints linetype 1 pointtype 10, \
ARG1 using 1:15 notitle with linespoints linetype 1 pointtype 6, \
0.2 notitle with lines linetype 2

#ARG1."_Averagerr.res" notitle with linespoints linetype 2

unset multiplot