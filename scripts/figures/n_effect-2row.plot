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

set out "out/scalability-2row.eps"
set terminal postscript portrait enhanced mono "Helvetica" 32 

set size 2.58000, 1.240000
set pointsize 3.000000


set multiplot layout 3,2
#The first plot (which generates the common legend)
#===========================================================
#Note the NaN variable.

set size 2.50000,0.100000;  
set origin -0.15,1.18;

set key center top horizontal samplen 2 font "Helvatica, 28"
set yrange [0.0000000001:0.0000000002]

unset border
unset tics
#unset label

plot NaN title 'FS(Fixed Scheme)' with linespoints linetype 1 pointtype 1, \
NaN title 'US(Unfixed Scheme)' with linespoints linetype 1 pointtype 4, \
NaN title 'Kalgo' with linespoints linetype 1 pointtype 8, \
NaN title 'SE-Oracle' with linespoints linetype 1 pointtype 10, \
NaN title 'EAR-Oracle' with linespoints linetype 1 pointtype 6#, \
#NaN title 'Theoretical bound ({/Symbol e})' with lines linetype 2
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
set origin 0.0,0.550;

set xlabel  "{data size}"
set ylabel  "Building Time (s)"
set key left

set xrange [0.020000: 0.28000000]
set yrange [1000: 100000]
set label 11 center at graph 0.5,-0.5, char 1 "(a) Building Time" 
set bmargin 5
set format x "%g"
set format y "10^{%T}"

# set xtics ( "1/400" 0.05, "1/225" 0.1, "1/100" 0.15, "1/25" 0.2, "1" 0.25)
set xtics ( "0.2M" 0.05, "0.4M" 0.1, "0.6M" 0.15, "0.8M" 0.2, "1.0M" 0.25)
# set xtics rotate by -90

set log y

plot ARG1."_indexTime.res" using 1:5  notitle with linespoints linetype 1 pointtype 10, \
ARG1."_indexTime.res" using 1:4 notitle with linespoints linetype 1 pointtype 6

#The 4th plot (notitle)
#=========================================
set size 1.400000,0.600000;  
set origin 0.0,0.0;

set xlabel  "{data size}"
set ylabel  "Size (MB)"
set key above

set xrange [0.020000: 0.28000000]
set yrange [10000: 1000000]
set label 11 center at graph 0.5,-0.5, char 1 "(b) Space Consumption" 
set bmargin 5
set format x "%g"
#set format y "10^{%T}"

# set xtics ( "0.05" 0.05, "0.1" 0.1, "0.15" 0.15, "0.2" 0.2, "0.25" 0.25)
# set xtics ( "1/400" 0.05, "1/225" 0.1, "1/100" 0.15, "1/25" 0.2, "1" 0.25)
set xtics ( "0.2M" 0.05, "0.4M" 0.1, "0.6M" 0.15, "0.8M" 0.2, "1.0M" 0.25)

set log y

plot ARG1."_indexSize.res"  using 1:5 notitle with linespoints linetype 1 pointtype 10, \
ARG1."_indexSize.res"  using 1:4 notitle with linespoints linetype 1 pointtype 6

#The 5th plot (notitle)
#=========================================
set size 1.400000,0.6000000;  
set origin 1.28,0.55;

set xlabel  "{data size}"
set ylabel  "Query Time (ms)"
set key above

set xrange [0.020000: 0.28000000]
set yrange [10: 10000000]
set label 11 center at graph 0.5,-0.5, char 1 "(c) Query Time" 
set bmargin 5
set format x "%g"
#set format y "10^{%T}"

# set xtics ( "0.05" 0.05, "0.1" 0.1, "0.15" 0.15, "0.2" 0.2, "0.25" 0.25)
# set xtics ( "1/400" 0.05, "1/225" 0.1, "1/100" 0.15, "1/25" 0.2, "1" 0.25)
set xtics ( "0.2M" 0.05, "0.4M" 0.1, "0.6M" 0.15, "0.8M" 0.2, "1.0M" 0.25)

set log y

plot ARG1."_queryTime.res"  using 1:2 notitle with linespoints linetype 1 pointtype 1, \
ARG1."_queryTime.res"  using 1:6 notitle with linespoints linetype 1 pointtype 4, \
ARG1."_queryTime.res"  using 1:3 notitle with linespoints linetype 1 pointtype 8, \
ARG1."_queryTime.res"  using 1:4 notitle with linespoints linetype 1 pointtype 6, \
ARG1."_queryTime.res"  using 1:5 notitle with linespoints linetype 1 pointtype 10

#The 6th plot (notitle)
#=========================================
set size 1.400000,0.6000000;  
set origin 1.28,0.0;

set xlabel  "{data size}"
set ylabel  "Relative Error"
set key above

unset log y
set xrange [0.020000: 0.28000000]
# set yrange [0.000: 0.25]
set yrange [0.000: 0.03001]
set label 11 center at graph 0.5,-0.5, char 1 "(d) Relative Error" 
set bmargin 5
set format x "%g"
#set format y "10^{%T}"
# set xtics ( "0.05" 0.05, "0.1" 0.1, "0.15" 0.15, "0.2" 0.2, "0.25" 0.25)
# set xtics ( "1/400" 0.05, "1/225" 0.1, "1/100" 0.15, "1/25" 0.2, "1" 0.25)
set xtics ( "0.2M" 0.05, "0.4M" 0.1, "0.6M" 0.15, "0.8M" 0.2, "1.0M" 0.25)
# set ytics ("0.01" 0.0, "0.05" 0.05, "0.10" 0.1, "0.15" 0.15, "0.20" 0.2, "0.25" 0.25)
# set ytics ("0.01" 0.01, "0.02" 0.02, "0.03" 0.03, "0.04" 0.04, "0.05" 0.05, "0.06" 0.06)
set ytics ("0.00" 0.00, "0.01" 0.01, "0.02" 0.02, "0.03" 0.03, "0.04" 0.04, "0.05" 0.05)

# set label 12 'Default Error Bound ({/Symbol e} = 0.2)' center at graph 0.5, 0.7 font "Helvetica, 32"

plot ARG1."_Averagerr.res"  using 1:2 notitle with linespoints linetype 1 pointtype 1, \
ARG1."_Averagerr.res"  using 1:6 notitle with linespoints linetype 1 pointtype 4, \
ARG1."_Averagerr.res"  using 1:3 notitle with linespoints linetype 1 pointtype 8, \
ARG1."_Averagerr.res"  using 1:4 notitle with linespoints linetype 1 pointtype 6, \
ARG1."_Averagerr.res"  using 1:5 notitle with linespoints linetype 1 pointtype 10, \
0.2 notitle with lines linetype 2

#ARG1."_Averagerr.res" notitle with linespoints linetype 2

unset multiplot