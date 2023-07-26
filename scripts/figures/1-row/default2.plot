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

set out "out/default.eps"
set terminal postscript portrait enhanced mono "Helvetica" 32 

# set size 3.0000, 1.3200000
set size 5.0000, 0.800000
set pointsize 3.000000

set multiplot layout 2,4
#The first plot (which generates the common legend)
#===========================================================
#Note the NaN variable.

set size 5.000000,0.100000;  
set origin -0.15, 0.62;
# set origin -1.550,0.62;

set key center top horizontal samplen 2
set yrange [0.0000000001:0.0000000002]

unset border
unset tics
#unset label

plot NaN title 'FS(Fixed Scheme)' lt 1 with boxes fs solid 0.25, \
NaN title 'US(Unfixed Scheme)' lt 1 with boxes fs solid 0.6, \
NaN title 'K-Algo' lt 1 with boxes fs pattern 1, \
NaN title 'SE-Oracle' lt 1 with boxes fs pattern 5, \
NaN title 'EAR-Oracle' lt 1 with boxes fs pattern 3#, \
# NaN title 'Theoretical bound ({/Symbol e})' with lines linetype 2
 
set border
set tics
set label

#Plot default building time
#=========================================
set size 1.20000,0.6000000;  
# set origin 0.0,0.6;
set origin 0.0,0.0;

set xlabel  "dataset"

# set xlabel font ", 30"
# set xtics font ", 25"
set autoscale xfix
set ylabel  "Building Time (s)"
# set key left
unset key

set yrange [10 : 100000]
set xrange [0 : 55]
set label 11 center at graph 0.5,-0.5, char 1 "(a) Building Time" 
set bmargin 5
set format x "%g"
set format y "10^{%T}"

# set xtics ("AL" 3, "" 6.5, "HE" 10, "" 13.5, "HO" 17, "" 20.5,  "BH" 24, "" 27.5,  "EP" 31, "" 34.5,  "SF" 38)
set xtics ("HM" 3, "" 6.5, "BM" 10, "" 13.5, "HL" 17, "" 20.5,  "RM" 24, "" 27.5,  "GF" 31, "" 34.5,  "LM" 38, "" 41.5, "BH" 45, "" 48.5, "EP" 52)
# set xtics rotate by 45

set log y

plot ARG1 using ($1-0.75):3:(1.5) title 'SE-Oracle' lt 1 with boxes fs pattern 5, \
ARG1 using ($1+0.75):4:(1.5) title 'EAR-Oracle' lt 1 with boxes fs pattern 3

#Plot default index size
#=========================================
set size 1.20000,0.6000000;  
# set origin 1.5,0.6;
set origin 1.1,0.0;

set xlabel  "dataset"
# set xlabel font ", 30"
# set xtics font ", 25"
set ylabel  "Index Size (MB)"
# set key left
unset key

set yrange [10 : 500000]
# set xrange [3.5 : 19]
set xrange [0 : 55]
set label 11 center at graph 0.5,-0.5,char 1 "(b) Space Consumption" 
set bmargin 5
set format x "%g"
set format y "10^{%T}"

# set xtics ("Alborz" 5, "" 7.5, "Heart" 10, "" 12.5, "Honolu" 15, "" 17.5,  "BH" 20, "" 22.5,  "EP" 25, "" 27.5,  "SF" 30)
# set xtics ("AL" 3, "" 6.5, "HE" 10, "" 13.5, "HO" 17, "" 20.5,  "BH" 24, "" 27.5,  "EP" 31, "" 34.5,  "SF" 38)
# set xtics ("HM" 3, "" 6.5, "BM" 10, "" 13.5, "HL" 17, "" 20.5,  "RM" 24, "" 27.5,  "GF" 31, "" 34.5,  "LM" 38)
set xtics ("HM" 3, "" 6.5, "BM" 10, "" 13.5, "HL" 17, "" 20.5,  "RM" 24, "" 27.5,  "GF" 31, "" 34.5,  "LM" 38, "" 41.5, "BH" 45, "" 48.5, "EP" 52)



set log y

plot ARG1 using ($1-0.75):5:(1.5) title 'SE-Oracle' lt 1 with boxes fs pattern 5, \
ARG1 using ($1+0.75):6:(1.5) title 'EAR-Oracle' lt 1 with boxes fs pattern 3


#Plot default query time
#=========================================
set size 1.50000,0.6000000;  
# set origin 0.0,0.0;
set origin 2.2,0.0;

set xlabel  "dataset"
# set xlabel font ", 30"
# set xtics font ", 25"
set ylabel  "Query Time (ms)"
# set key left
unset key

set yrange [0.5 : 1000000]
set xrange [0 : 80]
set label 11 center at graph 0.5,char 1 "(c) Query Time" 
set bmargin 5
set format x "%g"
set format y "10^{%T}"

# set xtics ("AL" 5, "" 10, "HE" 15, "" 20, "HO" 25, "" 30,  "BH" 35, "" 40,  "EP" 45, "" 50,  "SF" 55)
set xtics ("HM" 5, "" 10, "BM" 15, "" 20, "HL" 25, "" 30,  "RM" 35, "" 40,  "GF" 45, "" 50,  "LM" 55, "" 60, "BH" 65, "" 70, "EP" 75)
set log y

plot ARG1 using ($2 - 3):7:(1.5) title 'FS' lt 1 with boxes fs solid 0.25, \
ARG1 using ($2 - 1.5):8:(1.5) title 'US' lt 1 with boxes fs solid 0.6, \
ARG1 using ($2):9:(1.5) title 'K-Algo' lt 1 with boxes fs pattern 1, \
ARG1 using ($2 + 1.5):10:(1.5) title 'SE-Oracle' lt 1 with boxes fs pattern 5, \
ARG1 using ($2 + 3):11:(1.5) title 'EAR-Oracle' lt 1 with boxes fs pattern 3


#Plot default average error
#=========================================
set size 1.50000,0.6000000;  
# set origin 1.5,0.0;
set origin 3.6,0.0;

set xlabel  "dataset"
# set xlabel font ", 30"
# set xtics font ", 25"
set ylabel  "Relative Error"
# set key left
unset key
unset log y
set yrange [0.000 : 0.06001]
set xrange [0 : 80]
set label 11 center at graph 0.5,char 1 "(d) Relative Error" 
set bmargin 5
set format x "%g"
# set format y "10^{%T}"

# set xtics ("Alborz" 5, "Everest" 15, "Grand" 25, "BH-low" 35, "EP-low" 45, "SF-low" 55)
# set xtics ("AL" 5, "" 10, "HE" 15, "" 20, "HO" 25, "" 30,  "BH" 35, "" 40,  "EP" 45, "" 50,  "SF" 55)
# set xtics ("RM" 5, "" 10, "TM" 15, "" 20, "VR" 25, "" 30,  "SM" 35, "" 40,  "GF" 45, "" 50,  "LM" 55)
# set xtics ("HM" 5, "" 10, "BM" 15, "" 20, "HL" 25, "" 30,  "RM" 35, "" 40,  "GF" 45, "" 50,  "LM" 55)
set xtics ("HM" 5, "" 10, "BM" 15, "" 20, "HL" 25, "" 30,  "RM" 35, "" 40,  "GF" 45, "" 50,  "LM" 55, "" 60, "BH" 65, "" 70, "EP" 75)



# set ytics ("0.01" 0.0, "0.05" 0.05, "0.10" 0.1, "0.15" 0.15, "0.20" 0.2, "0.25" 0.25)
# set ytics ("0.01" 0.01, "0.02" 0.02, "0.03" 0.03, "0.04" 0.04, "0.05" 0.05, "0.06" 0.06)
set ytics ("0.00" 0.00, "0.02" 0.02, "0.04" 0.04, "0.06" 0.06, "0.08" 0.08, "0.10" 0.10)
# set label 12 'Default Error Bound ({/Symbol e} = 0.2)' center at graph 0.5, 0.7 font "Helvetica, 32"
# set log y

plot ARG1 using ($2 - 3):12:(1.5) title 'FS' lt 1 with boxes fs solid 0.25, \
ARG1 using ($2 - 1.5):13:(1.5) title 'US' lt 1 with boxes fs solid 0.6, \
ARG1 using ($2):14:(1.5) title 'K-Algo' lt 1 with boxes fs pattern 1, \
ARG1 using ($2 + 1.5):15:(1.5) title 'SE-Oracle' lt 1 with boxes fs pattern 5, \
ARG1 using ($2 + 3):16:(1.5) title 'EAR-Oracle' lt 1 with boxes fs pattern 3, \
0.2 title 'Default Theoretical Bound' with lines linetype 2

unset multiplot

