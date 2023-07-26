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

set out "out/default-2row.eps"
set terminal postscript portrait enhanced mono "Helvetica" 32 

# set size 3.0000, 1.3200000
set size 2.50000, 1.2400000
set pointsize 3.000000

set multiplot layout 3,2
#The first plot (which generates the common legend)
#===========================================================
#Note the NaN variable.

set size 2.400000,0.100000;  
set origin -0.15, 1.18;
# set origin -1.550,0.62;

set key center top horizontal samplen 2 font "Helvatica, 28"
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
set origin 0.0,0.55;

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
set origin 0.0,0.0;

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
set origin 1.1,0.55;

set xlabel  "dataset"
# set xlabel font ", 30"
# set xtics font ", 25"
set ylabel  "Query Time (ms)"
# set key left
unset key

set yrange [0.5 : 1000000]
set xrange [0 : 80]
set label 11 center at graph 0.5,-0.5, char 1 "(c) Query Time" 
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
set origin 1.1,0.0;

set xlabel  "dataset"
# set xlabel font ", 30"
# set xtics font ", 25"
set ylabel  "Relative Error"
# set key left
unset key
unset log y
set yrange [0.000 : 0.06001]
set xrange [0 : 80]
set label 11 center at graph 0.5,-0.5, char 1 "(d) Relative Error" 
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


# set out "IT.eps"

# set terminal postscript portrait enhanced mono "Helvetica" 22

# set multiplot layout 2, 2
# #Common legend
# set size 4.000000,0.100000;  
# set origin -1.550,0.62;

# set key center top horizontal samplen 2
# set yrange [0.0000000001:0.0000000002]

# unset border
# unset tics
# #unset label

# plot NaN title 'fixedS' lt 1 with boxes fs solid 0.25, \
# NaN title 'unfixedS' lt 1 with boxes fs solid 0.6, \
# NaN title 'K-Algo' lt 1 with boxes fs pattern 4, \
# NaN title 'SE-Oracle' lt 1 with boxes fs pattern 3, \
# NaN title 'EAR-Oracle' lt 1 with boxes fs pattern 5

# set size 1.100000, 0.60000
# # set size 4.00000, 0.690000
# set origin 0.0,0.0
# set ylabel "Building Time(s)"
# set pointsize 2
# set format y "10^{%L}"

# set key right top
# set xtics font ",20"
# set ytics font ",20"
# set xlabel font ",20"
# set ylabel font ",20"
# set ylabel offset -1
# set key font ",20"
# set key top horiz center
# #set key margin 12
# set key width 4

# set ytics 10
# set xtics 1.00
# set yrange [10 : 500000]
# set xrange [3.5 : 19]
# set xtics ("Alborz" 5, "Everest" 7.5, "Grand" 10, "BH-low" 12.5, "EP-low" 15, "SF-low" 17.5)
# set label 11 center at graph 0.5,char 1 "(a)"
# # set xrange [-5 : 91]
# # set xtics ("HP" -1, "" 3.25,  "CA" 7, "" 11.25, "AD" 15, "" 19.25, "WV" 23, "" 27.25, "WS" 31, "" 35.25, "WG" 39, "" 43.25, "DB" 47, "" 51.25, "UP" 55, "" 59.25, "WP" 63, "" 67.25, "LJ" 71, "" 75.25, "WD" 79, "" 83.25, "WF" 87)

# set logscale y
# set style fill solid border -1


# plot "indexTime.res" using ($1-0.5):5:(1) title 'SE-Oracle' lt 1 with boxes fs pattern 3, \
# "indexTime.res" using ($1+0.5):6:(1) title 'EAR-Oracle' lt 1 with boxes fs pattern 5

# # plot "indexTime.res" using ($1 - 5):2:(1) title 'fixedS' lt 1 with boxes fs solid 0.25, \
# # "indexTime.res" using ($1 - 2.5):3:(1) title 'unfixedS' lt 1 with boxes fs solid 0.6, \
# # "indexTime.res" using ($1):4:(1) title 'K-Algo' lt 1 with boxes fs pattern 4, \
# # "indexTime.res" using ($1 + 2.5):5:(1) title 'SE-Oracle' lt 1 with boxes fs pattern 3, \
# # "indexTime.res" using ($1 + 5):6:(1) title 'EAR-Oracle' lt 1 with boxes fs pattern 5

# # plot "QT.dat" using ($1*8-11.5):3:(1) title 'READS-D' lt 1 with boxes fs solid 0.25,  \
# # 	 "QT.dat" using ($1*8+1-11.5):4:(1) title 'READS-Rq' lt 1 with boxes fs solid 0.6, \
# # 	 "QT.dat" using ($1*8+2-11.5):7:(1) title 'L-TSF' lt 1  with  boxes fs pattern 4, \
# # 	 "QT.dat" using ($1*8+3-11):2:(1) title 'READS' lt 1 with boxes fs pattern 3, \
# # 	 "QT.dat" using ($1*8+4-11):6:(1) title 'TSF' lt 1 with boxes fs pattern 5, \
# # 	 "QT.dat" using ($1*8+5-11):5:(1) title 'SLING' lt 1 with boxes fs pattern 2

# unset multiplot
