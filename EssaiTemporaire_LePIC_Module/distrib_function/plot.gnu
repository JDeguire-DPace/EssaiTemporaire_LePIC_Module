 set style line 1 lt 1 lw 3 # red
 set style line 2 lt 2 lw 3 # green
 set style line 3 lt 3 lw 3 # blue
 set style line 4 lt 4 lw 3 # purple
 set style line 5 lt 5 lw 3 # light blue
 set style line 6 lt 6 lw 3 # yellow
 set style line 7 lt 7 lw 3 # black
 set style line 8 lt 7 lw 5 # thick black
  
 set size square
 set terminal postscript eps solid color enhanced "Helvetica" 24
 set output "fig1.eps"
  
 set xrange [-4:4]
 set xlabel "v_i/v_t"
 plot "pdf.dat" u 1:2 w histeps ls 1 noti,"pdf.dat" u 1:3 w histeps ls 2 noti, "pdf.dat" u 1:4 w histeps ls 3 noti,exp(-x**2.)/sqrt(pi) ls 7 noti
  
 x0=  8.358193365397583E-002
 xmax=   49.9916454241071     
 ymin=  1.218638706820055E-006
 Tp=   4.50635887066550     
 set logscale y
 set format y "%2.0t{/Symbol \327}10^{%L}"
 set xrange[0:xmax]
 set yrange[ymin:*]
 set xlabel "E_k(eV)"
 set ylabel "f (eV^{-3/2})"
 plot "pedf.dat" u 1:2 w l ls 1 noti, x0*exp(-x/Tp) ls 7 noti
