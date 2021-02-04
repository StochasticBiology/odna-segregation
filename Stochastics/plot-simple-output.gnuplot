reset

set multiplot

set xrange [0:10]
set yrange [-0.015:0.15]

# svg 1600,400

nu = 0.5
kappa = 0.001
nbig = 1000.
nsmall(x) = (x < 5 ? 500. : (x < 10 ? 800 : 200))
nc(x) = (x == 1 ? 1. : (x == 2 ? 2. : (x == 3 ? 4. : (x == 4 ? 16. : 1/0))))

# redo for subsampling
#f(a, b, c, d, x) = x*( (c == 1 ? (b%5 == 1 ? 1./nbig : (b%5 == 2 ? 2./nbig : (b%5 == 3 ? 4./nbig : (b%5 == 4 ? 16./nbig : 0)))) : 0) + (c == 2 ? 2.*(b%5 == 1 ? 1./nbig : (b%5 == 2 ? 2./nbig : (b%5 == 3 ? 4./nbig : (b%5 == 4 ? 16./nbig : 0)))) : 0) + (d == 1 ? (1./nsmall(b) - 1./nbig) : 0) + (a == 1 || a == 3 || a == 4 || a == 6 ? 2.*nu/nbig : 0) + (a == 2 || a == 3 || a == 5 || a == 6 ? 2*kappa : 0) )

# c == 1: hypergeometric sampling
# c == 2: binomial sampling
f(a, b, c, d, x) = x*( (c == 1 ? nc(b%5)/nsmall(b) * (nbig - nsmall(b))/(nbig) : 0) + (c == 2 ? nc(b%5)/nsmall(b) : 0) + (d == 1 ? (1./nsmall(b) - 1./nbig) : 0) + (a == 1 || a == 3 || a == 4 || a == 6 ? 2.*nu/nbig : 0) + (a == 2 || a == 3 || a == 5 || a == 6 ? 2*kappa : 0) )

set key outside left

# gillespie-divisions-rsample-ramplify

# divisions/5: 0 -- 500; 1 -- 800; 2 -- 200

set xlabel "Time"
set ylabel "V'(h)"

set size 0.33,1

set pointsize 0.5

set label 1 at 0.5,-0.0075 "n1 = 1000, n2 = 500"
set origin 0,0
plot "simple-output-0-0-0-0.txt" u 1:4 ls 1 title "DD/DA", f(0,0,0,0,x) ls 1 notitle, "simple-output-0-1-0-1.txt" u 1:4 ls 2 title "DD/RA", f(0,1,0,1,x) ls 2 notitle, "simple-output-0-1-1-0.txt" u 1:4 ls 3 title "HD/DA", f(0,1,1,0,x) ls 3 notitle, "simple-output-0-1-1-1.txt" u 1:4 ls 4 title "HD/RA", f(0,1,1,1,x) ls 4 notitle, "simple-output-0-1-2-0.txt" u 1:4 ls 5 title "BD/DA", f(0,1,2,0,x) ls 5 notitle, "simple-output-0-1-2-1.txt" u 1:4 ls 6 title "BD/RA", f(0,1,2,1,x) ls 6 notitle, "simple-output-0-2-1-1.txt" u 1:4 ls 7 title "HD2/RA", f(0,2,1,1,x) ls 7 notitle, "simple-output-0-3-1-1.txt" u 1:4 ls 8 title "HD4/RA", f(0,3,1,1,x) ls 8 notitle, "simple-output-0-4-1-1.txt" u 1:4 ls 9 title "HD16/RA", f(0,4,1,1,x) ls 9 notitle, "simple-output-0-2-1-0.txt" u 1:4 ls 10 title "HD2/DA", f(0,2,1,0,x) ls 10 notitle, "simple-output-0-3-1-0.txt" u 1:4 ls 11 title "HD4/DA", f(0,3,1,0,x) ls 11 notitle, "simple-output-0-4-1-0.txt" u 1:4 ls 12 title "HD16/DA", f(0,4,1,0,x) ls 12 notitle,    "simple-output-0-2-2-1.txt" u 1:4 ls 13 title "B2/RA", f(0,2,2,1,x) ls 13 notitle, "simple-output-0-3-2-1.txt" u 1:4 ls 14 title "B4/RA", f(0,3,2,1,x) ls 14 notitle, "simple-output-0-4-2-1.txt" u 1:4 ls 15 title "B16/RA", f(0,4,2,1,x) ls 15 notitle, "simple-output-0-2-2-0.txt" u 1:4 ls 16 title "B2/DA", f(0,2,2,0,x) ls 16 notitle, "simple-output-0-3-2-0.txt" u 1:4 ls 17 title "B4/DA", f(0,3,2,0,x) ls 17 notitle, "simple-output-0-4-2-0.txt" u 1:4 ls 18 title "B16/DA", f(0,4,2,0,x) ls 18 notitle 

set label 1 at 0.5,-0.0075 "n1 = 1000, n2 = 800"
set origin 0.33,0
plot "simple-output-0-5-0-0.txt" u 1:4 ls 1 title "DD/DA", f(0,5,0,0,x) ls 1 notitle, "simple-output-0-6-0-1.txt" u 1:4 ls 2 title "DD/RA", f(0,6,0,1,x) ls 2 notitle, "simple-output-0-6-1-0.txt" u 1:4 ls 3 title "HD/DA", f(0,6,1,0,x) ls 3 notitle, "simple-output-0-6-1-1.txt" u 1:4 ls 4 title "HD/RA", f(0,6,1,1,x) ls 4 notitle, "simple-output-0-6-2-0.txt" u 1:4 ls 5 title "BD/DA", f(0,6,2,0,x) ls 5 notitle, "simple-output-0-6-2-1.txt" u 1:4 ls 6 title "BD/RA", f(0,6,2,1,x) ls 6 notitle, "simple-output-0-7-1-1.txt" u 1:4 ls 7 title "HD2/RA", f(0,7,1,1,x) ls 7 notitle, "simple-output-0-8-1-1.txt" u 1:4 ls 8 title "HD4/RA", f(0,8,1,1,x) ls 8 notitle, "simple-output-0-9-1-1.txt" u 1:4 ls 9 title "HD16/RA", f(0,9,1,1,x) ls 9 notitle, "simple-output-0-7-1-0.txt" u 1:4 ls 10 title "HD2/DA", f(0,7,1,0,x) ls 10 notitle, "simple-output-0-8-1-0.txt" u 1:4 ls 11 title "HD4/DA", f(0,8,1,0,x) ls 11 notitle, "simple-output-0-9-1-0.txt" u 1:4 ls 12 title "HD16/DA", f(0,9,1,0,x) ls 12 notitle,    "simple-output-0-7-2-1.txt" u 1:4 ls 13 title "B2/RA", f(0,7,2,1,x) ls 13 notitle, "simple-output-0-8-2-1.txt" u 1:4 ls 14 title "B4/RA", f(0,8,2,1,x) ls 14 notitle, "simple-output-0-9-2-1.txt" u 1:4 ls 15 title "B16/RA", f(0,9,2,1,x) ls 15 notitle, "simple-output-0-7-2-0.txt" u 1:4 ls 16 title "B2/DA", f(0,7,2,0,x) ls 16 notitle, "simple-output-0-8-2-0.txt" u 1:4 ls 17 title "B4/DA", f(0,8,2,0,x) ls 17 notitle, "simple-output-0-9-2-0.txt" u 1:4 ls 18 title "B16/DA", f(0,9,2,0,x) ls 18 notitle 

set label 1 at 0.5,-0.0075 "n1 = 1000, n2 = 200"
set origin 0.66,0
plot "simple-output-0-10-0-0.txt" u 1:4 ls 1 title "DD/DA", f(0,10,0,0,x) ls 1 notitle, "simple-output-0-11-0-1.txt" u 1:4 ls 2 title "DD/RA", f(0,11,0,1,x) ls 2 notitle, "simple-output-0-11-1-0.txt" u 1:4 ls 3 title "HD/DA", f(0,11,1,0,x) ls 3 notitle, "simple-output-0-11-1-1.txt" u 1:4 ls 4 title "HD/RA", f(0,11,1,1,x) ls 4 notitle, "simple-output-0-11-2-0.txt" u 1:4 ls 5 title "BD/DA", f(0,11,2,0,x) ls 5 notitle, "simple-output-0-11-2-1.txt" u 1:4 ls 6 title "BD/RA", f(0,11,2,1,x) ls 6 notitle, "simple-output-0-12-1-1.txt" u 1:4 ls 7 title "HD2/RA", f(0,12,1,1,x) ls 7 notitle, "simple-output-0-13-1-1.txt" u 1:4 ls 8 title "HD4/RA", f(0,13,1,1,x) ls 8 notitle, "simple-output-0-14-1-1.txt" u 1:4 ls 9 title "HD16/RA", f(0,14,1,1,x) ls 9 notitle, "simple-output-0-12-1-0.txt" u 1:4 ls 10 title "HD2/DA", f(0,12,1,0,x) ls 10 notitle, "simple-output-0-13-1-0.txt" u 1:4 ls 11 title "HD4/DA", f(0,13,1,0,x) ls 11 notitle, "simple-output-0-14-1-0.txt" u 1:4 ls 12 title "HD16/DA", f(0,14,1,0,x) ls 12 notitle,    "simple-output-0-12-2-1.txt" u 1:4 ls 13 title "B2/RA", f(0,12,2,1,x) ls 13 notitle, "simple-output-0-13-2-1.txt" u 1:4 ls 14 title "B4/RA", f(0,13,2,1,x) ls 14 notitle, "simple-output-0-14-2-1.txt" u 1:4 ls 15 title "B16/RA", f(0,14,2,1,x) ls 15 notitle, "simple-output-0-12-2-0.txt" u 1:4 ls 16 title "B2/DA", f(0,12,2,0,x) ls 16 notitle, "simple-output-0-13-2-0.txt" u 1:4 ls 17 title "B4/DA", f(0,13,2,0,x) ls 17 notitle, "simple-output-0-14-2-0.txt" u 1:4 ls 18 title "B16/DA", f(0,14,2,0,x) ls 18 notitle 

