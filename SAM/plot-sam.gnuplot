reset
# svg 480,360
unset key
set border 2 lw 2
set grid
unset xtics
set ytics nomirror
set style fill solid
set ytics (0, 1, 2, 3, 4) rotate
set yrange [-1:5]
set ylabel "Overexpression in SAM relative to actin" offset 1,0
set boxwidth 0.8

set xrange [-0.8:*] 
strmatch(x) = (x eq "msh1" ? 3 : (x eq "reca2" ? 2 : (x eq "osb4" || x eq "reca3" || x eq "odb1" ? 1: 0)))
pval(mean, sem) = 1./(sqrt(1./(sem*sem))) - sem*erf((mean-1.)/(sqrt(2.)*sem))

plot "geneexp-summary.txt" u (strmatch(strcol(3)) == 3 ? $1: 1/0):4:5 w yerrorbars pt 1 lw 2 lc rgbcolor "#0000FF", "" u (strmatch(strcol(3)) == 2 ? $1: 1/0):4:5 w yerrorbars pt 1 lw 2 lc rgbcolor "#FF0000", "" u ((strmatch(strcol(3))) == 1 ? $1: 1/0):4:5 w yerrorbars pt 1 lw 2 lc rgbcolor "#AAAAFF", "" u (strmatch(strcol(3)) == 0 ? $1: 1/0):4:5 w yerrorbars pt 1 lw 2 lc rgbcolor "#FFAAAA",        "" u (strmatch(strcol(3)) == 3 ? $1: 1/0):4 w boxes lc rgbcolor "#0000FF", "" u (strmatch(strcol(3)) == 2 ? $1: 1/0):4 w boxes lc rgbcolor "#FF0000", "" u (strmatch(strcol(3)) == 1 ? $1: 1/0):4 w boxes lc rgbcolor "#AAAAFF", "" u (strmatch(strcol(3)) == 0 ? $1: 1/0):4 w boxes lc rgbcolor "#FFAAAA",         "" u 1:(-0.7):2 w labels rotate tc rgbcolor "#AAAAAA", "" u (strmatch(strcol(3)) % 2 == 1 ? $1: 1/0):(4.3):3 w labels rotate tc rgbcolor "#AAAAFF", "" u (strmatch(strcol(3)) % 2 == 0 ? $1: 1/0):(4.3):3 w labels rotate tc rgbcolor "#FFAAAA", 1 lc rgbcolor "#888888" lw 2     ,     "" u 1:(4.85):($5 == 0 ? "†" : (pval($4, $5) < 0.001/25 ? "***" : (pval($4, $5) < 0.01/25 ? "**" : (pval($4, $5) < 0.05/25 ? "*" : "")))) w labels rotate tc rgbcolor "#000000"
