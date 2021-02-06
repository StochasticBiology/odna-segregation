# gnuplot scripts for all gnuplottable results figures
# only run this after the simulations and analysis (run-me.sh) has been completed!
# the phylogeny (Fig 2) is produced by Mathematica -- see the Phylogenetics/ directory

cd Development
gnuplot -e "set term svg size 640,320; set output \"../plot-bottlenecks.svg\"; load \"plot-bottlenecks.gnuplot\"; quit;"
cd ..

cd Mouse
gnuplot -e "set term svg size 480,320; set output \"../plot-hb.svg\"; load \"plot-hb.gnuplot\"; quit;"
gnuplot -e "set term svg size 480,320; set output \"../plot-le.svg\"; load \"plot-le.gnuplot\"; quit;"
gnuplot -e "set term svg size 1024,800; set output \"../plot-model-comparison.svg\"; load \"plot-model-comparison.gnuplot\"; quit;"
cd ..

cd Partitioning
gnuplot -e "set term svg size 640,640; set output \"../plot-partitioning-simulation.svg\"; load \"plot-partitioning-simulation.gnuplot\"; quit;"
cd ..

cd SAM
gnuplot -e "set term svg size 480,360; set output \"../plot-sam.svg\"; load \"plot-sam.gnuplot\"; quit;"
cd ..

cd Stochastics
gnuplot -e "set term svg size 600,360; set output \"../plot-stochastic-results.svg\"; load \"plot-stochastic-results.gnuplot\"; quit;"
gnuplot -e "set term svg size 1600,400; set output \"../plot-simple-output.svg\"; load \"plot-simple-output.gnuplot\"; quit;"
cd ..
