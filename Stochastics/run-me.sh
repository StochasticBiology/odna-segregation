# compile and run stochastic simulations
gcc -o3 stochastic-simulations.c -lm -o stochastic-simulations.ce
gcc -o3 stochastic-simulations-simple.c -lm -o stochastic-simulations-simple.ce

./stochastic-simulations.ce > ss-out.tmp
./stochastic-simulations-simple.ce > sss-out.tmp

# when complete, use gnuplot scripts for visualisation

