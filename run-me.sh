# model fitting to mouse data
cd Mouse
chmod +x run-me.sh
./run-me.sh &
cd ..

# gene expression patterns in plant shoot apical meristems
cd SAM
chmod +x run-me.sh
./run-me.sh &
cd ..

# presence/absence of oDNA recombination genes across taxa
cd Phylogenetics
chmod +x run-me.sh
./run-me.sh &
cd ..

# genetic implications of physical oDNA partitioning at cell divisions
cd Partitioning
chmod +x run-me.sh
./run-me.sh &
cd ..

# stochastic simulations testing theory
cd Stochastics
chmod +x run-me.sh
./run-me.sh &
cd ..
